/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "smoothMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "DynamicList.H"
#include "cellQuality.H"
#include "fvCFD.H"
#include "solverPerformance.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(smoothMatrix, 0);
    addToRunTimeSelectionTable
    (
        option,
        smoothMatrix,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::smoothMatrix::applyNeighbourAveraging(fvMatrix<vector>& eqn)
{
    const volVectorField& U = eqn.psi();

    // Automatically detect bad cells based on mesh quality
    DynamicList<label> badCells;
    const scalarField& volumes = mesh_.V();

    // Get mesh quality metrics using OpenFOAM utilities
    //Can include more metrics if needed. Maybe residual based metrics?
    cellQuality meshQuality(mesh_);
    scalarField nonOrthogonality;
    scalarField skewness;

    if (nonOrthTol_ > 0)
    {
        nonOrthogonality = meshQuality.nonOrthogonality();
    }

    if (skewTol_ > 0)
    {
        skewness = meshQuality.skewness();
    }

    // Exclude boundary cells from consideration
    boolList isBoundaryCell(mesh_.nCells(), false);
    const fvBoundaryMesh& bMesh = mesh_.boundary();
    forAll(bMesh, patchi)
    {
        const labelUList& fc = bMesh[patchi].faceCells();
        forAll(fc, i)
        {
            isBoundaryCell[fc[i]] = true;
        }
    }

    const bool withinAvgWindow =
    (
        residualTol_ > 0
     && iterationCounter_ >= avgStartIter_
     && iterationCounter_ <= avgEndIter_
    );

    // Optional: compute linear-system residual vector (per cell) once,
    // track maximum over iterations if threshold requested
    tmp<Field<vector>> tResidual;
    const Field<vector>* residualVecPtr = nullptr;
    if (withinAvgWindow)
    {
        tResidual = eqn.residual();
        residualVecPtr = &tResidual();

        if (!residualSumFieldPtr_.valid())
        {
            residualSumFieldPtr_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        name_ + "_residualSum",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("zero", dimless, 0.0)
                )
            );
        }
        // One more sample in window this iteration
        ++windowSamples_;
    }

    forAll(mesh_.cells(), celli)
    {
        // Skip boundary cells
        if (excludeBoundary_ && isBoundaryCell[celli])
        {
            continue;
        }

        bool isBad = false;

        // Accumulate per-cell residual magnitude within window
        if (residualVecPtr && residualSumFieldPtr_.valid())
        {
            scalar& rsum = residualSumFieldPtr_().internalFieldRef()[celli];
            const scalar rmag = mag((*residualVecPtr)[celli]);
            rsum += rmag;
        }

        // Mesh quality criteria: always active during the simulation
        if ((smallVolTol_ > 0 && volumes[celli] < smallVolTol_)
         || (skewTol_ > 0 && skewness[celli] > skewTol_)
         || (nonOrthTol_ > 0 && nonOrthogonality[celli] > nonOrthTol_))
        {
            isBad = true;
        }

        // After averaging window ends: apply residual-based post-filter mask
        if (!isBad && residualPostApplied_ && residualTol_ > 0)
        {
            if (celli < residualBadCells_.size() && residualBadCells_[celli])
            {
                isBad = true;
            }
        }

        if (isBad)
        {
            badCells.append(celli);
        }
    }

    if (badCells.empty())
    {
        return;
    }

    labelList badCellsList(badCells);
    Field<vector> values(badCellsList.size(), vector::zero);

    const labelListList& cellCells = mesh_.cellCells();
    const vectorField& cellCentres = mesh_.cellCentres();

    const scalarField* rhoCellsPtr = nullptr;
    if (useDensityWeighting_)
    {
        if (mesh_.foundObject<volScalarField>(rhoName_))
        {
            const volScalarField& rhoFld = mesh_.lookupObject<volScalarField>(rhoName_);
            rhoCellsPtr = &rhoFld.internalField();
        }
        else
        {
            Info<< type() << ' ' << name_ << " density field '" << rhoName_ << "' not found. Disabling density weighting." << endl;
        }
    }

    forAll(badCellsList, i)
    {
        const label celli = badCellsList[i];
        const labelList& neighbours = cellCells[celli];

        vector weightedSum(vector::zero);
        scalar totalWeight = 0.0;

        for (const label n : neighbours)
        {
            if (n >= 0 && n < mesh_.nCells() && n != celli)
            {
                scalar weight = 1.0;  // Default uniform weight

                if (weightingMethod_ == "inverseDistance")
                {
                    vector delta = cellCentres[n] - cellCentres[celli];
                    scalar distance = mag(delta);
                    weight = 1.0 / max(distance, VSMALL);
                }
                else if (weightingMethod_ == "inverseDistanceSquared")
                {
                    vector delta = cellCentres[n] - cellCentres[celli];
                    scalar distance = mag(delta);
                    weight = 1.0 / max(sqr(distance), VSMALL);
                }
                else if (weightingMethod_ == "uniform")
                {
                    weight = 1.0;
                }
                else if (weightingMethod_ == "volumeWeighted")
                {
                    weight = max(volumes[n], VSMALL);
                }
                // Add more weighting methods here as needed

                if (rhoCellsPtr)
                {
                    const scalar rho_n = (*rhoCellsPtr)[n];
                    weightedSum += weight * (rho_n * U[n]);
                    totalWeight += weight* rho_n;
                    //Info << "rho_n: " << rho_n << " weight: " << weight << " totalWeight: " << totalWeight << endl;

                }
                else
                {
                    weightedSum += weight * U[n];
                    totalWeight += weight;
                }
            }
        }

        if (totalWeight > VSMALL)
        {
            //if (rhoCellsPtr)
           // {Info << "rhoWeighted average"<< endl;}

            values[i] = weightedSum / totalWeight;
        }
        else
        {
            values[i] = U[celli];
        }

        // Mark bad cells in the field if required
        if (writeBadCells_ && badCellFieldPtr_)
        {
            badCellFieldPtr_().internalFieldRef()[celli] = 1.0;
        }
    }

    // Enforce values for the selected cells directly in the matrix
    if (constraintMethod_ == "softBlend")
    {
        // Use setReference method - sets reference values with relaxation
            eqn.setReferences(badCellsList, values, true);

        Info<< type() << ' ' << name_ << " soft blended values for "
            << badCellsList.size() << " bad cells from neighbour average" << endl;
    }
    else // Default to setValues
    {
        // Use setValues method - directly enforces values
        eqn.setValues(badCellsList, values);

        Info<< type() << ' ' << name_ << " hard blended for "
            << badCellsList.size() << " bad cells from neighbour average" << endl;
    }

    // Finalize residual-based mask once when the window has ended
    if (residualTol_ > 0 && avgEndIter_ > 0 && iterationCounter_ > avgEndIter_ && !residualPostApplied_)
    {
        label nMarked = 0;
        if (residualSumFieldPtr_.valid() && windowSamples_ > 0)
        {
            const scalarField& rsum = residualSumFieldPtr_().internalField();
            residualBadCells_.setSize(mesh_.nCells(), false);
            forAll(rsum, celli)
            {
                const scalar ravg = rsum[celli]/windowSamples_;
                const bool markBad = (ravg > residualTol_);
                residualBadCells_[celli] = markBad;
                if (markBad) ++nMarked;
            }
        }
        residualPostApplied_ = true;
        Info<< type() << ' ' << name_ << " residual post-filter applied after iteration "
            << avgEndIter_ << ". Marked " << nMarked << " cells (threshold "
            << residualTol_ << ", samples " << windowSamples_ << ')' << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::smoothMatrix::smoothMatrix
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    smallVolTol_(0),
    skewTol_(0),
    nonOrthTol_(0),
    residualTol_(0),
    avgStartIter_(1),
    avgEndIter_(0),
    residualPostApplied_(false),
    iterationCounter_(0),
    windowSamples_(0),
    residualSumFieldPtr_(nullptr),
    residualBadCells_(mesh.nCells(), false),
    weightingMethod_("inverseDistance"),
    constraintMethod_("hardBlend"),
    useDensityWeighting_(false),
    rhoName_("rho"),
    excludeBoundary_(false),
    writeBadCells_(false),
    badCellFieldPtr_(nullptr),
    writeResiduals_(false),
    residualFieldPtr_(nullptr)
{
    read(dict);

    if (writeBadCells_)
    {
        badCellFieldPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    name_ + "_badCells",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::smoothMatrix::constrain
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    // Count constrain() calls as iterations
    ++iterationCounter_;

    // Apply neighbour averaging if any thresholds are configured
    if (smallVolTol_ > 0 || skewTol_ > 0 || nonOrthTol_ > 0 || residualTol_ > 0)
    {
        applyNeighbourAveraging(eqn);
    }

}





void Foam::fv::smoothMatrix::writeData(Ostream& os) const
{
    dict_.writeEntry(name_, os);
}


bool Foam::fv::smoothMatrix::read(const dictionary& dict)
{
    if (!fv::cellSetOption::read(dict))
    {
        return false;
    }

    // Read thresholds for automatic bad cell detection (all optional)
    dict.readIfPresent("smallVolTol", smallVolTol_);
    dict.readIfPresent("skewTol", skewTol_);
    dict.readIfPresent("nonOrthTol", nonOrthTol_);
    // Residual threshold (optional). Activated only when provided (> 0)
    if (!dict.readIfPresent("maxResidual", residualTol_))
    {
        dict.readIfPresent("residualTol", residualTol_);
    }

    // Read weighting method (optional)
    dict.readIfPresent("weightingMethod", weightingMethod_);
    // Read constraint method (optional)
    dict.readIfPresent("constraintMethod", constraintMethod_);

    // Averaging iteration window (optional, but required if residualTol_ > 0)
    bool haveStart = false, haveEnd = false;
    haveStart = dict.readIfPresent("avgStartIter", avgStartIter_);
    haveEnd   = dict.readIfPresent("avgEndIter", avgEndIter_);
    // Alternative keys
    bool altStart = dict.readIfPresent("averageStartIter", avgStartIter_);
    bool altEnd   = dict.readIfPresent("averageEndIter", avgEndIter_);
    haveStart = haveStart || altStart;
    haveEnd   = haveEnd   || altEnd;

    if (residualTol_ > 0)
    {
        if (!haveStart || !haveEnd)
        {
            FatalIOErrorInFunction(dict)
                << "Option 'maxResidual' specified but averaging window not provided. "
                << "Please set 'avgStartIter' and 'avgEndIter' (or 'averageStartIter'/'averageEndIter')."
                << exit(FatalIOError);
        }
        if (avgEndIter_ < avgStartIter_)
        {
            FatalIOErrorInFunction(dict)
                << "'avgEndIter' (" << avgEndIter_ << ") is less than 'avgStartIter' (" << avgStartIter_ << ")."
                << exit(FatalIOError);
        }
    }

    // Density weighting options (optional)
    dict.readIfPresent("useDensityWeighting", useDensityWeighting_);
    dict.readIfPresent("rhoName", rhoName_);
    dict.readIfPresent("excludeBoundary", excludeBoundary_);
    if (!dict.readIfPresent("UNames", fieldNames_))
    {
        fieldNames_.resize(1);
        fieldNames_.first() = coeffs_.getOrDefault<word>("U", "U");
    }

    // Write bad cell list to file (optional)
    dict.readIfPresent("writeBadCells", writeBadCells_);
    fv::option::resetApplied();

    return true;
}


// ************************************************************************* //
