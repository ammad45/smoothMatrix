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

Class
    Foam::fv::smoothMatrix

Group
    grpFvOptionsConstraints

Description
    Smooth velocity fields in cells with poor mesh quality by replacing
    values with neighbour-averaged values.

    This constraint automatically detects cells with poor mesh quality
    (small volume, high skewness, high non-orthogonality) and optionally
    high linear-system residuals, and replaces
    their velocity values with the average of their neighbouring cells.

    Constraints applied to:
    \verbatim
      U    | Velocity                   [m/s]
    \endverbatim

    Required fields:
    \verbatim
      U    | Velocity                   [m/s]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    smoothMatrix1
    {
        // Mandatory entries (unmodifiable)
        type            smoothMatrix;
        selectionMode   all;

        // Optional entries (runtime modifiable)
        smallVolTol     1e-15;      // Volume threshold
        skewTol         0.5;        // Skewness threshold
        nonOrthTol      70.0;       // Non-orthogonality threshold (degrees)
        UNames          (<Uname1> <Uname2> ... <UnameN>);

        // Conditional optional entries (runtime modifiable)

            // when UNames entry is not present
            U               U;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                           | Type     | Reqd  | Dflt
      type      | Type name: smoothMatrix               | word     | yes   | -
      selectionMode | Cell selection mode               | word     | yes   | -
      smallVolTol | Small volume threshold              | scalar   | no    | 0
      skewTol   | Skewness threshold                    | scalar   | no    | 0
      nonOrthTol | Non-orthogonality threshold (degrees) | scalar | no    | 0
      maxResidual | Residual magnitude threshold (|R|)    | scalar | no    | 0
      avgStartIter | Start iteration for residual averaging | label  | no    | 1
      avgEndIter | End iteration for residual averaging     | label  | no    | 0 (disabled)
      UNames    | Names of operand velocity fields      | wordList | no    | -
      U         | Name of operand velocity field        | word     | cndtnl | U
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink
      - \link writeFile.H \endlink

Note
  - When active, this constraint manipulates the system of equations.
    Users should ensure that it is not active when the case is converged
    (steady-state) or during the period of interest (transient) to ensure
    that its presence does not pollute the results.

See also
  - Foam::fv::limitVelocity

SourceFiles
    smoothMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef fv_smoothMatrix_H
#define fv_smoothMatrix_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class smoothMatrix Declaration
\*---------------------------------------------------------------------------*/

class smoothMatrix
:
    public fv::cellSetOption
{
protected:

    // Protected Data

        //- Small volume threshold for bad cell detection
        scalar smallVolTol_;

        //- Skewness threshold for bad cell detection
        scalar skewTol_;

        //- Non-orthogonality threshold for bad cell detection
        scalar nonOrthTol_;

        //- Residual magnitude threshold (activate criterion when > 0)
        scalar residualTol_;

        //- Averaging iteration window [avgStartIter_, avgEndIter_]
        label avgStartIter_;
        label avgEndIter_;

        //- Whether residual-based post-filter has been applied
        bool residualPostApplied_;
        //- Internal counters for iteration-based averaging
        label iterationCounter_;
        label windowSamples_;

        //- Accumulated per-cell residual magnitude across the window
        autoPtr<volScalarField> residualSumFieldPtr_;

        //- Persistent residual-based bad-cell mask after window completes
        boolList residualBadCells_;

        //- Weighting method for neighbour averaging
        word weightingMethod_;

        //- Constraint method: setValues or setReference
        word constraintMethod_;

        //- Enable density-weighted (momentum) averaging
        bool useDensityWeighting_;

        //- Name of density field for weighting
        word rhoName_;

	//- Exclude boundaries (good in lobe pump cases where bad cells lie at the surface)
        bool excludeBoundary_;
        //- Write bad cell list to file
        bool writeBadCells_;

        //- Field of bad cells (1 = bad, 0 = good)
        autoPtr<volScalarField> badCellFieldPtr_;

        //- Write residual field (delta applied to U) to file
        bool writeResiduals_;

        //- Linear-solver residual magnitude field (scalar at cell centres)
        autoPtr<volScalarField> residualFieldPtr_;


        //- Apply neighbour-averaged values to the supplied cells
        void applyNeighbourAveraging(fvMatrix<vector>& eqn);


public:

    //- Runtime type information
    TypeName("smoothMatrix");


    // Constructors

        //- Construct from components
        smoothMatrix
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        smoothMatrix(const smoothMatrix&) = delete;

        //- No copy assignment
        void operator=(const smoothMatrix&) = delete;


    //- Destructor
    virtual ~smoothMatrix() = default;


    // Member Functions

        // Set values directly

            //- Constrain vector matrix
            virtual void constrain(fvMatrix<vector>& eqn, const label fieldi);

     // I-O

            //- Write data
            virtual void writeData(Ostream& os) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
