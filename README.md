
# smoothMatrix

**smoothMatrix** is an OpenFOAM `fvOption` developed by [IANUS Simulation GmbH](https://ianus-simulation.de/en/) and designed to enhance numerical stability in CFD simulations. It automatically detects problematic cells based on mesh quality and applies neighbor-averaged corrections to prevent solution divergence while maintaining physical consistency.

## 🚀 Overview

In complex CFD simulations, mesh deficiencies—such as high skewness, non-orthogonality, or extreme aspect ratios—often lead to local instabilities that can crash a solver. `smoothMatrix` intervenes by identifying these "bad" cells and applying a weighted averaging technique to the velocity field, ensuring the source terms are treated consistently to preserve local flow behavior.

## ✨ Key Features

* **Different cell detection criterion:** Identifies problematic regions using mesh quality metrics (volume, skewness, non-orthogonality, initial residuals).
* **Multiple Weighting Schemes:**
* `inverseDistance`: Weights based on proximity to neighbors.
* `volumeWeighted`: Accounts for cell volume keeping momentum conserve.
* `uniform`: Simple arithmetic averaging.


* **Physical Consistency:** Includes optional **density-weighted averaging** for compressible flows and source-term-consistent corrections.
* **Constraint Methods:** Features both **hardBlend** (strict correction) and **softBlend** (gradual) blending to minimize numerical dissipation.

---

## 🛠 Installation

1. **Clone the repository:**
```bash
git clone git@bitbucket.org:hmarschall/smoothmatrix.git

```


2. **Compile the library:**
With a sourced working OpenFoam-v2412 navigate to the source directory and run :
```bash
cd smoothMatrix
wmake libso

```



## 📖 Usage

To use `smoothMatrix`, add it to your `constant/fvOptions` file:

```cpp
smoothMatrix1
{
    type         smoothMatrix;    
    active       yes;
    selectionMode    all;
    
    rhoName	rho; //name of the density for intrinsic quantities
    useDensityWeighting true; //default is false
    
    smallVolTol  1e-18; 
    skewTol      2.5;  //face averaged skewness value at the cell center
    nonOrthTol  60.0;  //face averaged non-Orthogonal value at the cell center
    U            U;
    excludeBoundary	false; //exclude boundary cells in selection process
   
    // Residual-based averaging window (iterations)
    avgStartIter    10;        // start from first iteration

    avgEndIter      1000;       // end of averaging window

    // Threshold on AVERAGE residual over the window
    maxResidual     1e-3;

    constraintMethod	hardBlend; //softBlend 
    writeBadCells	false; //write out cells where correction is applied
    weightingMethod	volumeWeighted;//inverseDistance;   //inverseDistanceSquared; uniform;
}
```

---

## 🗺 Roadmap: The Future of smoothMatrix

We are actively working on moving beyond purely geometric indicators. The next phase of development includes:

* **Hybrid Error Metrics:** Integrating solution-based error indicators (e.g., high velocity gradients or pressure residuals) alongside mesh quality.
* **Dynamic Adaptation:** Enabling the tool to "wake up" only when local solution errors exceed a specific threshold, even in "good" mesh areas.
* **Expanded Field Support:** Extending stabilization logic to scalars and turbulence model equations.

## 🤝 Contributing

Contributions are welcome! Whether it's bug reports, feature requests, or pull requests, feel free to open an issue or submit a contribution to help make `smoothMatrix` more robust.

---
