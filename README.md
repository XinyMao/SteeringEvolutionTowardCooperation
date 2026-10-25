# SteeringEvolutionTowardCooperation
Code for paper "Steering Evolution Toward Cooperation Despite Unfavorable Selection"

The repository contains:

- analytical calculations for sparse and finite steering;
- Python simulations on random-regular graphs and square lattices;
- bootstrap and threshold-crossing analysis;
- plot code for MATLAB figures

## 1. Overview

The repository provides code associated with the numerical analyses and figure preparation presented in the manuscript.

The code is organized into two categories:

- **Python:** Entry points for numerical simulations of sparse, finite, and adaptive intervention.
- **MATLAB:** Scripts for generating and exporting figures from previously computed simulation results.

The MATLAB plotting scripts use processed numerical data as inputs. They do not independently regenerate the underlying simulation trajectories.

## 2. Software requirements

The computational workflow uses:

- MATLAB R2026a (MathWorks)
- Python 3.14.0
- Python scientific-computing libraries, including NumPy, SciPy, pandas, and Matplotlib

These are the software versions reported for the manuscript. Compatibility with other versions has not been systematically verified.

Some MATLAB functions use features for table processing, graph visualization, and vector graphics export.

## 3. Repository contents

### Python simulation entry points

| File | Description |
|---|---|
| `sparse.py` | Entry point for sparse-intervention simulations |
| `finite.py` | Entry point for finite-intervention simulations |
| `adaptive.py` | Entry point for adaptive-relocation simulations associated with Fig. 5 |

The Python files delegate execution to the corresponding simulation implementations in the main project.

### MATLAB figure scripts

| File | Description |
|---|---|
| `SparseIntervention.m` | Visualization of sparse-intervention results |
| `FiniteInervention.m` | Visualization of finite-intervention results |
| `InterventionPlacement.m` | Visualization of intervention-placement and exposure results |
| `AdaptiveInterventionExposure.m` | Visualization of adaptive exposure landscapes |
| `AdaptiveInterventionThreshold.m` | Visualization of adaptive intervention thresholds |
| `AdaptiveInterventionSpatialCorrelation.m` | Visualization of adaptive dynamics and spatial configurations |

The filenames above follow the supplied files. In particular, `FiniteInervention.m` retains its existing spelling.

## 4. Numerical simulation workflows

### 4.1 Sparse intervention

The sparse-intervention analysis investigates the effects of a small intervention fraction on evolutionary cooperation dynamics.

The Python entry point is:

`python/sparse.py`

It delegates execution to:

`scripts/run_sparse_drift.py`

The corresponding MATLAB plotting script generates the sparse-intervention figure from processed simulation and bootstrap results.

### 4.2 Finite intervention

The finite-intervention analysis examines evolutionary drift and cooperation thresholds beyond the sparse-intervention regime.

The Python entry point is:

`python/finite.py`

It delegates execution to:

`scripts/run_finite_steering.py`

The corresponding MATLAB script uses processed numerical results to construct the figure associated with finite intervention.

### 4.3 Intervention placement and exposure

The intervention-placement analysis examines how the spatial arrangement of intervention sites affects their exposure to defection and the resulting evolutionary drift.

The MATLAB plotting script reads previously computed configuration-level results and statistical summaries.

### 4.4 Adaptive intervention

The adaptive-intervention analysis investigates the effects of relocating intervention sites in response to the evolving population configuration.

The Python entry point is:

`python/adaptive.py`

It supports three experiment selectors:

- `fig5`: adaptive exposure analysis
- `fig6`: adaptive threshold analysis
- `fig7`: phase-space and spatial-configuration analysis

The corresponding MATLAB scripts generate the figure panels from previously computed results.

## 5. Running the code

### Python

The Python entry points are intended to be run from the main project repository, where the canonical simulation scripts and their dependencies are available.

For example:

```bash
python python/sparse.py --help
python python/finite.py --help
python python/adaptive.py --experiment fig5c --help
```

These commands require the corresponding simulation runners to be present in the expected project structure.

The experiment parameters and available command-line options are defined by the underlying simulation runners.

### MATLAB

The MATLAB scripts can be called after adding their directory to the MATLAB search path.

For example, the sparse-intervention plotting function can be invoked using its original function name:

```matlab
plot_sparse_final_figure2
```

These commands assume that the required processed input files are available at the paths expected by the scripts.

For scripts that use fixed project-relative paths, the original directory structure must be preserved or the input paths must be adjusted.

## 6. Input data and figure reproduction

The MATLAB scripts generally read processed simulation results stored in CSV, JSON, or MAT files.

Examples of expected inputs include:

- Sparse-intervention pooled and cell-level bootstrap tables
- Finite-intervention drift and threshold estimates
- Intervention-placement configuration summaries
- Adaptive exposure summaries
- Adaptive local and global threshold estimates
- Phase-space trajectories and network snapshots

Several plotting scripts generate editable MATLAB figures and vector PDF outputs. Some scripts also export PNG images, CSV summaries, or statistical information.

**Important:** The MATLAB plotting scripts require their corresponding processed input data. The plotting code alone is not sufficient to reproduce the manuscript figures if these input files are absent.

## 7. Reproducibility

The code is provided to document the computational procedures associated with the manuscript.

The numerical experiments and theoretical comparisons are described in the main text and Supplementary Information.

To reproduce the computational results, users should have access to:

1. The relevant simulation implementations and dependencies.
2. The parameter settings described in the manuscript and Supplementary Information.
3. The processed input data required by the figure-generation scripts.

The Python entry points do not contain independent copies of the full simulation implementations. They invoke the corresponding runners from the main project.

The MATLAB scripts are intended primarily for visualization and post-processing, rather than rerunning the numerical simulations.

For questions concerning the computational methods, numerical data, or code, please contact the corresponding author of the associated manuscript.
