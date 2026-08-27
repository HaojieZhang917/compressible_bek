# TwoDisk paper revision archive

Archived on 2026-08-26 before starting a new project.

This directory contains the scripts, generated data, plots, logs, reports, and temporary inspection files produced during the manuscript revision. Files were moved here without deleting their contents. Their original names and internal result-directory structures were retained wherever possible.

## Main contents

- Reviewer-specific checks and response support: `r3_c5_validation_results`, `c7_validation_results`, `REFEREE3_C6_C7_*`, and related validation scripts.
- C6 and roughness-width studies: `c_6_*`, `c6_*`, `cr_width_*`, and the associated plotting and report files.
- Configuration and amplitude comparisons: `configuration_*`, `fig21_*`, `six_point_mass_flux_sensitivity_results`, `total_amplitude.jl`, and `typeii_mass_flux_sensitivity.jl`.
- Boussinesq and dynamical-singularity studies: `boussinesq_*`, `three_solution_dynamics`, `dynamical_singularity_comparison`, and their analysis scripts.
- von Karman and spatial-growth studies: `vonkarman_*`, `zarf_spatial_growth_*`, and the associated merge, repair, validation, and plotting scripts.
- Other revision evidence: adjoint-equivalence checks, grid-convergence outputs, exported eigenfunction and mass-flux data, reports, logs, and preview images.
- Temporary manuscript inspection material: `tmp` and `__pycache__`.

## Files intentionally left in the parent directory

The reusable project core remains in `TwoDisk`: the base-flow and stability solvers, the main notebooks, `baseflow_Res1000.npz`, and the long-term `data` directory. The locally modified `receptivity.ipynb` was also retained as a core working notebook.

Some archived scripts may use paths relative to the former `TwoDisk` root. If an old calculation must be rerun, either update those paths or temporarily place the required core solver/data files beside the script.

`legacy_empty_directory` is an empty directory that was originally named with a single space. It was renamed during archiving so that it is visible and manageable on Windows.
