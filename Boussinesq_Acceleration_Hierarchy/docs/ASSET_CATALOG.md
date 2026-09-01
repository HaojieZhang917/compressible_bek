# Asset catalog and provenance

## 1. Saddle-node baseline

Destination: `baselines/saddle_node/`

Source: `/home/zhj/Rotating-Flow-ToolKit/Boussinesq_SaddleNode/`

Copied as a complete self-contained project. It contains:

- native Julia common modules and solvers;
- finite-domain and rational-Chebyshev von Karman continuation;
- rotor--stator continuation and three-solution dynamics;
- zero-mode, temporal-spectrum and square-root-law comparisons;
- reports, source/data catalogs and manuscript introduction material.

This is the preferred saddle-node baseline. The older Python fold scripts in `Vonkarmen_bone/scripts/` were not copied into the new rotating-disc baseline because the native Julia implementation and true-infinity mapping supersede them.

## 2. Rotating-disc/Blackburn baseline

Destination: `baselines/rotating_disk/`

Source: `/home/zhj/Rotating-Flow-ToolKit/compress/BEK/Vonkarmen_bone/`

### Copied code

- local Blackburn operator and neutral-curve runner;
- Blackburn batch, extended-temperature and Type-II convergence drivers;
- operator and runner tests;
- compatibility wrappers for the maintained `RotatingDiskFlow` package;
- active neutral-curve, benchmark, grid and model-comparison scripts;
- derivation and workflow documentation;
- the main Boussinesq notebook as a historical interactive reference.

### Copied primary results

- `blackburn_neutral_curve_batch/`;
- `blackburn_neutral_results/`;
- `blackburn_frozen_neutral_results/`;
- `zero_frequency_complete_error_analysis/`;
- `zero_frequency_extended_error_analysis/`;
- `typeII_grid_convergence_results/`;
- `grid_independence/`;
- compressible neutral-curve reference files from `neutral_curve_batch/`.

### Copied secondary results

- selected Lopez-versus-traditional base-flow summaries and figures;
- selected Lopez-versus-compressible summaries and figures;
- selected compact base-flow feature/error tables;
- Chapman-versus-Sutherland comparison;
- one Sutherland radial-marching case.

### Deliberately not copied

- `archive/`, generated duplicates and paused positive-frequency work;
- old finite-domain fold scripts superseded by `baselines/saddle_node/`;
- large duplicate raw profile/interpolation tables that can be regenerated;
- the integrated Lopez neutral-curve dataset whose own report states that the Lopez curves predate removal of an obsolete axial thermal-feedback entry;
- unrelated TwoDisk manuscript-revision, receptivity and roughness-response assets.

The omitted files remain unchanged at their original locations.

## 3. Important limitations attached to copied results

1. The `baseflow_comparison_data` report describes a SciPy reimplementation and a compressible construction that does not fully recompute temperature-coupled radial and azimuthal momentum; use it as a diagnostic, not as the final reference.
2. Lopez and Blackburn use a linear density law. Results at large wall-temperature ratios are failure diagnostics rather than quantitative Boussinesq predictions.
3. Frozen-property and frozen-base-flow cases isolate mechanisms but are not self-consistent model validation.
4. The Type-II transition is bracketed only by sampled temperatures and parameter resolution; it is not an exact bifurcation boundary.
5. Temporal stability in the saddle-node project is restricted to the axisymmetric similarity subspace.

## 4. Original directories preserved

The following sources were read and copied without moving or deleting their contents:

- `/home/zhj/Rotating-Flow-ToolKit/Boussinesq_SaddleNode/`
- `/home/zhj/Rotating-Flow-ToolKit/compress/BEK/Vonkarmen_bone/`
- `/home/zhj/Rotating-Flow-ToolKit/TwoDisk/paper_revision_archive_2026-08-26/`

The TwoDisk archive was audited but not copied separately because its saddle-node material is already reorganised and provenance-tracked in `baselines/saddle_node/`.

