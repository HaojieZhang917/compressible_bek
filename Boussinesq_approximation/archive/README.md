# Archive

This directory contains material retained for provenance but not used by the
current workflow.

## `legacy_scripts/`

- The three `compute_lopez_*` and `resume_lopez_*` programs predate
  `NeutralCurveRunner.jl` and are retained to reproduce older continuation runs.
- `Stability.jl` is the superseded stability implementation.
- `ShootBone.py` is the superseded shooting implementation.
- `generate_baseflow_comparison.jl` depends on the historical `Bone.py`, which
  is not present in this repository. Use
  `scripts/generate_boussinesq_compressible_report.py` for the maintained path.

## Auxiliary Artifacts

- `generated_duplicates/` contains output variants retained because they were
  not byte-identical to the primary result.
- `tecplot_layouts/` contains optional Tecplot layout files moved out of the
  scientific data directory.

The one-off diagnostics `audit_crd_sta.jl`, `diagnose_crd_type2.jl`, and
`scan_lopez_modes.jl` were removed during cleanup. Their useful checks are
covered by the maintained benchmark, grid-independence, and package test tools.
