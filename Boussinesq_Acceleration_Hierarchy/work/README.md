# Active work area

Use this directory for all new calculations.

- `src/`: shared Julia modules for acceleration decomposition, model coefficients, continuation, fold sensitivity and modal sensitivity.
- `scripts/`: executable studies and verification cases.
- `data/`: new immutable inputs or externally generated reference data.
- `results/`: outputs from new runs, grouped by script and parameter set.

The contents of `baselines/` are evidence snapshots. Import or read them, but do not overwrite them.

Before running a new production sweep, add a record under `docs/` containing the governing model endpoint, any diagnostic `eta_j` path, `Ro` and `Tw` samples, numerical resolution and tolerances, source-data dependencies and the intended `work/results/` output directory.

Current implementation priority:

1. BEK dimensional/nondimensional derivation and regression cases;
2. one self-consistent Julia steady residual for the traditional and canonical endpoints;
3. pseudo-arclength fold continuation and cusp checks;
4. left/right fold-null-vector sensitivity;
5. stability calculations only after the base-flow mechanism passes validation.

The traditional-model full-family base-flow foundation is now indexed by
`results/full_bek_traditional_availability/production_v2_N120_a2.0_b0.6_c0.5/final_applicability.csv`.
Use that file and its sibling `profiles/` directory for subsequent work; the
`smoke_*`, `exploratory_*`, and interrupted `production_N120_*` directories
are diagnostic evidence rather than authoritative inputs.
