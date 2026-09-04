# Boussinesq Acceleration Hierarchy

This workspace supports a new study of how Boussinesq density--acceleration closures alter bifurcation and stability predictions across the Bodewadt--Ekman--von Karman (BEK) family of rotating boundary layers. The intended research chain is

```text
BEK differential rotation
    -> frame-consistent density--acceleration closures
    -> self-consistent base-flow bifurcation maps
    -> fold/cusp mechanism and adjoint sensitivity
    -> similarity-subspace and Type-I/Type-II stability
    -> critical amplification and observable error estimates
```

The workspace was assembled on 2026-08-31 by copying selected scripts and results from existing directories. No original scientific data were moved or deleted.

## Layout

- `baselines/saddle_node/`: self-contained Julia project copied from `Boussinesq_SaddleNode`. This is the primary source for the traditional-model folds, rational-Chebyshev semi-infinite mapping, zero modes, temporal spectra and rotor--stator comparison.
- `baselines/rotating_disk/`: curated rotating-disc workspace copied from `compress/BEK/Vonkarmen_bone`. It contains the local Blackburn operator, neutral-curve runners, tests, active analysis scripts and selected results.
- `work/src/`: new reusable Julia implementation for the density--acceleration hierarchy.
- `work/scripts/`: new derivation checks, parameter sweeps and validation drivers.
- `work/data/`: immutable inputs generated specifically for the new study.
- `work/results/`: regenerated and new outputs. Do not write new results into `baselines/`.
- `docs/CURRENT_STATUS.md`: evidence-based audit of what has and has not been completed.
- `docs/RESEARCH_PROGRESS_LOG.md`: append-only master log of every completed
  calculation, its evidence boundary and its role in the paper's research logic.
- `docs/ASSET_CATALOG.md`: provenance, selection decisions and known limitations.
- `docs/MODEL_DEFINITIONS.md`: agreed Rossby-number convention, closure endpoints and diagnostic coefficients.
- `docs/RESEARCH_PLAN.md`: staged execution plan for the new paper.

## Baseline checks

Saddle-node project integrity:

```bash
cd /home/zhj/Rotating-Flow-ToolKit/Boussinesq_Acceleration_Hierarchy/baselines/saddle_node
python3 check_project.py
```

Blackburn operator and runner:

```bash
cd /home/zhj/Rotating-Flow-ToolKit/Boussinesq_Acceleration_Hierarchy/baselines/rotating_disk
julia --project=. test_blackburn_stability.jl
julia --project=. test_blackburn_runner_smoke.jl
```

The compatibility wrappers in `baselines/rotating_disk/` still resolve the maintained package at `/home/zhj/Rotating-Flow-ToolKit/RotatingDiskFlow` because the copied workspace preserves the required directory depth.

## Immediate development target

The first new implementation should not start with a full neutral-curve sweep. It should:

1. reproduce Lingwood's BEK convention, in which `Ro=-1`, `0` and `1` denote the von Karman, Ekman and Bodewadt limits and `Co=2-Ro-Ro^2`;
2. derive the thermal traditional and canonical Boussinesq equations from one dimensional rotating-frame formulation;
3. implement and verify the self-consistent endpoint models at `Tw=1` and small temperature differences;
4. continue the traditional and canonical fold loci in the `(Ro,Tw)` plane;
5. use term coefficients `eta_j` only as diagnostic homotopy and sensitivity coordinates to identify how folds move or disappear.

Only after this base-flow test succeeds should the same term decomposition be added to the Blackburn perturbation operator and its adjoint sensitivity analysis.
