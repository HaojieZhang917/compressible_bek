# Workspace Scripts

Run these commands from the `Vonkarmen_bone` workspace. Every script resolves
input and output paths from the workspace root, so invocation from another
working directory is also supported.

## Neutral Curves

```bash
julia --project=. scripts/RecomputeBenchmarkCurves.jl
julia --project=. scripts/ComputeLopezTypeI.jl 1.08
julia --project=. scripts/ComputeLopezTypeII.jl 1.18 1.20
julia --project=. scripts/MergeLopezNeutralBranches.jl
julia --project=. scripts/IntegrateNeutralCurves.jl
julia --project=. scripts/AnalyzeNeutralCriticalErrors.jl
```

`AnalyzeBoussinesqRange.jl` performs the model and property-switch comparison.
`NeutralCurveRunner.jl` in the workspace root remains the preferred programmatic
entry point for new continuation cases.

## Validation

```bash
julia --project=. scripts/CheckMalikBenchmarks.jl
julia --project=. scripts/LopezGridIndependence.jl
```

The first command checks the Malik Type-I and Type-II benchmark locations. The
second performs the more expensive collocation-grid study.

## Base Flow And Post-Processing

```bash
python scripts/compare_lopez_boussinesq.py
python scripts/analyze_compressible_lopez.py
python scripts/error_convergence.py
python scripts/investigate_boussinesq_fold.py
python scripts/compute_domain_branches.py
python scripts/generate_boussinesq_compressible_report.py
python scripts/compare_chapman_sutherland.py
python scripts/convert_lopez_csv_to_tecplot.py
```

These scripts write into the existing `*_comparison`, `*_analysis`, and
base-flow data directories. They do not place generated files in `scripts/`.

For parameter examples, output formats, branch merging, and the recommended
calculation order, see [`../docs/workflow-guide-cn.md`](../docs/workflow-guide-cn.md).
