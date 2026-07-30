# Vonkarmen_bone Research Workspace

This directory is the daily research workspace for rotating-disk base-flow,
stability, and neutral-curve studies. The reusable Julia package is maintained
separately at:

```text
/home/zhj/Rotating-Flow-ToolKit/RotatingDiskFlow
```

## Layout

| Path | Purpose |
|---|---|
| `Bounssinesq.ipynb` | Main interactive notebook |
| `CRD_STA.jl`, `LopezBaseflow.jl`, `LopezStability.jl` | Compatibility entry points for existing notebooks |
| `NeutralContinuation.jl`, `NeutralCurveRunner.jl` | Compatibility entry points for neutral-curve calculations |
| `SutherlandMarching.jl` | Compatibility entry point for the marching base-flow solver |
| `scripts/` | Reusable analysis, continuation, validation, and export commands |
| `docs/` | Mathematical derivations and workflow notes |
| `archive/` | Superseded scripts and preserved auxiliary artifacts |
| `*_comparison/`, `*_analysis/`, `neutral_curve_*` | Scientific output and post-processing data |

The compatibility files load the maintained implementations from
`../../../RotatingDiskFlow/src/`, so existing notebook cells do not need to be
rewritten.

## Julia Setup

For package-oriented work:

```julia
using Pkg
Pkg.activate("/home/zhj/Rotating-Flow-ToolKit/RotatingDiskFlow")
using RotatingDiskFlow
```

Existing notebook code may continue to use:

```julia
include("NeutralCurveRunner.jl")
using .NeutralCurveRunner
```

For scripts in this workspace, activate the local environment:

```bash
julia --project=. scripts/CheckMalikBenchmarks.jl
```

See [`scripts/README.md`](scripts/README.md) for the available commands and
[`archive/README.md`](archive/README.md) for the retention policy.

A step-by-step Chinese guide for generating base flows and neutral curves is
available in [`docs/workflow-guide-cn.md`](docs/workflow-guide-cn.md).

For continuing this research in a new Codex conversation, load
[`CONVERSATION_HANDOFF.txt`](CONVERSATION_HANDOFF.txt) first.

## Output Policy

Computed data remain in their existing result directories so that notebooks and
paper figures keep stable paths. Runtime logs, Python bytecode, notebook caches,
and regenerated neutral-curve directories are ignored by Git.

Release-package documentation is available in:

```text
/home/zhj/Rotating-Flow-ToolKit/RotatingDiskFlow/README.md
/home/zhj/Rotating-Flow-ToolKit/RotatingDiskFlow/docs/neutral-curves.md
```
