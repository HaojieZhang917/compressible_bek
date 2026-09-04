# Project policy

- New numerical programs should be written in Julia unless the user explicitly requests another language.
- Existing Python post-processing and data-conversion utilities may be retained, but do not add Python dependencies to new numerical solvers.
- Treat everything under `baselines/` as preserved evidence. Do not overwrite baseline data; direct regenerated output to a new directory under `work/results/`.
- Put new shared numerical code in `work/src/` and executable studies in `work/scripts/`.
- Record every new parameter sweep, model definition and source-data dependency in `docs/` before using it in a manuscript claim.
- After every completed numerical calculation, convergence test, continuation,
  model cross-check or mechanism diagnostic, append its conclusion and evidence
  scope to `docs/RESEARCH_PROGRESS_LOG.md` before using it downstream. Preserve
  old log entries; supersede them with a new dated entry instead of rewriting
  the research history.
- Distinguish fully self-consistent model comparisons from frozen-base-flow or term-toggle diagnostics.
- Do not describe similarity-subspace stability as full three-dimensional stability.
- Preserve the original source directories. This workspace was created by copying selected assets, not by moving them.
