# Project programming policy

- All newly generated programs in this project should be written in Julia unless the user explicitly requests another language.
- Existing Python-only data conversion, verification and summary utilities are exempt and should not be proactively rewritten.
- Do not introduce Python scripts, PythonCall, PyCall or a Python runtime dependency.
- Reuse `src/BoussinesqSaddleNode.jl` for shared numerics and I/O instead of duplicating solver code.
- Run programs from the project root with `julia --project=. path/to/script.jl`.
- Preserve existing scientific data unless a task explicitly requests regeneration.
