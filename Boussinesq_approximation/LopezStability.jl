# Compatibility entry point. New code should use `RotatingDiskFlow.LopezStability`.
include(joinpath(@__DIR__, "..", "..", "..", "RotatingDiskFlow", "src", "LopezStability.jl"))
