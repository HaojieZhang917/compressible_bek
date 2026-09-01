# Compatibility entry point for notebooks and existing scripts.
include(joinpath(
    @__DIR__, "..", "..", "..", "RotatingDiskFlow", "src",
    "NeutralCurveRunner.jl",
))

if abspath(PROGRAM_FILE) == @__FILE__
    using .NeutralCurveRunner
    NeutralCurveRunner.main()
end
