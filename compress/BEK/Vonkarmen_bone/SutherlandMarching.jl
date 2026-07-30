# Compatibility entry point for notebooks and existing scripts.
include(joinpath(
    @__DIR__, "..", "..", "..", "RotatingDiskFlow", "src",
    "SutherlandMarching.jl",
))

if abspath(PROGRAM_FILE) == @__FILE__
    using .SutherlandMarching
    SutherlandMarching.main()
end
