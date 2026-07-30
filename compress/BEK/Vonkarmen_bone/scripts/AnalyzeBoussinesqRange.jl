using Printf
using LinearAlgebra
using Statistics

const WORKSPACE_ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTPUT_DIR = get(
    ENV,"NEUTRAL_CURVE_DIR",joinpath(WORKSPACE_ROOT,"neutral_curve_batch"),
)
const TEMPERATURES = (1.00,1.02,1.04,1.06,1.10,1.12,1.14,1.16,1.18,1.20)

function load_curve(path)
    rows = Vector{Vector{Float64}}()
    for line in Iterators.drop(eachline(path),2)
        isempty(strip(line)) && continue
        push!(rows,parse.(Float64,split(strip(line))))
    end
    isempty(rows) && return zeros(0,7)
    return reduce(vcat,permutedims.(rows))
end

function interpolate_R(data,beta)
    order = sortperm(data[:,3])
    x = data[order,3]
    y = data[order,2]
    first(x) <= beta <= last(x) || return NaN
    index = searchsortedlast(x,beta)
    index == length(x) && return y[end]
    index == 0 && return y[1]
    fraction = (beta-x[index])/(x[index+1]-x[index])
    return y[index]+fraction*(y[index+1]-y[index])
end

function interval_error(lopez,compressible,beta_low,beta_high)
    lower = max(beta_low,minimum(lopez[:,3]),minimum(compressible[:,3]))
    upper = min(beta_high,maximum(lopez[:,3]),maximum(compressible[:,3]))
    lower < upper || return nothing
    beta = range(lower,upper; length=301)
    R_lopez = interpolate_R.(Ref(lopez),beta)
    R_compressible = interpolate_R.(Ref(compressible),beta)
    difference = R_lopez-R_compressible
    return (
        beta_low=lower,beta_high=upper,
        relative_l2=norm(difference)/norm(R_compressible),
        mean_absolute=mean(abs.(difference)),
        max_absolute=maximum(abs.(difference)),
        mean_signed=mean(difference),
    )
end

function critical_point(data; beta_low=0.055)
    indices = findall(data[:,3] .>= beta_low)
    isempty(indices) && return nothing
    local_index = argmin(data[indices,2])
    index = indices[local_index]
    return (R=data[index,2],beta=data[index,3])
end

function turning_points(data)
    order = sortperm(data[:,3])
    beta = data[order,3]
    R = data[order,2]
    result = NamedTuple[]
    for index in 2:length(R)-1
        left = R[index]-R[index-1]
        right = R[index+1]-R[index]
        left*right < 0 || continue
        push!(result,(
            kind=left < 0 ? :minimum : :maximum,
            R=R[index],beta=beta[index],
        ))
    end
    return result
end

function model_paths(Tw)
    tag = string(round(Tw; digits=2))
    lopez = joinpath(OUTPUT_DIR,"ome=0.0_Tw=$(tag)_model=lopez.dat")
    branch = joinpath(
        OUTPUT_DIR,"ome=0.0_Tw=$(tag)_model=lopez_branch=typeII.dat",
    )
    compressible = joinpath(
        OUTPUT_DIR,
        "ome=0.0_Tw=$(tag)_model=compressible_Mr=0.3_" *
        "propPert=on_baseProp=variable.dat",
    )
    no_perturbation = joinpath(
        OUTPUT_DIR,
        "ome=0.0_Tw=$(tag)_model=compressible_Mr=0.3_" *
        "propPert=off_baseProp=variable.dat",
    )
    frozen = joinpath(
        OUTPUT_DIR,
        "ome=0.0_Tw=$(tag)_model=compressible_Mr=0.3_" *
        "propPert=off_baseProp=frozen.dat",
    )
    return (
        lopez=lopez,branch=branch,compressible=compressible,
        no_perturbation=no_perturbation,frozen=frozen,
    )
end

println("Property-switch decomposition relative to the full compressible model")
for Tw in (1.02,1.04,1.06,1.10,1.12,1.14,1.16,1.18,1.20)
    paths = model_paths(Tw)
    full = load_curve(paths.compressible)
    no_perturbation = load_curve(paths.no_perturbation)
    frozen = load_curve(paths.frozen)
    no_perturbation_error = interval_error(no_perturbation,full,0.03,0.12)
    frozen_error = interval_error(frozen,full,0.03,0.12)
    critical_full = critical_point(full)
    critical_no_perturbation = critical_point(no_perturbation)
    critical_frozen = critical_point(frozen)
    @printf(
        "Tw=%.2f noPert relL2=%6.3f%% dRc=%+7.3f; frozen relL2=%6.3f%% dRc=%+7.3f\n",
        Tw,100no_perturbation_error.relative_l2,
        critical_no_perturbation.R-critical_full.R,
        100frozen_error.relative_l2,critical_frozen.R-critical_full.R,
    )
end

function print_error(label,error)
    if error === nothing
        @printf("  %-8s unavailable\n",label)
        return
    end
    @printf(
        "  %-8s beta=[%.4f, %.4f] relL2=%7.3f%% mean|dR|=%7.2f max|dR|=%7.2f mean(dR)=%+7.2f\n",
        label,error.beta_low,error.beta_high,100error.relative_l2,
        error.mean_absolute,error.max_absolute,error.mean_signed,
    )
end

for Tw in TEMPERATURES
    paths = model_paths(Tw)
    lopez = load_curve(paths.lopez)
    compressible = load_curve(paths.compressible)
    low_lopez = isfile(paths.branch) ? load_curve(paths.branch) : lopez
    critical_lopez = critical_point(lopez)
    critical_compressible = critical_point(compressible)
    @printf("Tw=%.2f\n",Tw)
    @printf(
        "  critical Lopez: R=%8.3f beta=%.5f; compressible: R=%8.3f beta=%.5f; dR=%+8.3f\n",
        critical_lopez.R,critical_lopez.beta,
        critical_compressible.R,critical_compressible.beta,
        critical_lopez.R-critical_compressible.R,
    )
    print_error("Type-I",interval_error(lopez,compressible,0.055,0.12))
    print_error("Type-II",interval_error(low_lopez,compressible,0.030,0.050))
    for beta in (0.040,0.045,0.070,0.100)
        source = beta <= 0.05 ? low_lopez : lopez
        R_lopez = interpolate_R(source,beta)
        R_compressible = interpolate_R(compressible,beta)
        if isfinite(R_lopez) && isfinite(R_compressible)
            @printf(
                "  beta=%.3f R_L=%8.2f R_C=%8.2f dR=%+8.2f (%+6.2f%%)\n",
                beta,R_lopez,R_compressible,R_lopez-R_compressible,
                100(R_lopez-R_compressible)/R_compressible,
            )
        end
    end
    points = turning_points(lopez)
    if !isempty(points)
        print("  Lopez R(beta) turns:")
        for point in points
            @printf(" %s(R=%.2f,beta=%.5f)",point.kind,point.R,point.beta)
        end
        println()
    end
end
