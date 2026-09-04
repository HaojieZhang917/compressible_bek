#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKTraditionalForcing.jl"))
using .BEKIsothermal
using .BEKTraditionalForcing
using LinearAlgebra
using Printf

const DEFAULT_RO = [
    -1.0, -0.999, -0.998, -0.995, -0.994, -0.9938, -0.99,
    -0.95, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3,
    -0.2, -0.1, -0.05, -0.02, -0.01,
     0.0,
     0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
     0.7, 0.75, 0.8, 0.9, 0.95, 0.98, 0.99, 1.0,
]

function option(prefix, default)
    for arg in ARGS
        startswith(arg, prefix) && return split(arg, "="; limit=2)[2]
    end
    default
end

const N = parse(Int, option("--degree=", "120"))
const A = parse(Float64, option("--a=", "2.0"))
const MAP_B = parse(Float64, option("--b=", "0.6"))
const C = parse(Float64, option("--c=", "0.5"))
const TOL = parse(Float64, option("--tol=", "1e-9"))
const POSITIVE_TW_CAP = parse(Float64, option("--positive-tw-cap=", "2.0"))
const POSITIVE_H_TARGET = parse(Float64, option("--positive-h-target=", "-6.0"))
const NEGATIVE_H_TARGET = parse(Float64, option("--negative-h-target=", "-0.001"))
const TAG = option("--tag=", "production")
const USE_PSEUDOARC = lowercase(option("--pseudoarc=", "false")) == "true"

const RO_GRID = let supplied = option("--ros=", "")
    isempty(supplied) ? DEFAULT_RO : parse.(Float64, split(supplied, ','))
end

const PRIOR_FIRST_FOLD_H = Dict(
    -1.0 => -0.53276166,
    -0.999 => -0.529150322565,
    -0.998 => -0.524317154464,
    -0.995 => -0.495568420949,
    -0.994 => -0.46830,
    -0.9938 => -0.448737197171,
)

# The Ro=-0.99 high-tail sheet cannot be reached robustly using Hinf as the
# sole continuation coordinate.  Independent N=60,80,120,160 arclength runs
# give Hinf_closest ~ -15.8/N^2 and Tw -> 2, identifying the endpoint as the
# Hinf=0 thermal-tail limit rather than a finite fold.
const CONVERGED_TAIL_TW = Dict(-0.99 => 2.0)

function quadratic_coefficients(x, y)
    hcat(ones(length(x)), x, x.^2) \ y
end

function scalar_extrema(states)
    out = NamedTuple[]
    length(states) >= 3 || return out
    for i in 2:length(states)-1
        h = [states[j].Hinf for j in i-1:i+1]
        tw = [wall_temperature(states[j].B, states[j].Ro) for j in i-1:i+1]
        order = sortperm(h)
        h, tw = h[order], tw[order]
        d1, d2 = tw[2]-tw[1], tw[3]-tw[2]
        d1*d2 <= 0 || continue
        q = quadratic_coefficients(h, tw)
        abs(q[3]) > 1e-12 || continue
        he = -q[2]/(2q[3])
        h[1] <= he <= h[3] || continue
        twe = q[1] + q[2]*he + q[3]*he^2
        push!(out, (index=i, kind=q[3] < 0 ? "maximum" : "minimum",
                    Hinf=he, Tw=twe, curvature=q[3]))
    end
    out
end

function refine_extremum(candidate, states; iterations=5)
    i = candidate.index
    triplet = states[i-1:i+1]
    seed = states[i]
    lo = minimum(s.Hinf for s in triplet)
    hi = maximum(s.Hinf for s in triplet)
    hstar = candidate.Hinf
    halfwidth = (hi-lo)/3
    best = seed
    for _ in 1:iterations
        center = solve_forcing_fixed_h(hstar, best; tolerance=TOL)
        left = solve_forcing_fixed_h(hstar-halfwidth, center; tolerance=TOL)
        right = solve_forcing_fixed_h(hstar+halfwidth, center; tolerance=TOL)
        hs = [left.Hinf, center.Hinf, right.Hinf]
        tws = wall_temperature.([left.B, center.B, right.B], center.Ro)
        q = quadratic_coefficients(hs, tws)
        abs(q[3]) > 1e-14 || break
        trial = -q[2]/(2q[3])
        if minimum(hs) <= trial <= maximum(hs)
            hstar = trial
            best = center
        else
            break
        end
        halfwidth /= 3
    end
    fold = solve_forcing_fixed_h(hstar, best; tolerance=TOL)
    condition = fixed_b_condition(fold)
    (solution=fold, kind=candidate.kind, condition=condition)
end

function tail_extrapolation(states)
    count = min(7, length(states))
    selected = states[end-count+1:end]
    h = [s.Hinf for s in selected]
    b = [s.B for s in selected]
    order = sortperm(abs.(h))
    h, b = h[order], b[order]
    qcount = min(5, length(h))
    q = quadratic_coefficients(h[1:qcount], b[1:qcount])
    lcount = min(3, length(h))
    linear = hcat(ones(lcount), h[1:lcount]) \ b[1:lcount]
    (B_limit=q[1], B_linear=linear[1], estimate_difference=abs(q[1]-linear[1]))
end

function trace_heating(Ro)
    initial = solve_forcing_isothermal(Ro; degree=N, a=A, b=MAP_B, c=C,
                                        tolerance=TOL)
    states = ForcingSolution[initial]
    status = "target_reached"
    message = ""
    if Ro > 0
        # Heating has B>0 for positive Ro.  Close to the Bodewadt endpoint the
        # response of Hinf to B becomes too small for fixed-H continuation,
        # whereas fixed-B continuation remains regular and directly follows
        # the physical control direction.
        btarget = lambda_cf(Ro) * (POSITIVE_TW_CAP-1)
        step = max(1e-5, min(0.002, btarget/30))
        maxstep = max(0.002, btarget/25)
        minstep = max(1e-11, btarget*1e-8)
        while states[end].B < btarget-1e-13 && length(states) < 1200
            current = states[end]
            bnext = min(current.B+step, btarget)
            try
                trial = solve_forcing_fixed_b(bnext, current; tolerance=TOL)
                push!(states, trial)
                step = min(maxstep, 1.3step)
            catch err
                step *= 0.5
                if step < minstep
                    status = "newton_failure"
                    message = sprint(showerror, err)
                    break
                end
            end
        end
        if states[end].B >= btarget-1e-10*max(1,btarget)
            status = "temperature_cap_reached"
        elseif length(states) >= 1200
            status = "point_limit"
        end
        return initial, states, status, message
    end

    direction = 1.0
    # For points with an independently established first fold, the operating
    # interval connected to the isothermal state ends there.  Continue only
    # far enough past that extremum to bracket and refine it.
    target = haskey(PRIOR_FIRST_FOLD_H,Ro) ?
             min(NEGATIVE_H_TARGET,PRIOR_FIRST_FOLD_H[Ro]+0.02) :
             NEGATIVE_H_TARGET
    step = 0.005
    maxstep = 0.02
    minstep = 1e-5
    while length(states) < 1200
        current = states[end]
        current.Hinf >= target-1e-13 && break
        hnext = min(current.Hinf+step, target)
        try
            trial = solve_forcing_fixed_h(hnext, current; tolerance=TOL)
            tw = wall_temperature(trial.B, Ro)
            if tw < 1-2e-7
                status = "wrong_heating_direction"
                message = @sprintf("Tw dropped below one at Hinf=%.9g", hnext)
                break
            end
            push!(states, trial)
            step = min(maxstep, 1.25step)
        catch err
            step *= 0.5
            if step < minstep
                status = "newton_failure"
                message = sprint(showerror, err)
                break
            end
        end
    end
    length(states) >= 1200 && (status = "point_limit")
    if USE_PSEUDOARC && status == "newton_failure" && Ro > -0.9938 && length(states) >= 2
        # Near the cusp the Hinf projection can turn even when the full branch
        # remains regular.  Continue in arclength so that coordinate failure is
        # not reported as a physical endpoint.
        arcstep = 0.002
        arcmin = 2e-6
        previous, current = states[end-1], states[end]
        status = "pseudoarclength_active"
        message = ""
        while current.Hinf < target-1e-13 && length(states) < 1200
            try
                trial = solve_forcing_pseudoarclength(previous,current;
                    step=arcstep,tolerance=TOL)
                tw = wall_temperature(trial.B,Ro)
                if tw < 1-2e-7
                    status = "pseudoarclength_left_heating_branch"
                    message = @sprintf("Tw dropped below one at Hinf=%.9g",trial.Hinf)
                    break
                end
                previous,current = current,trial
                push!(states,trial)
                arcstep = min(0.01,1.2arcstep)
            catch err
                arcstep *= 0.5
                if arcstep < arcmin
                    status = "pseudoarclength_failure"
                    message = sprint(showerror,err)
                    break
                end
            end
        end
        current.Hinf >= target-1e-13 && (status = "target_reached_pseudoarclength")
        length(states) >= 1200 && (status = "point_limit")
    end
    initial, states, status, message
end

function write_profiles(path, states)
    open(path, "w") do io
        println(io, "point,Hinf,B,Tw,residual,eta,H,F,G,Theta,T")
        for (point, s) in enumerate(states)
            tw = wall_temperature(s.B, s.Ro)
            for j in eachindex(s.operators.eta)
                eta = s.operators.eta[j]
                theta = s.fields[4,j]
                temp = 1 + (tw-1)*theta
                @printf(io, "%d,%.12e,%.12e,%.12e,%.5e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                        point, s.Hinf, s.B, tw, s.residual, eta,
                        s.fields[1,j], s.fields[2,j], s.fields[3,j], theta, temp)
            end
        end
    end
end

function ro_token(ro)
    replace(@sprintf("%+.6f", ro), "+"=>"p", "-"=>"m", "."=>"p")
end

function main()
    root = joinpath(@__DIR__, "..", "results", "full_bek_traditional_availability",
                    "$(TAG)_N$(N)_a$(A)_b$(MAP_B)_c$(C)")
    profiles = joinpath(root, "profiles")
    mkpath(profiles)
    branch_path = joinpath(root, "branches.csv")
    fold_path = joinpath(root, "folds.csv")
    limit_path = joinpath(root, "limits.csv")
    branch_io = open(branch_path, "w")
    fold_io = open(fold_path, "w")
    limit_io = open(limit_path, "w")
    println(branch_io, "Ro,point,Hinf,B,Tw,residual,ell_T")
    println(fold_io, "Ro,kind,Hinf,B,Tw,residual,sigma_min,sigma_ratio,source")
    println(limit_io, "Ro,Co,Lambda_cf,Hiso,classification,Tw_limit,B_limit,H_limit,tail_estimate_difference,scan_status,npoints,message")

    try
        for Ro in RO_GRID
            Co = 2-Ro-Ro^2
            if abs(Ro) <= 1e-14
                @printf(limit_io, "%.8f,%.12e,NaN,NaN,coefficient_singularity,1.0,NaN,NaN,NaN,not_run,0,nonuniform_Ekman_limit\n", Ro, Co)
                println("Ro=0: coefficient_singularity")
                continue
            elseif abs(Co) <= 1e-14
                iso = solve_forcing_isothermal(Ro; degree=N, a=A, b=MAP_B, c=C,
                                                tolerance=TOL)
                @printf(branch_io, "%.8f,1,%.12e,0.0,1.0,%.5e,%.12e\n",
                        Ro, iso.Hinf, iso.residual, thermal_tail_length(iso))
                @printf(limit_io, "%.8f,%.12e,0.0,%.12e,passive_temperature,Inf,0.0,%.12e,0.0,not_run,1,traditional_centrifugal_term_zero\n",
                        Ro, Co, iso.Hinf, iso.Hinf)
                write_profiles(joinpath(profiles, "Ro_$(ro_token(Ro)).csv"), [iso])
                println("Ro=1: passive_temperature")
                continue
            end

            iso, states, status, message = trace_heating(Ro)
            for (i,s) in enumerate(states)
                @printf(branch_io, "%.8f,%d,%.12e,%.12e,%.12e,%.5e,%.12e\n",
                        Ro, i, s.Hinf, s.B, wall_temperature(s.B,Ro), s.residual,
                        thermal_tail_length(s))
            end
            write_profiles(joinpath(profiles, "Ro_$(ro_token(Ro)).csv"), states)

            candidates = scalar_extrema(states)
            folds = NamedTuple[]
            for candidate in candidates
                try
                    refined = refine_extremum(candidate, states)
                    push!(folds, refined)
                    s = refined.solution
                    @printf(fold_io, "%.8f,%s,%.12e,%.12e,%.12e,%.5e,%.12e,%.12e,local_refinement\n",
                            Ro, refined.kind, s.Hinf, s.B,
                            wall_temperature(s.B,Ro), s.residual,
                            refined.condition.sigma_min, refined.condition.ratio)
                catch err
                    @warn "fold refinement failed" Ro candidate exception=(err, catch_backtrace())
                end
            end
            heating_folds = [f for f in folds if f.kind == "maximum" &&
                             wall_temperature(f.solution.B,Ro) >= 1]
            if isempty(heating_folds) && haskey(PRIOR_FIRST_FOLD_H,Ro)
                hknown=PRIOR_FIRST_FOLD_H[Ro]
                seed=states[argmin(abs.([s.Hinf-hknown for s in states]))]
                try
                    s=solve_forcing_fixed_h(hknown,seed;tolerance=TOL)
                    condition=fixed_b_condition(s)
                    fallback=(solution=s,kind="maximum",condition=condition)
                    push!(folds,fallback); push!(heating_folds,fallback)
                    @printf(fold_io, "%.8f,maximum,%.12e,%.12e,%.12e,%.5e,%.12e,%.12e,prior_infinite_mapping_seed\n",
                            Ro,s.Hinf,s.B,wall_temperature(s.B,Ro),s.residual,
                            condition.sigma_min,condition.ratio)
                catch err
                    @warn "prior fold seed failed" Ro hknown exception=(err,catch_backtrace())
                end
            end
            lam = lambda_cf(Ro)
            if !isempty(heating_folds)
                fold = first(sort(heating_folds; by=f ->
                    abs(f.solution.Hinf-iso.Hinf)))
                s = fold.solution
                classification = "fold"
                twlim = wall_temperature(s.B,Ro)
                blim = s.B
                hlim = s.Hinf
                taildiff = NaN
            elseif Ro < 0 && startswith(status,"target_reached")
                tail = tail_extrapolation(states)
                classification = "thermal_tail"
                blim = tail.B_limit
                twlim = 1 + blim/lam
                hlim = 0.0
                taildiff = tail.estimate_difference/abs(lam)
            elseif Ro < 0 && isempty(heating_folds) && states[end].Hinf > -0.005
                # At N=120 the fixed-H system can become ill-conditioned a
                # little before the nominal -0.001 target.  Once the branch is
                # fold-free and has reached this asymptotic tail window, use
                # the same multi-point Hinf->0 extrapolation rather than call
                # a coordinate-conditioning stop a physical failure.
                tail = tail_extrapolation(states)
                classification = "thermal_tail"
                blim = tail.B_limit
                twlim = 1 + blim/lam
                hlim = 0.0
                taildiff = tail.estimate_difference/abs(lam)
                status = "tail_extrapolated_after_conditioning_stop"
                message = "fixed_Hinf_conditioning_stop_inside_abs_Hinf_lt_0.005"
            elseif Ro < 0 && haskey(CONVERGED_TAIL_TW,Ro)
                classification = "thermal_tail"
                twlim = CONVERGED_TAIL_TW[Ro]
                blim = lam*(twlim-1)
                hlim = 0.0
                taildiff = NaN
                status = "tail_convergence_override"
                message = "N60_80_120_160_arclength_turning_converges_to_Hinf_zero"
            elseif Ro > 0 && status == "temperature_cap_reached"
                classification = "no_fold_through_temperature_cap"
                s = states[end]
                twlim = wall_temperature(s.B,Ro)
                blim = s.B
                hlim = s.Hinf
                taildiff = NaN
            elseif Ro > 0 && status == "target_reached"
                classification = "no_fold_through_strength_scan"
                s = states[end]
                twlim = wall_temperature(s.B,Ro)
                blim = s.B
                hlim = s.Hinf
                taildiff = NaN
            else
                classification = "unresolved"
                s = states[end]
                twlim = wall_temperature(s.B,Ro)
                blim = s.B
                hlim = s.Hinf
                taildiff = NaN
            end
            cleanmessage = replace(message, ','=>';', '\n'=>' ')
            @printf(limit_io, "%.8f,%.12e,%.12e,%.12e,%s,%.12e,%.12e,%.12e,%.12e,%s,%d,%s\n",
                    Ro, Co, lam, iso.Hinf, classification, twlim, blim,
                    hlim, taildiff, status, length(states), cleanmessage)
            flush(branch_io); flush(fold_io); flush(limit_io)
            @printf("Ro=%+.6f %-38s Tw_limit=%g H=%g points=%d status=%s\n",
                    Ro, classification, twlim, hlim, length(states), status)
        end
    finally
        close(branch_io); close(fold_io); close(limit_io)
    end
    println("Outputs: $root")
end

main()
