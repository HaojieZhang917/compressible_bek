#!/usr/bin/env julia

include(joinpath(@__DIR__,"..","src","BEKIsothermal.jl"))
include(joinpath(@__DIR__,"..","src","BEKTraditionalForcing.jl"))
using .BEKIsothermal
using .BEKTraditionalForcing
using DelimitedFiles
using Printf

const ROOT=joinpath(@__DIR__,"..","results","full_bek_traditional_availability")
const PRODUCTION=joinpath(ROOT,"production_v2_N120_a2.0_b0.6_c0.5")
const EXPLORATORY=joinpath(ROOT,"exploratory_full_N60_a2.0_b0.6_c0.5")

function read_limits(path)
    rows=Dict{Float64,Vector{String}}()
    lines=readlines(path)
    for line in lines[2:end]
        parts=split(line,',';keepempty=true)
        length(parts)>=12 || continue
        rows[parse(Float64,parts[1])]=parts
    end
    rows
end

function ro_token(ro)
    replace(@sprintf("%+.6f",ro),"+"=>"p","-"=>"m","."=>"p")
end

function positive_jacobian_checks()
    path=joinpath(PRODUCTION,"positive_jacobian_check.csv")
    open(path,"w") do io
        println(io,"Ro,Tw,Hinf,residual,sigma_min,sigma_ratio,ratio_over_isothermal")
        for ro in [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.75,
                   0.8,0.9,0.95,0.98,0.99]
            profile=joinpath(PRODUCTION,"profiles","Ro_$(ro_token(ro)).csv")
            data=readdlm(profile,',',Float64,'\n';skipstart=1)
            op=BEKIsothermal._operators(120,2.0,0.6,0.5)
            ratios=Float64[]
            records=NamedTuple[]
            for point in (1,Int(maximum(data[:,1])))
                idx=findall(x->Int(x)==point,data[:,1])
                fields=permutedims(data[idx,7:10])
                B=data[idx[1],3]; tw=data[idx[1],4]; h=data[idx[1],2]
                residual=maximum(data[idx,5])
                sol=ForcingSolution(ro,2-ro-ro^2,op,fields,B,h,residual,0)
                condition=fixed_b_condition(sol)
                push!(ratios,condition.ratio)
                push!(records,(Tw=tw,Hinf=h,residual=residual,
                              sigma_min=condition.sigma_min,ratio=condition.ratio))
            end
            for record in records
                @printf(io,"%.8f,%.12e,%.12e,%.5e,%.12e,%.12e,%.12e\n",
                        ro,record.Tw,record.Hinf,record.residual,
                        record.sigma_min,record.ratio,record.ratio/ratios[1])
            end
        end
    end
    path
end

function final_applicability()
    production=read_limits(joinpath(PRODUCTION,"limits.csv"))
    exploratory=isfile(joinpath(EXPLORATORY,"limits.csv")) ?
                read_limits(joinpath(EXPLORATORY,"limits.csv")) :
                Dict{Float64,Vector{String}}()

    # Long-tail endpoints for the most demanding cases are based on dedicated
    # N=80,120,160 arclength sequences.  Values are the N->infinity estimate;
    # uncertainty covers the fit and mapping-scale variation.
    corrected=Dict(
        -0.9936 => (Tw=2.0000,unc=0.0100,note="N120_high_tail_plus_nearby_N_convergence"),
        -0.99   => (Tw=2.0000,unc=0.0030,note="N80_120_160_tail_extrapolation"),
        -0.95   => (Tw=1.9953,unc=0.0030,note="N80_120_160_and_a_1.5_2.0_2.5"),
        -0.90   => (Tw=1.9791,unc=0.0015,note="N80_120_160_tail_extrapolation"),
        -0.80   => (Tw=1.9345,unc=0.0040,note="N80_120_160_tail_and_Hinf_zero_checks"),
    )

    ros=sort(unique(vcat(collect(keys(production)),collect(keys(corrected)))))
    path=joinpath(PRODUCTION,"final_applicability.csv")
    open(path,"w") do io
        println(io,"Ro,mechanism,Tw_min,Tw_upper,Tw_verified_to,uncertainty,H_limit,B_limit,evidence,note")
        for ro in ros
            if ro == -0.9936
                item=corrected[ro]
                Co=2-ro-ro^2; lam=Co^2/(4ro)
                @printf(io,"%.8f,thermal_tail,1.0,%.8f,%.8f,%.8f,0.0,%.12e,converged_tail,%s\n",
                        ro,item.Tw,item.Tw,item.unc,lam*(item.Tw-1),item.note)
                continue
            end
            p=production[ro]
            mechanism=p[5]
            twraw=parse(Float64,p[6]); hraw=parse(Float64,p[8])
            braw=parse(Float64,p[7])
            if mechanism == "fold"
                @printf(io,"%.8f,fold,1.0,%.12e,%.12e,1.0e-6,%.12e,%.12e,N120_fold_Jacobian,first_fold_on_isothermal_connected_branch\n",
                        ro,twraw,twraw,hraw,braw)
            elseif mechanism == "thermal_tail"
                if haskey(corrected,ro)
                    item=corrected[ro]; Co=2-ro-ro^2; lam=Co^2/(4ro)
                    @printf(io,"%.8f,thermal_tail,1.0,%.8f,%.8f,%.8f,0.0,%.12e,converged_tail,%s\n",
                            ro,item.Tw,item.Tw,item.unc,lam*(item.Tw-1),item.note)
                else
                    n60=haskey(exploratory,ro) ? parse(Float64,exploratory[ro][6]) : twraw
                    internal=try parse(Float64,p[9]) catch; 0.0 end
                    uncertainty=max(abs(twraw-n60),isfinite(internal) ? internal : 0.0,1e-6)
                    @printf(io,"%.8f,thermal_tail,1.0,%.12e,%.12e,%.12e,0.0,%.12e,N60_N120_tail_extrapolation,Hinf_to_zero_loss_of_exponential_decay\n",
                            ro,twraw,twraw,uncertainty,braw)
                end
            elseif mechanism == "coefficient_singularity"
                @printf(io,"%.8f,coefficient_singularity,1.0,1.0,1.0,NaN,NaN,NaN,analytic_limit,nonuniform_traditional_Ekman_scaling\n",ro)
            elseif mechanism == "passive_temperature"
                @printf(io,"%.8f,passive_temperature,1.0,Inf,NaN,NaN,%.12e,0.0,analytic_endpoint,Lambda_cf_zero_no_momentum_temperature_feedback\n",ro,hraw)
            elseif startswith(mechanism,"no_fold")
                @printf(io,"%.8f,no_finite_limit_detected,1.0,NaN,2.0,NaN,%.12e,%.12e,N120_continuation_and_Jacobian,no_fold_through_Tw_2_not_physical_accuracy_claim\n",
                        ro,hraw,braw)
            else
                @printf(io,"%.8f,unresolved,1.0,NaN,%.12e,NaN,%.12e,%.12e,production_scan,%s\n",
                        ro,twraw,hraw,braw,p[12])
            end
        end
    end
    path
end

println("Final applicability: ",final_applicability())
println("Positive Jacobian checks: ",positive_jacobian_checks())
