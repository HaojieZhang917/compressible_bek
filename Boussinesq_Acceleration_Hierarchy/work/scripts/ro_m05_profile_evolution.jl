#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKTraditionalForcing.jl"))
using .BEKIsothermal
using .BEKTraditionalForcing
using Printf

const RO = -0.5
const TW_VALUES = collect(1.0:0.1:1.6)
const N = 120
const A = 2.0
const MAP_B = 0.6
const C = 0.5
const TOL = 1e-10
const OUTDIR = joinpath(@__DIR__, "..", "results", "ro_m05_profile_evolution",
                        "N$(N)_a$(A)_b$(MAP_B)_c$(C)")

function solve_sweep()
    states = ForcingSolution[]
    current = solve_forcing_isothermal(RO; degree=N, a=A, b=MAP_B, c=C,
                                        tolerance=TOL)
    push!(states, current)
    lam = lambda_cf(RO)
    for tw in TW_VALUES[2:end]
        current = solve_forcing_fixed_b(lam * (tw - 1), current; tolerance=TOL)
        push!(states, current)
    end
    states
end

function write_results(states)
    mkpath(OUTDIR)
    open(joinpath(OUTDIR, "summary.csv"), "w") do io
        println(io, "Tw,B,Hinf,ell_T,residual,iterations")
        for (tw, s) in zip(TW_VALUES, states)
            @printf(io, "%.12e,%.12e,%.12e,%.12e,%.5e,%d\n",
                    tw, s.B, s.Hinf, thermal_tail_length(s), s.residual, s.iterations)
        end
    end
    open(joinpath(OUTDIR, "profiles.csv"), "w") do io
        println(io, "Tw,eta,H,F,G,Theta,T")
        for (tw, s) in zip(TW_VALUES, states)
            for j in eachindex(s.operators.eta)
                eta = s.operators.eta[j]
                theta = s.fields[4,j]
                temp = 1 + (tw - 1) * theta
                @printf(io, "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                        tw, eta, s.fields[1,j], s.fields[2,j], s.fields[3,j], theta, temp)
            end
        end
    end
end

xml_escape(s) = replace(string(s), "&"=>"&amp;", "<"=>"&lt;", ">"=>"&gt;")

function polyline_points(x, y, xmin, xmax, ymin, ymax, left, top, width, height)
    points = String[]
    for (xx, yy) in zip(x, y)
        isfinite(xx) && isfinite(yy) || continue
        xmin <= xx <= xmax || continue
        px = left + width * (xx - xmin) / (xmax - xmin)
        py = top + height * (1 - (yy - ymin) / (ymax - ymin))
        push!(points, @sprintf("%.2f,%.2f", px, py))
    end
    join(points, " ")
end

function panel(io, states; xfield=:eta, yfield, title, xlabel="eta", ylabel,
               xmin=0.0, xmax=12.0, ymin=nothing, ymax=nothing,
               left, top, width=390.0, height=245.0)
    series = Vector{Tuple{Vector{Float64},Vector{Float64}}}()
    for (tw, s) in zip(TW_VALUES, states)
        eta = s.operators.eta
        y = if yfield == :H
            vec(s.fields[1,:])
        elseif yfield == :F
            vec(s.fields[2,:])
        elseif yfield == :G
            vec(s.fields[3,:])
        elseif yfield == :Theta
            vec(s.fields[4,:])
        elseif yfield == :T
            1 .+ (tw - 1) .* vec(s.fields[4,:])
        elseif yfield == :Hinf
            TW_VALUES
        elseif yfield == :ellT
            TW_VALUES
        else
            error("unknown yfield")
        end
        x = yfield in (:Hinf, :ellT) ? TW_VALUES : eta
        yy = if yfield == :Hinf
            [q.Hinf for q in states]
        elseif yfield == :ellT
            [thermal_tail_length(q) for q in states]
        else
            y
        end
        push!(series, (Float64.(x), Float64.(yy)))
        yfield in (:Hinf, :ellT) && break
    end
    visible = Float64[]
    for (x,y) in series, (xx,yy) in zip(x,y)
        isfinite(xx) && isfinite(yy) && xmin <= xx <= xmax && push!(visible,yy)
    end
    lo = ymin === nothing ? minimum(visible) : Float64(ymin)
    hi = ymax === nothing ? maximum(visible) : Float64(ymax)
    pad = 0.06 * max(hi-lo, 1e-6)
    ymin2 = ymin === nothing ? lo-pad : Float64(ymin)
    ymax2 = ymax === nothing ? hi+pad : Float64(ymax)

    println(io, "<g>")
    @printf(io, "<text x=\"%.1f\" y=\"%.1f\" class=\"title\">%s</text>\n",
            left+width/2, top-18, xml_escape(title))
    @printf(io, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" class=\"frame\"/>\n",
            left,top,width,height)
    for k in 0:4
        xx = left + width*k/4
        value = xmin + (xmax-xmin)*k/4
        @printf(io,"<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" class=\"grid\"/>\n",xx,top,xx,top+height)
        @printf(io,"<text x=\"%.1f\" y=\"%.1f\" class=\"tick\">%.2g</text>\n",xx,top+height+20,value)
        yy = top + height*(1-k/4)
        yvalue = ymin2 + (ymax2-ymin2)*k/4
        @printf(io,"<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" class=\"grid\"/>\n",left,yy,left+width,yy)
        @printf(io,"<text x=\"%.1f\" y=\"%.1f\" class=\"ytick\">%.3g</text>\n",left-9,yy+4,yvalue)
    end
    @printf(io,"<text x=\"%.1f\" y=\"%.1f\" class=\"axis\">%s</text>\n",left+width/2,top+height+44,xml_escape(xlabel))
    @printf(io,"<text x=\"%.1f\" y=\"%.1f\" transform=\"rotate(-90 %.1f %.1f)\" class=\"axis\">%s</text>\n",
            left-58,top+height/2,left-58,top+height/2,xml_escape(ylabel))
    colors = ["#1f77b4","#17becf","#2ca02c","#bcbd22","#ff7f0e","#d62728","#9467bd"]
    for (i,(x,y)) in enumerate(series)
        pts = polyline_points(x,y,xmin,xmax,ymin2,ymax2,left,top,width,height)
        @printf(io,"<polyline points=\"%s\" fill=\"none\" stroke=\"%s\" stroke-width=\"2.1\"/>\n",pts,colors[i])
        if yfield in (:Hinf,:ellT)
            for (xx,yy) in zip(x,y)
                px=left+width*(xx-xmin)/(xmax-xmin); py=top+height*(1-(yy-ymin2)/(ymax2-ymin2))
                @printf(io,"<circle cx=\"%.2f\" cy=\"%.2f\" r=\"3.2\" fill=\"%s\"/>\n",px,py,colors[1])
            end
        end
    end
    println(io,"</g>")
end

function write_svg(states)
    ellmax = maximum(thermal_tail_length.(states))
    thermal_xmax = max(16.0, ceil(6ellmax/5)*5)
    path = joinpath(OUTDIR, "ro_m05_profile_evolution.svg")
    open(path,"w") do io
        println(io,"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1400\" height=\"1050\" viewBox=\"0 0 1400 1050\">")
        println(io,"<style>.bg{fill:#fff}.frame{fill:none;stroke:#222;stroke-width:1}.grid{stroke:#d9d9d9;stroke-width:1}.title{font-family:sans-serif;font-size:18px;font-weight:bold;text-anchor:middle;fill:#111}.axis{font-family:sans-serif;font-size:15px;text-anchor:middle;fill:#111}.tick{font-family:sans-serif;font-size:12px;text-anchor:middle;fill:#333}.ytick{font-family:sans-serif;font-size:12px;text-anchor:end;fill:#333}.legend{font-family:sans-serif;font-size:14px;fill:#111}.main{font-family:sans-serif;font-size:24px;font-weight:bold;text-anchor:middle;fill:#111}.note{font-family:sans-serif;font-size:14px;fill:#333}</style>")
        println(io,"<rect class=\"bg\" width=\"1400\" height=\"1050\"/>")
        println(io,"<text x=\"700\" y=\"34\" class=\"main\">Ro=-0.5: traditional-Boussinesq BEK base-flow evolution</text>")
        panel(io,states;yfield=:F,title="(a) Radial velocity F",ylabel="F",left=95.0,top=90.0)
        panel(io,states;yfield=:G,title="(b) Azimuthal velocity G",ylabel="G",left=540.0,top=90.0)
        panel(io,states;yfield=:H,title="(c) Axial velocity H",ylabel="H",left=985.0,top=90.0)
        panel(io,states;yfield=:T,title="(d) Temperature T",ylabel="T",xmax=thermal_xmax,ymin=1.0,ymax=1.6,left=95.0,top=440.0)
        panel(io,states;yfield=:Hinf,title="(e) Far-field axial velocity",xlabel="Tw",ylabel="Hinf",xmin=1.0,xmax=1.6,left=540.0,top=440.0)
        panel(io,states;yfield=:ellT,title="(f) Thermal-tail length",xlabel="Tw",ylabel="ell_T",xmin=1.0,xmax=1.6,left=985.0,top=440.0)
        colors = ["#1f77b4","#17becf","#2ca02c","#bcbd22","#ff7f0e","#d62728","#9467bd"]
        println(io,"<g transform=\"translate(195 815)\">")
        for (i,tw) in enumerate(TW_VALUES)
            x=155*(i-1)
            @printf(io,"<line x1=\"%d\" y1=\"0\" x2=\"%d\" y2=\"0\" stroke=\"%s\" stroke-width=\"3\"/><text x=\"%d\" y=\"5\" dx=\"34\" class=\"legend\">Tw=%.1f</text>\n",x,x+26,colors[i],x,tw)
        end
        println(io,"</g>")
        @printf(io,"<text x=\"95\" y=\"875\" class=\"note\">N=%d, a=%.1f, b=%.1f, c=%.1f, Pr=0.72; exact fixed-Tw solutions. Thermal panel extends to eta=%.0f.</text>\n",N,A,MAP_B,C,thermal_xmax)
        println(io,"</svg>")
    end
    path
end

function main()
    states = solve_sweep()
    write_results(states)
    figure = write_svg(states)
    for (tw,s) in zip(TW_VALUES,states)
        @printf("Tw=%.1f B=% .8f Hinf=% .10f ellT=%.8f residual=%.3e\n",
                tw,s.B,s.Hinf,thermal_tail_length(s),s.residual)
    end
    println("figure=",figure)
end

main()
