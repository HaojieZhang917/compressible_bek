" Designed by hj_zhang on Nov 16,2024,TJU "
module FDqScheme

    using Polynomials, GenericLinearAlgebra, SparseArrays, LinearAlgebra

    struct LagrangePolynomial
        order::Int64
        p_0
    end

    function construct_lagrange_polynomial(x)
        p_0 = Array{Polynomial, 1}(undef, length(x))
        for i in 1:length(x)
            p_0[i] = Polynomial([-x[i], 1])
        end
        lagrange_polynomial=LagrangePolynomial(length(x)-1, p_0)
        return lagrange_polynomial
    end

    function polyval(poly::LagrangePolynomial, x)
        return prod([poly.p_0[i](x) for i in 1:length(poly.p_0)])
    end

    #返回x0到x1中的极值
    function find_extreme(poly::LagrangePolynomial, x0, x1)

        xrange=range(x0, x1, length=11)
        idx = 0
        while (xrange[end]-xrange[begin])>eps(1.0e0)*20.0e0
            y=[polyval(poly, x) for x in xrange]
            value, idx=findmax(abs.(y))
            xrange=range(xrange[max(idx-1, 1)], xrange[min(idx+1, 11)], length=11)
        end
        return xrange[idx]
    end


    function lagrange_pi(poly::LagrangePolynomial, x, xi)
        # lagrange_pi=LagrangePolynomial(poly.order, poly.p_0)
        lagrange_pi = LagrangePolynomial(poly.order, [poly.p_0[i] for i in 1:poly.order+1])
        lagrange_pi.p_0[findall(x->x==xi, x)].=Polynomial([1.0e0])
        return lagrange_pi
    end

    function derivative_lagrange_pi(poly::LagrangePolynomial, x)

        lagrange_pi=LagrangePolynomial(poly.order, [poly.p_0[i] for i in 1:poly.order+1])
        derivative_lagrange_pi = 0.0e0
        for i in 1:poly.order+1
            derivative_langrange_pi_local = 1.0e0
            for j in 1:poly.order+1
                if j!=i
                    derivative_langrange_pi_local = derivative_langrange_pi_local*poly.p_0[j](x)
                    # @show j, poly.p_0[j](x)
                else
                    derivative_langrange_pi_local = derivative_langrange_pi_local*derivative(poly.p_0[j])(x)
                    # @show j, derivative(poly.p_0[j])(x)
                end
            end
            # @show derivative_langrange_pi_local
            derivative_lagrange_pi=derivative_lagrange_pi + derivative_langrange_pi_local
        end
        return derivative_lagrange_pi
    end

    function LagrangeScheme(x, x0)
        pi0=construct_lagrange_polynomial(x)
        iloc = findall(x->x==x0, x)[1]
        LagrangeScheme=[]
        for i in 1:length(x)
            push!(LagrangeScheme, 
                derivative_lagrange_pi(lagrange_pi(pi0, x, x[i]), x[iloc])/polyval(lagrange_pi(pi0, x, x[i]), x[i]))
        end
        return LagrangeScheme
    end

    #construct q order interpolation of degree q-1
    function construct_si(N, q)
        N = N - 1
        si=[]
        if isodd(q)
            for i in 1:Int32((q-1)/2)
                push!(si, 0)
            end
            for i in 0:N-q
                push!(si, i)
            end
            for i in 1:Int32((q+1)/2)
                push!(si, N-q)
            end
        else
            for i in 1:Int32(q/2)
                push!(si, 0)
            end
            for i in 0:N-q
                push!(si, i)
            end
            for i in 1:Int32(q/2)
                push!(si, N-q)
            end
        end
        return reshape(Int64.(si), N+1, :)
    end

    # return a matrix for a distributed nodes in a section
    function construct_schemes(x::Vector, q)
        n = length(x)
        si = construct_si(n, q) .+ 1
        scheme_matrix = zeros(n, n)
        for i in 1:n
            scheme_matrix[i, si[i]:si[i]+q]=LagrangeScheme(x[si[i]:si[i]+q], x[i])
        end
        return sparse(scheme_matrix)
    end

    #construct_polynomial factor pi_y
    function construct_polynomial(y, N, q, si)
        pi_y = Array{LagrangePolynomial}(undef, N-1)
        for i in 1:N-1
            pi_y[i] = construct_lagrange_polynomial(y[si[i]+1:si[i]+q+1])
        end
        return pi_y
    end

    # search_extrema_xi
    function search_extrema_xi(y, pi_y, N)
        x=Array{Float64}(undef, N+1)
        x[1] = -1#*Float64(N)
        x[end] = 1#*Float64(N)
        for i in 1:N-1
            x[i+1] = find_extreme(pi_y[i], y[i], y[i+1])
        end
        return x
    end

    #get_f
    function get_f(y, N, q, si, x, flg=true)
        #construct_polynomial factor pi_y
            pi_y = construct_polynomial(y, N, q-1, si)
        # search extrema xi
        if flg
            x_extrema = search_extrema_xi(y, pi_y, N)
        else
            x_extrema = x
        end
        #return abs(pi(x)) at xi
        abs_pi_x = []

        push!(abs_pi_x, abs(polyval(pi_y[1], x_extrema[1])))
        for i in 1:N-1
            push!(abs_pi_x, abs(polyval(pi_y[i], x_extrema[i+1])))
        end
        push!(abs_pi_x, abs(polyval(pi_y[N-1], x_extrema[N+1])))

        f= []
        for i in 1:N
            push!(f, abs_pi_x[i+1]-abs_pi_x[i])
        end
        return Float64.(f), Float64.(x_extrema), pi_y
    end

    function construct_J(y, x, f, N, q, si)
        epsilon = sqrt(eps(1.0e0))
        ymin = 1.0e-6
        J = zeros(N, N)

        for i in 1:N
            h = zeros(N)
            if abs(y[i])>ymin
                h[i] = epsilon*y[i]
            else
                h[i] = ymin*sign(y[i]+eps(1.0e0))
            end
            y_new=y+h
            f_new, = get_f(y_new, N, q, si, x, false)
            J[:, i] = (f_new-f)/h[i]
        end
        return sparse(J)
    end


    function FD_q_Scheme(n, q_target)

        N = n - 1

        # guess y_i using q=N, equal space
        y = Array{Float64}(undef, N)
        for i in 1:N
            y[i] = (2.0e0/Float64(2N)-1.0e0+2.0e0/Float64(N)*Float64(i-1))
        end

        x=zeros(N+1)

        for q in 2:1:q_target
            si = construct_si(N, q-1)

            f, x, pi_y= get_f(y, N, q, si, x)

            maxdy = 1.0e0
            while maxdy > 2.2e-6
                J = construct_J(y, x, f, N, q, si)
                dy = - J\f

                if(y+dy!=sort(vec(y+dy))) 
                    @show J
                    @show J[1, 1], J[end, end]
                    sleep()
                end
                y = y+ dy
                f, x= get_f(y, N, q, si, x)
                maxdy = maximum(abs.(dy))
            end
        end

        return x, construct_schemes(collect(x), q_target)
    end

end
module CRD_BF

    export sol_baseflowODE,velocity,Cheb,phi_var,f_q,T_var

    using LinearAlgebra
    using BoundaryValueDiffEq
    using DifferentialEquations
    using BSplineKit
    using IterativeSolvers

 function In_Su(Num)
        T0 = 273
        S = 114
        sigma = 0.7
        Tw = 0.5
        function ODE!(du, u , p, t)
            
            U = u[1]
            V = u[2]
            W = u[3]
            T = u[4]
            dU = u[5]
            dV = u[6]
            dT = u[7]
            du[1] = u[5]
            du[2] = u[6]
            du[3] = -2 * u[1]
            du[4] = u[7]
            du[5] = ((u[4])^(1/2)*((T0+S)/(u[4].*273+S)))^(-1)*(u[3]*u[5]+u[1]^2-(u[2]+1)^2)
            du[6] = ((u[4])^(1/2)*((T0+S)/(u[4].*273+S)))^(-1)*(u[3]*u[6]+2*u[1]*(u[2]+1))
            du[7] = ((u[4])^(1/2)*((T0+S)/(u[4].*273+S)))^(-1)*sigma*(u[3]*u[7])

        end
            function bc!(residual, u, p, t)
                residual[1] = u[begin][1] - 0
                residual[2] = u[end][1] - 0
                residual[3] = u[begin][2] - 0
                residual[4] = u[end][2] + 1
                residual[5] = u[begin][3] - 0
                residual[6] = u[begin][4] - Tw
                residual[7] = u[end][4] - 1
            end
            tspan = (0,20)
        ini = [0,0,0,0.5,0.5103341351120374,-0.6151547026271073,0]
        prob = BVProblem(ODE!, bc! ,ini ,tspan)
        sol = solve(prob, Shooting(Vern7()))
        t=range(0.0, 20, Num)
        u=sol(t)
        return u , t
        end
 function sol_baseflowODE(tspan,Num)

        function oneDiskODE!(du, u , p, t)
            
            U = u[1]
            dU = u[2]
            V = u[3]
            dV = u[4]
            W = u[5]
            du[1] = dU
            ddU = U^2 + W*dU - (V+1.0e0)^2
            du[2] = ddU
            ddV = 2.0e0*U*(V + 1.0e0) + W*dV
            du[3] = dV
            du[4] = ddV                          
            du[5] = -2.0e0*U

        end
        function oneDiskbc!(residual, u , p, t)

            residual[1] = u[begin][1] 
            residual[2] = u[begin][3] 
            residual[3] = u[begin][5] 
            residual[4] = u[end][1] 

            residual[5] = u[end][3] + 1.0e0
        end
            prob = BVProblem(oneDiskODE!, oneDiskbc!, [0.0, 0.5103341351120374, 0.0, -0.6151547026271073, 0.0] ,tspan, dtmax=0.01)
            sol = solve(prob, Shooting(Vern7()), dt=0.01)
            t=range(0.0, 20, Num)
            u=sol(t)
        
        return u , t

        end
 function velocity(u,t,phi)

    U = u[1 , :]
    dU = u[2 , :]
    V = u[3 , :]
    dV = u[4 , :]
    W = u[5 , :]
    F_U = itp = interpolate(t, U , BSplineOrder(4))
    F_dU = itp = interpolate(t, dU , BSplineOrder(4))
    F_dV = itp = interpolate(t, dV , BSplineOrder(4))
    F_W= itp = interpolate(t, W , BSplineOrder(4))
    F_phi = interpolate(t, phi , BSplineOrder(4))
    return U,dU,V,dV,W,F_U,F_dU,F_dV,F_W,F_phi

  end

 function Cheb(u,v,w,phi,t,N)
        θ = range(0,length=N+1,stop=pi)
        x = reshape(-cos.(θ), N+1, 1)
        c = [2; ones(N-1, 1) ; 2] .* (-1) .^ (0:N)
        X = repeat(x, 1, N+1);
        dX = X - X';
        D = (c * (1 ./ c)') ./ (dX .+ I(N+1));
        D = D - diagm(vec(sum(D, dims=2))); 
        # for i=1:N+1
        #     x[i] = 10 * x[i] .+ 10
        # end
        # D = 0.1 * D
        # D2 = 0.01 * D^2
        for i=1:N+1
            D[i,:]=D[i,:].*((2*x[i]^3-x[i]^2+3*x[i]-4)^2/(20*(6*x[i]^2-2*x[i]+3)))
        end
        for i=1:N+1
            x[i]=(4*x[i]^3-2*x[i]^2+6*x[i]+12)/(-2*x[i]^3+x[i]^2-3*x[i]+4)
            if x[i]>20
                x[i]=20
            end
        end
        D2 = D^2
        itpw = itp = BSplineKit.interpolate(t, w , BSplineOrder(4))
        itpu = itp = BSplineKit.interpolate(t, u , BSplineOrder(4))
        itpv = itp = BSplineKit.interpolate(t, v , BSplineOrder(4))
        itpphi = itp = BSplineKit.interpolate(t, phi , BSplineOrder(4))
        u = zeros(N+1,1)
        v = zeros(N+1,1)
        w = zeros(N+1,1)
        phi = zeros(N+1,1)
        for i=1:N+1
            u[i,1] = itpu(x[i])
            v[i,1] = itpv(x[i])
            w[i,1] = itpw(x[i])
            phi[i,1] = itpphi(x[i])
        end
        return D,D2,x,u,v,w,phi
        end


 function Physical_Interpretation(T,delt,Num)
        z = zeros(Num,1)
        integral_T = delt * T
        for i = 1 : 1 : Num
            z[i,1] = sum(integral_T[1:i])
        end
        return z
     end
    end
 function phi_var(u,del,N)
    phi = zeros(N)
    for i = 1 : N
        z = 0
        for j = 1 : i
            z = z + u[j,1]*del
        end
        phi[i] = z
    end
    return phi
  end
 function f_q(sigma,F_du,F_dv,F_u,F_phi,tspan,t)
    function ODE_f!(du,u,p,t)
        q = u[1]
        dq = u[2]
        du[1] = dq
        du[2] = 2*sigma*(F_du(t)^2+F_dv(t)^2+F_u(t)*u[1]-F_phi(t)*u[2])
    end
    function BC_f!(residual, u, p, t)
        residual[1] = u[begin][1]
        residual[2] = u[end][1]
    end     
    prob = DifferentialEquations.BVProblem(ODE_f!, BC_f!, [0.0, 0.0], tspan)
    sol = solve(prob, Shooting(Vern7()), dt=0.01)
    f = sol(t)
    function ODE_q!(du,u,p,t)
        q = u[1]
        dq = u[2]
        du[1] = dq
        du[2] = -2*sigma*F_phi(t)*u[2]
    end
    function BC_q!(residual, u, p, t)
        residual[1] = u[begin][1] - 1
        residual[2] = u[end][1]
    end     
    prob = DifferentialEquations.BVProblem(ODE_q!, BC_q!, [1,0], tspan)
    sol1 = solve(prob,MIRK4(), dt=0.01)
    q = sol1(t)
    f = f[1,:]
    q = q[1,:]
    return f,q
  end
 function T_var(Mx,f,q,Tw,gamma)
    T = 1 .- ( (gamma-1)/2 )*Mx^2 * f + (Tw - 1) * q
    return T
  end

import .CRD_BF
import .FDqScheme
using DifferentialEquations
using BSplineKit
using LinearAlgebra
function baseflow_var(N_cheb)

    N = 1001
    tspan = (0,20)
    t = range(0,20,N)
    sigma = 0.72
    u,z = CRD_BF.sol_baseflowODE(tspan,N)
    u0 = u[1,:]
    PHI = phi_var(u0,t.step.hi,N)
    u0,du0,v0,dv0,w0,F_u,F_du,F_dv,F_w,F_phi = CRD_BF.velocity(u,t,PHI)
    f,q = f_q(sigma,F_du,F_dv,F_u,F_phi,tspan,t)
    D,D2,x,F,G,W,phi = CRD_BF.Cheb(u0,v0,w0,PHI,t,N_cheb)
    # itpw = itp = BSplineKit.interpolate(t, w0 , BSplineOrder(4))
    # itpu = itp = BSplineKit.interpolate(t, u0 , BSplineOrder(4))
    # itpv = itp = BSplineKit.interpolate(t, v0 , BSplineOrder(4))
    # itpphi = itp = BSplineKit.interpolate(t, PHI , BSplineOrder(4))

    # F = zeros(N_cheb+1,1)
    # G = zeros(N_cheb+1,1)
    # W = zeros(N_cheb+1,1)
    # phi = zeros(N_cheb+1,1)

    # for i=1:N_cheb
    
    #     F[i,1] = itpu(x[i])
    #     G[i,1] = itpv(x[i])
    #     W[i,1] = itpw(x[i])
    #     phi[i,1] = itpphi(x[i])
    
    # end

    return F,G,W,f,q,D,D2,x,t,phi
 end

function T_ca(R,Ma,f,q,W,gamma,Tw,x,t,N)
    Mx = Ma * R
    T = T_var(Mx,f,q,Tw,gamma)
    itpT = BSplineKit.interpolate(t,T,BSplineOrder(4)) 
    T = zeros(N+1,1)
    for i=1:N+1
        T[i,1] = itpT(x[i])
    end
    RHO = 1 ./ T
    H = W .* T
    return RHO,H,T
 end
function Time_mode(sigma,gamma,Re,Ma,Tw,al,be)
    F,G,W,f,q,D,x,t,PHI = baseflow_var(N_cheb)
    rho,H,T = T_ca(Re,Ma,f,q,W,gamma,Tw,x,t,N_cheb)
    lam = -(2/3) * T
    kappa = (1/sigma) * T
    # phi_hat = [u,v,w,rho,p,t]

    A11 = ((im * al * R + 1) .* rho) .* I(N_cheb + 1)
    A12 = (im * be * R).*rho .* I(N_cheb + 1)
    A13 = R * rho .*(D*rho.* I(N_cheb + 1) + rho .* D)
    A14 = (im * R .* (al * F + be * G)) .* I(N_cheb + 1)  + 2 * F .* I(N_cheb + 1)  + rho .* (D * H) .* I(N_cheb + 1) + rho .* H .* D 
    A15 = zeros(N_cheb + 1, N_cheb + 1)
    A16 = zeros(N_cheb + 1, N_cheb + 1)

    A21 = im * R * rho .* (al *  F + be * G) .* I(N_cheb + 1) + rho .* F .* I(N_cheb + 1) + al^2 .* (lam + 2  * T) .* I(N_cheb + 1) + be^2 * T .* I(N_cheb + 1) + rho.^2 .* H .* D - rho .* D2
    A22 = -2 * rho .* (G .+ 1) .* I(N_cheb + 1) + al * be * (lam .+ T) .* I(N_cheb + 1)
    A23 = R * rho.^2 .* (D*F) .* I(N_cheb + 1) - im * al *(rho .* (D*T) .* I(N_cheb + 1) + (1 .+ lam .* rho) .* D)
    A24 = rho .* (D2 * F) .* I(N_cheb + 1)
    A25 = im * al * R .* I(N_cheb + 1)
    A26 = -rho .* D * (rho .* (D*F)) .* I(N_cheb + 1) - rho.^2 .* (D*F) .* D

    A31 = 2 * rho .* (G .+ 1) .* I(N_cheb + 1) + al * be * (lam .+ T) .* I(N_cheb + 1)
    A32 = im * R * rho .* (al *  F + be * G) .* I(N_cheb + 1) + rho .* F .* I(N_cheb + 1) + be^2 * (lam + 2 * T) .* I(N_cheb + 1) + al^2 * T .* I(N_cheb + 1) + rho.^2 .* H .* D - rho .* D2
    A33 = -im * be * (rho .* (D*T) .* I(N_cheb + 1) + (1 .+ lam .* rho) .* D) + R * rho.^2 .* (D*G) .* I(N_cheb + 1)
    A34 = 2 * F .* (G .+ 1) .* I(N_cheb + 1) + rho .* H .* (D*G) .* I(N_cheb + 1)
    A35 = im * be * R .* I(N_cheb + 1)
    A36 = -rho .* D * (rho .* (D*G)) .* I(N_cheb + 1) - rho.^2 .* (D*G) .* D

    A41 = -im * al *(rho .* (D*lam).* I(N_cheb + 1) + (1 .+ lam .* rho).*D)
    A42 = -im * be * (rho .* (D*lam) .* I(N_cheb + 1) + (1 .+ lam .* rho) .* D)
    A43 = im * R * rho .* (al * F + be * G) .* I(N_cheb + 1) + rho.^2 .* (D * H .* I(N_cheb + 1) + H .* D) 
        - rho .* (2 .+ lam .* rho) .* D2 + (al^2 + be^2) .* T .* I(N_cheb + 1)
    A44 = zeros(N_cheb + 1, N_cheb + 1)
    A45 = R * rho .* D 
    A46 = - im * rho .* D * (al * F + be * G) .* I(N_cheb + 1)

    A51 = -2 * (gamma - 1) * Ma^2 * rho .* (D*F) .* D 
    A52 = -2 * (gamma - 1) * Ma^2 * rho .* (D*G) .* D
    A53 = -2 * im * (gamma - 1) * Ma^2 * D * (al * F + be * G) .* I(N_cheb + 1) + R * rho.^2 .* D * T .* I(N_cheb + 1)
    A54 = rho .* H .* (D * T) .* I(N_cheb + 1)
    A55 = -(gamma - 1) * Ma^2 .* (im * R * (al * F + be * G) .* I(N_cheb + 1) + rho .* H .* D)
    A56 = im * R * rho .* (al * F + be * G) .* I(N_cheb + 1) + (al^2 + be^2) * kappa .* I(N_cheb + 1) + rho.^2 .* (H .* D - kappa .* D2)
        + (1/sigma) * (-rho .* ((D*rho .* D*T + rho .* D2 * T ).* I(N_cheb + 1) - (rho .* (D * T) .* D)))
        + (-(gamma - 1) * Ma^2 * rho.^2 .* ((D * F).^2 + (D * G).^2 )).* I(N_cheb + 1)

    A61 = zeros(N_cheb + 1, N_cheb + 1)
    A62 = zeros(N_cheb + 1, N_cheb + 1)
    A63 = zeros(N_cheb + 1, N_cheb + 1)
    A64 = - T .* I(N_cheb + 1)
    A65 = gamma * Ma^2 .* I(N_cheb + 1)
    A66 = -rho .* I(N_cheb + 1)

    A = [A11 A12 A13 A14 A15 A16; A21 A22 A23 A24 A25 A26; A31 A32 A33 A34 A35 A36; A41 A42 A43 A44 A45 A46; A51 A52 A53 A54 A55 A56; A61 A62 A63 A64 A65 A66]

    B11 = zeros(N_cheb + 1, N_cheb + 1)
    B12 = zeros(N_cheb + 1, N_cheb + 1)
    B13 = zeros(N_cheb + 1, N_cheb + 1)
    B14 = im * R .* I(N_cheb + 1)
    B15 = zeros(N_cheb + 1, N_cheb + 1)
    B16 = zeros(N_cheb + 1, N_cheb + 1)

    B21 = im * R * rho .* I(N_cheb + 1)
    B22 = zeros(N_cheb + 1, N_cheb + 1)
    B23 = zeros(N_cheb + 1, N_cheb + 1)
    B24 = zeros(N_cheb + 1, N_cheb + 1)
    B25 = zeros(N_cheb + 1, N_cheb + 1)
    B26 = zeros(N_cheb + 1, N_cheb + 1)

    B31 = zeros(N_cheb + 1, N_cheb + 1)
    B32 = im * R * rho .* I(N_cheb + 1)
    B33 = zeros(N_cheb + 1, N_cheb + 1)
    B34 = zeros(N_cheb + 1, N_cheb + 1)
    B35 = zeros(N_cheb + 1, N_cheb + 1)
    B36 = zeros(N_cheb + 1, N_cheb + 1)

    B41 = zeros(N_cheb + 1, N_cheb + 1)
    B42 = zeros(N_cheb + 1, N_cheb + 1)
    B43 = im * R * rho .* I(N_cheb + 1)
    B44 = zeros(N_cheb + 1, N_cheb + 1)
    B45 = zeros(N_cheb + 1, N_cheb + 1)
    B46 = zeros(N_cheb + 1, N_cheb + 1)

    B51 = zeros(N_cheb + 1, N_cheb + 1)
    B52 = zeros(N_cheb + 1, N_cheb + 1)
    B53 = zeros(N_cheb + 1, N_cheb + 1)
    B54 = zeros(N_cheb + 1, N_cheb + 1)
    B55 = -(gamma-1) * Ma^2 * im * R  .* I(N_cheb + 1)
    B56 = im * R * rho .* I(N_cheb + 1)

    B61 = zeros(N_cheb + 1, N_cheb + 1)
    B62 = zeros(N_cheb + 1, N_cheb + 1)
    B63 = zeros(N_cheb + 1, N_cheb + 1)
    B64 = zeros(N_cheb + 1, N_cheb + 1)
    B65 = zeros(N_cheb + 1, N_cheb + 1)
    B66 = zeros(N_cheb + 1, N_cheb + 1)

    B = [B11 B12 B13 B14 B15 B16; B21 B22 B23 B24 B25 B26; B31 B32 B33 B34 B35 B36; B41 B42 B43 B44 B45 B46; B51 B52 B53 B54 B55 B56; B61 B62 B63 B64 B65 B66]

    A = A[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,5N_cheb + 5,5N_cheb + 6,6N_cheb + 6)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,5N_cheb + 5,5N_cheb + 6,6N_cheb + 6))]
    B = B[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,5N_cheb + 5,5N_cheb + 6,6N_cheb + 6)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,5N_cheb + 5,5N_cheb + 6,6N_cheb + 6))]


    return A,B
 end
function Spatial_mode(F,G,H,rho,lam,kappa,T,sigma,gamma,R,Ma,omega,be)

    # (A0 +A1*alpha +A2*alpha^2)ϕ=0 
    # phi_hat = [u,v,w,rho,p,t]
    A0_11 = rho .* I(N_cheb + 1)
    A0_12 = im * be * R * rho .* I(N_cheb + 1)
    A0_13 = R * rho .* (D*rho .* I(N_cheb + 1) + rho .* D)
    A0_14 = im * R * (be * G .- omega ) .* I(N_cheb + 1)  + 2 * F .* I(N_cheb + 1)  + rho .* (D * H) .* I(N_cheb + 1) + rho .* H .* D 
    A0_15 = zeros(N_cheb + 1, N_cheb + 1)
    A0_16 = zeros(N_cheb + 1, N_cheb + 1)

    A0_21 = im * R * rho .* (be * G .- omega ) .* I(N_cheb + 1) + rho .* F .* I(N_cheb + 1) + be^2 * T .* I(N_cheb + 1) + rho.^2 .* H .* D - rho .* D2
    A0_22 = -2 * rho .* (G .+ 1) .* I(N_cheb + 1)
    A0_23 = R * rho.^2 .* D*F .* I(N_cheb + 1)
    A0_24 = rho .* (D2 * F) .* I(N_cheb + 1)
    A0_25 = zeros(N_cheb + 1, N_cheb + 1)
    A0_26 = -rho .* (D * (rho .* (D * F))) .* I(N_cheb + 1) - rho.^2 .* (D*F) .* D

    A0_41 = zeros(N_cheb + 1, N_cheb + 1)
    A0_42 = -im * be * (rho .* (D*lam) + (1 .+ lam .* rho)).*D
    A0_43 = im * R * rho .* (be * G .- omega) .* I(N_cheb + 1) + rho.^2 .* ((D*H) .* I(N_cheb + 1) + H .* D) - rho .* (2 .+ lam .* rho) .* D + be^2 * T .* I(N_cheb + 1)
    A0_44 = zeros(N_cheb + 1, N_cheb + 1)
    A0_45 = R * rho .* D
    A0_46 = -im * rho .* (be .* (D*G)) .* I(N_cheb + 1)

    A0_31 = 2 * rho .* (G .+ 1) .* I(N_cheb + 1)
    A0_32 = im * R * rho .* (be * G .- omega) .* I(N_cheb + 1) + rho .* F .* I(N_cheb + 1) + be^2 * (lam + 2 * T) .* I(N_cheb + 1) + rho.^2 .* H .* D - rho .* D
    A0_33 = R * rho.^2  .* (D*G) .* I(N_cheb + 1) - im * be * (rho .* (D*lam) .* I(N_cheb + 1) + (1 .+ lam .* rho) .* D)
    A0_34 = 2 * F .* (G .+ 1) .* I(N_cheb + 1) + rho .* H .* (D*G) .* I(N_cheb + 1)
    A0_35 = im * be * R .* I(N_cheb + 1)
    A0_36 = -rho .* (D * (rho .* (D * G))) .* I(N_cheb + 1) - rho.^2 .* (D*G) .* D

    A0_51 = -2 * (gamma - 1) * Ma^2 * rho .* (D*F) .* D
    A0_52 = -2 * (gamma - 1) * Ma^2 * rho .* (D*G) .* D
    A0_53 = -2 * im * (gamma - 1) * Ma^2 * (be * (D*G))  .* I(N_cheb + 1 ) + R * rho.^2 .* (D*T) .* I(N_cheb + 1)
    A0_54 = rho .* H .* (D*T) .* I(N_cheb + 1)
    A0_55 = -(gamma - 1) * Ma^2 * (im * R * (be * G .- omega) .* I(N_cheb + 1) + rho .* H .* D)
    A0_56 = im * R * rho .* (be * G .- omega) .* I(N_cheb + 1)+ be^2 * kappa .* I(N_cheb + 1) + rho.^2 .* (H .* D - kappa .* D2) + (1/sigma) * (-rho .* ((D*rho) .* T .* I(N_cheb + 1) + rho .* (D2 * T) .* I(N_cheb + 1) - rho .* (D*T) .* D))
        + (-(gamma - 1) * Ma^2 * rho.^2 .* ((D*F).^2 + (D*G).^2) .* I(N_cheb + 1))

    A0_61 = zeros(N_cheb + 1, N_cheb + 1)
    A0_62 = zeros(N_cheb + 1, N_cheb + 1)
    A0_63 = zeros(N_cheb + 1, N_cheb + 1)
    A0_64 = -T .* I(N_cheb + 1)
    A0_65 = gamma * Ma^2 .* I(N_cheb + 1)
    A0_66 = -rho .* I(N_cheb + 1)

    A1_11 = im * R * rho .* I(N_cheb + 1)
    A1_12 = zeros(N_cheb + 1, N_cheb + 1)
    A1_13 = zeros(N_cheb + 1, N_cheb + 1)
    A1_14 = im * R * F .* I(N_cheb + 1)
    A1_15 = zeros(N_cheb + 1, N_cheb + 1)
    A1_16 = zeros(N_cheb + 1, N_cheb + 1)

    A1_21 = im * R * rho .* F .* I(N_cheb + 1)
    A1_22 = be * (lam .+ T) .* I(N_cheb + 1)
    A1_23 = -im * ((rho .* (D*T) + (1 .+ lam .* rho)).* D)
    A1_24 = zeros(N_cheb + 1, N_cheb + 1)
    A1_25 = im * R .* I(N_cheb + 1)
    A1_26 = zeros(N_cheb + 1, N_cheb + 1)

    A1_31 = be * (lam + T) .* I(N_cheb + 1)
    A1_32 = im * R * rho .* F .* I(N_cheb + 1)
    A1_33 = zeros(N_cheb + 1, N_cheb + 1)
    A1_34 = zeros(N_cheb + 1, N_cheb + 1)
    A1_35 = zeros(N_cheb + 1, N_cheb + 1)
    A1_36 = zeros(N_cheb + 1, N_cheb + 1)

    A1_41 = -im * (rho .* (D*lam) .* I(N_cheb + 1) + (1 .+ lam .* rho) .* D)
    A1_42 = zeros(N_cheb + 1, N_cheb + 1)
    A1_43 = im * R * rho .* F .* I(N_cheb + 1)
    A1_44 = zeros(N_cheb + 1, N_cheb + 1)
    A1_45 = zeros(N_cheb + 1, N_cheb + 1)
    A1_46 = zeros(N_cheb + 1, N_cheb + 1)

    A1_51 = zeros(N_cheb + 1, N_cheb + 1)
    A1_52 = zeros(N_cheb + 1, N_cheb + 1)
    A1_53 = zeros(N_cheb + 1, N_cheb + 1)
    A1_54 = zeros(N_cheb + 1, N_cheb + 1)
    A1_55 = -(gamma - 1) * Ma^2 * im * R * F .* I(N_cheb + 1)
    A1_56 = im * R * rho .* F .* I(N_cheb + 1)

    A1_61 = zeros(N_cheb + 1, N_cheb + 1)
    A1_62 = zeros(N_cheb + 1, N_cheb + 1)
    A1_63 = zeros(N_cheb + 1, N_cheb + 1)
    A1_64 = zeros(N_cheb + 1, N_cheb + 1)
    A1_65 = zeros(N_cheb + 1, N_cheb + 1)
    A1_66 = zeros(N_cheb + 1, N_cheb + 1)

    A2_11 = zeros(N_cheb + 1, N_cheb + 1)
    A2_12 = zeros(N_cheb + 1, N_cheb + 1)
    A2_13 = zeros(N_cheb + 1, N_cheb + 1)
    A2_14 = zeros(N_cheb + 1, N_cheb + 1)
    A2_15 = zeros(N_cheb + 1, N_cheb + 1)
    A2_16 = zeros(N_cheb + 1, N_cheb + 1)

    A2_21 = (lam + 2 * T) .* I(N_cheb + 1)
    A2_22 = zeros(N_cheb + 1, N_cheb + 1)
    A2_23 = zeros(N_cheb + 1, N_cheb + 1)
    A2_24 = zeros(N_cheb + 1, N_cheb + 1)
    A2_25 = zeros(N_cheb + 1, N_cheb + 1)
    A2_26 = zeros(N_cheb + 1, N_cheb + 1)

    A2_31 = zeros(N_cheb + 1, N_cheb + 1)
    A2_32 = T .* I(N_cheb + 1)
    A2_33 = zeros(N_cheb + 1, N_cheb + 1)
    A2_34 = zeros(N_cheb + 1, N_cheb + 1)
    A2_35 = zeros(N_cheb + 1, N_cheb + 1)
    A2_36 = zeros(N_cheb + 1, N_cheb + 1)

    A2_41 = zeros(N_cheb + 1, N_cheb + 1)
    A2_42 = zeros(N_cheb + 1, N_cheb + 1)
    A2_43 = T .* I(N_cheb + 1)
    A2_44 = zeros(N_cheb + 1, N_cheb + 1)
    A2_45 = zeros(N_cheb + 1, N_cheb + 1)
    A2_46 = zeros(N_cheb + 1, N_cheb + 1)

    A2_51 = zeros(N_cheb + 1, N_cheb + 1)
    A2_52 = zeros(N_cheb + 1, N_cheb + 1)
    A2_53 = zeros(N_cheb + 1, N_cheb + 1)
    A2_54 = zeros(N_cheb + 1, N_cheb + 1)
    A2_55 = zeros(N_cheb + 1, N_cheb + 1)
    A2_56 = kappa .* I(N_cheb + 1)

    A2_61 = zeros(N_cheb + 1, N_cheb + 1)
    A2_62 = zeros(N_cheb + 1, N_cheb + 1)
    A2_63 = zeros(N_cheb + 1, N_cheb + 1)
    A2_64 = zeros(N_cheb + 1, N_cheb + 1)
    A2_65 = zeros(N_cheb + 1, N_cheb + 1)
    A2_66 = zeros(N_cheb + 1, N_cheb + 1) 

    A0 = [A0_11 A0_12 A0_13 A0_14 A0_15 A0_16 ; A0_21 A0_22 A0_23 A0_24 A0_25 A0_26 ; A0_31 A0_32 A0_33 A0_34 A0_35 A0_36 ; A0_41 A0_42 A0_43 A0_44 A0_45 A0_46 ; A0_51 A0_52 A0_53 A0_54 A0_55 A0_56 ; A0_61 A0_62 A0_63 A0_64 A0_65 A0_66]
    A1 = [A1_11 A1_12 A1_13 A1_14 A1_15 A1_16 ; A1_21 A1_22 A1_23 A1_24 A1_25 A1_26 ; A1_31 A1_32 A1_33 A1_34 A1_35 A1_36 ; A1_41 A1_42 A1_43 A1_44 A1_45 A1_46 ; A1_51 A1_52 A1_53 A1_54 A1_55 A1_56 ; A1_61 A1_62 A1_63 A1_64 A1_65 A1_66]
    A2 = [A2_11 A2_12 A2_13 A2_14 A2_15 A2_16 ; A2_21 A2_22 A2_23 A2_24 A2_25 A2_26 ; A2_31 A2_32 A2_33 A2_34 A2_35 A2_36 ; A2_41 A2_42 A2_43 A2_44 A2_45 A2_46 ; A2_51 A2_52 A2_53 A2_54 A2_55 A2_56 ; A2_61 A2_62 A2_63 A2_64 A2_65 A2_66]

    A0 = A0[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,5N_cheb + 5,5N_cheb + 6,6N_cheb + 6)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,5N_cheb + 5,5N_cheb + 6,6N_cheb + 6))]
    A1 = A1[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,5N_cheb + 5,5N_cheb + 6,6N_cheb + 6)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,5N_cheb + 5,5N_cheb + 6,6N_cheb + 6))]
    A2 = A2[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,5N_cheb + 5,5N_cheb + 6,6N_cheb + 6)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,5N_cheb + 5,5N_cheb + 6,6N_cheb + 6))]

    return A0,A1,A2
end