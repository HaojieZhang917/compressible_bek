" 
    ∊ Designed by hj_zhang on Nov 16,2024,TJU 
    ∊ Add Stability Analysis on Jan 16,2025,TJU
                                                "
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
        return (scheme_matrix)
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

    export sol_baseflowODE,velocity,Cheb,phi_var,f_q,T_var,Physical_Interpretation
    using LinearAlgebra
    using BSplineKit
    using DifferentialEquations
    using IterativeSolvers
    using PyCall
    
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
        t=range(0.0, 30, Num)
        u=sol(t)
        return u , t
        end
 function sol_baseflowODE(Ro)
        py"""
        import numpy as np
        from scipy.integrate import solve_bvp
        import matplotlib.pyplot as plt
        kappa = $Ro
        def oneDiskODE(z, y):
        
                # Y0 = H, Y1 = F', Y2 = F, Y3 = G', Y4 = G
                dydz = np.zeros((5, len(z)))
                dydz = np.array([-
                                2.0*
                                y[2], kappa * 
                                (y[2] *
                                y[2] +
                                    y[0] *
                                    y[1] -
                                    (y[4] *
                                    y[4]- 
                                    1.0)) -
                                (2.0 -
                                    kappa - 
                                    kappa**2) *
                                (y[4] -
                                    1.0), y[1], kappa *
                                (2.0 * 
                                    y[2] *
                                    y[4] +
                                    y[0] *
                                    y[3]) +
                                (2.0 -
                                    kappa -
                                    kappa**2) *
                                y[2], y[3]])
                return dydz 
        
        def oneDiskBC(ya, yb):
                resa = np.array([ya[0],
                                ya[2],
                                ya[4]])
                
                resb = np.array([yb[2],
                                yb[4] - 1.0])
                
                return np.concatenate((resa, resb))
        
        
        z = np.linspace(0, 30, 20000)
        y = np.zeros((5, len(z)))
        y_guess = np.zeros((5, z.size))
        if kappa == 1:
                y_guess[0] = 1.2
                y_guess[1] = 0
                y_guess[2] = 0
                y_guess[3] = 0
                y_guess[4] = 1
        elif kappa == -1:
                y_guess[0] = 1.2
                y_guess[1] = 0
                y_guess[2] = 0
                y_guess[3] = 0
                y_guess[4] = 1
        else:
                y_guess[0] = 1.2
                y_guess[1] = 0
                y_guess[2] = 0
                y_guess[3] = 0
                y_guess[4] = 1
        
        
        
        
        solution = solve_bvp(oneDiskODE, oneDiskBC, z, y_guess,tol=1e-10,max_nodes=5000000)
        
        x_plot = np.linspace(0, 30, 20000)
        
        
        y1_plot = solution.sol(x_plot)[0]
        y2_plot = solution.sol(x_plot)[2]
        y3_plot = solution.sol(x_plot)[4]
        y4_plot = solution.sol(x_plot)[1]
        y5_plot = solution.sol(x_plot)[3]
        
        """
        w0 = py"y1_plot"
        u0 = py"y2_plot"
        v0 = py"y3_plot"
        du0 = py"y4_plot"
        dv0 = py"y5_plot"
        x = py"x_plot"
        
        return -u0,-v0,-w0,-du0,-dv0,x

        end
 function velocity(u0,v0,w0,du0,dv0,t,phi)
    U = u0
    V = v0
    W = w0
    dU = du0
    dV = dv0
    F_U = itp = interpolate(t, U , BSplineOrder(4))
    F_dU = itp = interpolate(t, dU , BSplineOrder(4))
    F_dV = itp = interpolate(t, dV , BSplineOrder(4))
    F_W= itp = interpolate(t, W , BSplineOrder(4))
    F_phi = interpolate(t, phi , BSplineOrder(4))
    return U,dU,V,dV,W,F_U,F_dU,F_dV,F_W,F_phi

  end

 function Cheb(N,mode)
        θ = range(0,length=N+1,stop=pi)
        x = reshape(-cos.(θ), N+1, 1)
        c = [2; ones(N-1, 1) ; 2] .* (-1) .^ (0:N)
        X = repeat(x, 1, N+1);
        dX = X - X';
        D = (c * (1 ./ c)') ./ (dX .+ I(N+1));
        D = D - diagm(vec(sum(D, dims=2)));
    return D,x
 end


 function Physical_Interpretation(T,delt,Num)
        z = zeros(Num)
        z[1] = 0
        integral_T = delt * T
        for i = 2 : 1 : Num
                z[i] = sum(integral_T[1:i])
        end
        z = vec(z)
        return z
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
end
import .CRD_BF
using BSplineKit
using LinearAlgebra
function baseflow_var(N_cheb,Ro,mode)

    N = 20000
    tspan = (0,30)
    t = range(0,30,N)
    sigma = 0.72
    u0,v0,w0,du0,dv0,x = CRD_BF.sol_baseflowODE(Ro)
    PHI = CRD_BF.phi_var(u0,t.step.hi,N)
    u0,du0,v0,dv0,w0,F_u,F_du,F_dv,F_w,F_phi = CRD_BF.velocity(u0,v0,w0,du0,dv0,t,PHI)
    if mode == "cheb"
        D,x = CRD_BF.Cheb(N_cheb,mode)
    elseif mode =="fdq"
        x,D = x,D = FDqScheme.FD_q_Scheme(N_cheb+1, 20)
    else
        error("no method")
    end
    if Ro == -1
        a = 2
        b = 0.6
        c = 0.5
    else
        a = 1
        b = 0.6
        c = 0.5
    end
    for i=1:N_cheb + 1
        D[i,:]=D[i,:].* (1-b*x[i]-(1-b)*(x[i]^3+c*(1-x[i]^2)))^2/(2a*(b .+ 3 * (1-b)*x[i]^2 - 2 * c * (1-b) * x[i]))
    end
    for i=1:N_cheb + 1
        x[i] = a * (1+b*x[i]+(1-b)*(x[i]^3+c*(1-x[i]^2)))/(1-b*x[i]-(1-b)*(x[i]^3+c*(1-x[i]^2)))
        if x[i]>30
            x[i]=30
        end
    end
    
    D2 = D^2;

    f,q = CRD_BF.f_q(sigma,F_du,F_dv,F_u,F_phi,tspan,t)
    
    u0 = -1 * u0
    v0 = -1 * v0
    w0 = -1 * w0
    
    return u0,v0,w0,f,q,D,D2,x
 end

function T_ca(Mr,f,q,W,gamma,Tw)
    T = CRD_BF.T_var(Mr,f,q,Tw,gamma)
    H = W .* T
    return H,T

 end
function interp(u,v,w,T,x,N,mode)
    if mode == "sim"
        z = range(0,30,20000)
        itu = BSplineKit.interpolate(z, u , BSplineOrder(4))
        itv = BSplineKit.interpolate(z, v , BSplineOrder(4))
        itw = BSplineKit.interpolate(z, w , BSplineOrder(4))
        itT = BSplineKit.interpolate(z, T , BSplineOrder(4))
        F = zeros(N+1,1)
        G = zeros(N+1,1)
        H = zeros(N+1,1)
        T = zeros(N+1,1)
        for i = 1 : N + 1
            F[i,1] = itu(x[i])
            G[i,1] = itv(x[i])
            H[i,1] = itw(x[i])
            T[i,1] = itT(x[i])
        end
        rho = 1 ./ T
    end
    if mode == "phy"
        z = CRD_BF.Physical_Interpretation(T,0.02,1001)
        itu = BSplineKit.interpolate(z, u , BSplineOrder(4))
        itv = BSplineKit.interpolate(z, v , BSplineOrder(4))
        itw = BSplineKit.interpolate(z, w , BSplineOrder(4))
        itT = BSplineKit.interpolate(z, T , BSplineOrder(4))
        F = zeros(N+1,1)
        G = zeros(N+1,1)
        H = zeros(N+1,1)
        T = zeros(N+1,1)
        for i = 1 : N + 1
            F[i,1] = itu(x[i])
            G[i,1] = itv(x[i])
            H[i,1] = itw(x[i])
            T[i,1] = itT(x[i])
        end
        rho = 1 ./ T
    end
        return F,G,H,T,rho,z
  end

function Timemode(F,G,H,rho,lam,kappa,T,sigma,gamma,R,Ma,al,be,N_cheb,Ro,Co,D,D2)
        A0_11 = rho .* I(N_cheb + 1) + al * im * R * rho .* I(N_cheb + 1)
        A0_12 = im * be * R * rho .* I(N_cheb + 1)
        A0_13 = R * rho .* (D*rho .* I(N_cheb + 1) + rho .* D)
        A0_14 = im * R * (be * G ) .* I(N_cheb + 1)  + 2 * F .* I(N_cheb + 1)  + rho .* (D * H) .* I(N_cheb + 1) + rho .* H .* D + al * im * R * F .* I(N_cheb + 1)
        A0_15 = zeros(N_cheb + 1, N_cheb + 1)
    
        A0_21 = im * R * rho .* (be * G ) .* I(N_cheb + 1) + Ro * rho .*  F .* I(N_cheb + 1) + be^2 * T .* I(N_cheb + 1) + Ro * rho.^2 .* H .* D - rho .* D2 + al*(im * R * rho .* F .* I(N_cheb + 1)) + al^2 * ((lam + 2 * T) .* I(N_cheb + 1))
        A0_22 = -1 * rho .* (2*Ro * G .+ Co) .* I(N_cheb + 1) + al *(be * (lam .+ T) .* I(N_cheb + 1))
        A0_23 = R * rho.^2 .* D*F .* I(N_cheb + 1) + al * (-im * ((rho .* (D*T) + (1 .+ lam .* rho)).* D))
        A0_24 = rho .* (D2 * F) .* I(N_cheb + 1) + al * (im * R .* I(N_cheb + 1) * (gamma*Ma^2)^(-1) * T .* I(N_cheb + 1))
        A0_25 = -rho .* (D * rho .* D * F + rho .* (D2 * F)) .* I(N_cheb + 1) - rho.^2 .* (D*F) .* D + al * (im * R .* I(N_cheb + 1) * (gamma*Ma^2)^(-1) * rho .* I(N_cheb + 1))
    
        A0_31 =  rho .* (2 * Ro * G .+ Co) .* I(N_cheb + 1) + al * (be * (lam + T) .* I(N_cheb + 1))
        A0_32 = im * R * rho .* (be * G ) .* I(N_cheb + 1) + Ro * rho .* F .* I(N_cheb + 1) + be^2 * (lam + 2 * T) .* I(N_cheb + 1) + Ro * rho.^2 .* H .* D - rho .* D2 + al * (im * R * rho .* F .* I(N_cheb + 1)) + al^2 * (T .* I(N_cheb + 1))
        A0_33 = R * rho.^2  .* (D*G) .* I(N_cheb + 1) - im * be * (rho .* (D*T) .* I(N_cheb + 1) + (1 .+ lam .* rho) .* D)
        A0_34 = F .* (2 * Ro * G .+ Co) .* I(N_cheb + 1) + Ro * rho .* H .* (D*G) .* I(N_cheb + 1) + im * be * R * (gamma*Ma^2)^(-1) * T .* I(N_cheb + 1)
        A0_35 = -rho .* (D * rho .* D * G + rho .* (D2 * G)) .* I(N_cheb + 1) - rho.^2 .* (D*G) .* D  + im * be * R * (gamma*Ma^2)^(-1) * rho .* I(N_cheb + 1)
    
        A0_41 = zeros(N_cheb + 1, N_cheb + 1) + al * (-im * (rho .* (D*lam) .* I(N_cheb + 1) + (1 .+ lam .* rho) .* D))
        A0_42 = -im * be * (rho .* (D*lam) + (1 .+ lam .* rho)).*D
        A0_43 = im * R * rho .* (be * G) .* I(N_cheb + 1) + Ro * rho.^2 .*  ((D*H) .* I(N_cheb + 1) + H .* D) - rho .* (2 .+ lam .* rho) .* D2 + be^2 * T .* I(N_cheb + 1) + al * (im * R * rho .* F .* I(N_cheb + 1)) + al^2 * (T .* I(N_cheb + 1))
        A0_44 = (gamma*Ma^2)^(-1) * R * rho .* ((D*T) .* I(N_cheb + 1) + T .* D)
        A0_45 = -im * rho .* (be .* (D*G)) .* I(N_cheb + 1) + (gamma*Ma^2)^(-1) * R * rho .* (D*rho.* I(N_cheb + 1) + rho .* D) + al *  (-im * rho .* ( (D*F)) .* I(N_cheb + 1))
    
        A0_51 = -2 * (gamma - 1) * Ma^2 * rho .* (D*F) .* D
        A0_52 = -2 * (gamma - 1) * Ma^2 * rho .* (D*G) .* D
        A0_53 = -2 * im * (gamma - 1) * Ma^2 * (be * (D*G))  .* I(N_cheb + 1 ) + R * rho.^2 .* (D*T) .* I(N_cheb + 1) + al * ( -2 * im * (gamma - 1) * Ma^2 * D * F .* I(N_cheb + 1))
        A0_54 = rho .* H .* (D*T) .* I(N_cheb + 1) - ((gamma - 1) *  (gamma)^(-1) * ( (im * R * (be * G ) .* T .* I(N_cheb + 1) )+ rho .* H .* D*T .* I(N_cheb + 1) + rho .* H .* T .* D)) + al * (-(gamma - 1) * (gamma)^(-1) *  im * R * F  .* T .* I(N_cheb + 1))
        A0_55 = im * R * rho .* (be * G) .* I(N_cheb + 1)+ be^2 * kappa .* I(N_cheb + 1) + rho.^2 .* (H .* D - kappa .* D2) + (1/sigma) * (-rho .* ((D*rho) .* (D*T) .* I(N_cheb + 1) + rho .* (D2 * T) .* I(N_cheb + 1) - rho .* (D*T) .* D))+ (-(gamma - 1) * Ma^2 * rho.^2 .* ((D*F).^2 + (D*G).^2) .* I(N_cheb + 1)) - ((gamma - 1) * (gamma)^(-1) * ((im * R * (be * G) .* rho .* I(N_cheb + 1)) + rho .* H .* D*rho .* I(N_cheb + 1) + rho .* H .* rho .* D))
        + al * (im * R * rho .* F .* I(N_cheb + 1) - (gamma - 1) * (gamma)^(-1) * im * R * F .* rho .* I(N_cheb + 1)) + al^2 * (kappa .* I(N_cheb + 1))
        A0_omega_coeff_11 = zeros(N_cheb + 1, N_cheb + 1)
        A0_omega_coeff_12 = zeros(N_cheb + 1, N_cheb + 1)
        A0_omega_coeff_13 = zeros(N_cheb + 1, N_cheb + 1)
        A0_omega_coeff_14 = im * R .* I(N_cheb + 1)
        A0_omega_coeff_15 = zeros(N_cheb + 1, N_cheb + 1)
    
        A0_omega_coeff_21 = im * R * rho .* I(N_cheb + 1)
        A0_omega_coeff_22 = zeros(N_cheb + 1, N_cheb + 1)
        A0_omega_coeff_23 = zeros(N_cheb + 1, N_cheb + 1)
        A0_omega_coeff_24 = zeros(N_cheb + 1, N_cheb + 1)
        A0_omega_coeff_25 = zeros(N_cheb + 1, N_cheb + 1)
    
        A0_omega_coeff_31 = zeros(N_cheb + 1, N_cheb + 1)
        A0_omega_coeff_32 = im * R * rho .* I(N_cheb + 1)
        A0_omega_coeff_33 = zeros(N_cheb + 1, N_cheb + 1)
        A0_omega_coeff_34 = zeros(N_cheb + 1, N_cheb + 1)
        A0_omega_coeff_35 = zeros(N_cheb + 1, N_cheb + 1)
    
        A0_omega_coeff_41 = zeros(N_cheb + 1, N_cheb + 1)
        A0_omega_coeff_42 = zeros(N_cheb + 1, N_cheb + 1)
        A0_omega_coeff_43 = im * R * rho .* I(N_cheb + 1)
        A0_omega_coeff_44 = zeros(N_cheb + 1, N_cheb + 1)
        A0_omega_coeff_45 = zeros(N_cheb + 1, N_cheb + 1)
    
        A0_omega_coeff_51 = zeros(N_cheb + 1, N_cheb + 1)
        A0_omega_coeff_52 = zeros(N_cheb + 1, N_cheb + 1)
        A0_omega_coeff_53 = zeros(N_cheb + 1, N_cheb + 1)
        A0_omega_coeff_54 = ((gamma - 1) * (gamma)^(-1) * im * R .* T .* I(N_cheb + 1))
        A0_omega_coeff_55 = im * R * rho .* I(N_cheb + 1) + ((gamma - 1) * (gamma)^(-1) * im * R .* rho .* I(N_cheb + 1))
        A0 = [A0_11 A0_12 A0_13 A0_14 A0_15 ; A0_21 A0_22 A0_23 A0_24 A0_25 ; A0_31 A0_32 A0_33 A0_34 A0_35 ; A0_41 A0_42 A0_43 A0_44 A0_45 ; A0_51 A0_52 A0_53 A0_54 A0_55  ]
        A0_omega_coeff = [A0_omega_coeff_11 A0_omega_coeff_12 A0_omega_coeff_13 A0_omega_coeff_14 A0_omega_coeff_15;
        A0_omega_coeff_21 A0_omega_coeff_22 A0_omega_coeff_23 A0_omega_coeff_24 A0_omega_coeff_25;
        A0_omega_coeff_31 A0_omega_coeff_32 A0_omega_coeff_33 A0_omega_coeff_34 A0_omega_coeff_35;
        A0_omega_coeff_41 A0_omega_coeff_42 A0_omega_coeff_43 A0_omega_coeff_44 A0_omega_coeff_45;
        A0_omega_coeff_51 A0_omega_coeff_52 A0_omega_coeff_53 A0_omega_coeff_54 A0_omega_coeff_55]
        B0 = A0
        B1 = A0_omega_coeff
        B0 = B0[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5))]
        B1 = B1[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5))]
     return B0,B1
end
function Spatial_mode(F,G,H,rho,lam,kappa,T,sigma,gamma,R,Ma,omega,be,N_cheb,D,D2)
    Ro = 1
    Co = 2
    # (A0 +A1*alpha +A2*alpha^2)ϕ=0 
    # phi_hat = [u,v,w,rho,p,t]
    A0_11 = rho .* I(N_cheb + 1)
    A0_12 = im * be * R * rho .* I(N_cheb + 1)
    A0_13 = R * rho .* (D*rho .* I(N_cheb + 1) + rho .* D)
    A0_14 = im * R * (be * G .- omega ) .* I(N_cheb + 1)  + 2 * F .* I(N_cheb + 1)  + rho .* (D * H) .* I(N_cheb + 1) + rho .* H .* D 
    A0_15 = zeros(N_cheb + 1, N_cheb + 1)

    A0_21 = im * R * rho .* (be * G .- omega ) .* I(N_cheb + 1) + Ro * rho .* F .* I(N_cheb + 1) + be^2 * T .* I(N_cheb + 1) + Ro * rho.^2 .* H .* D - rho .* D2
    A0_22 = -1 * rho .* (2*Ro * G .+ Co) .* I(N_cheb + 1)
    A0_23 = Ro * R * rho.^2 .* D*F .* I(N_cheb + 1)
    A0_24 = rho .* (D2 * F) .* I(N_cheb + 1)
    A0_25 = -rho .* (D * rho .* D * F + rho .* (D2 * F)) .* I(N_cheb + 1) - rho.^2 .* (D*F) .* D

    A0_31 =  rho .* (2 * Ro * G .+ Co) .* I(N_cheb + 1)
    A0_32 = im * R * rho .* (be * G .- omega) .* I(N_cheb + 1) + Ro * rho .* F .* I(N_cheb + 1) + be^2 * (lam + 2 * T) .* I(N_cheb + 1) + Ro * rho.^2 .* H .* D - rho .* D2
    A0_33 = R * rho.^2  .* (D*G) .* I(N_cheb + 1) - im * be * (rho .* (D*T) .* I(N_cheb + 1) + (1 .+ lam .* rho) .* D)
    A0_34 = F .* (2 * Ro * G .+ Co) .* I(N_cheb + 1) + Ro * rho .* H .* (D*G) .* I(N_cheb + 1) + im * be * R * (gamma*Ma^2)^(-1) * T .* I(N_cheb + 1)
    A0_35 = -rho .* (D * rho .* D * G + rho .* (D2 * G)) .* I(N_cheb + 1) - rho.^2 .* (D*G) .* D  + im * be * R * (gamma*Ma^2)^(-1) * rho .* I(N_cheb + 1)

    A0_41 = zeros(N_cheb + 1, N_cheb + 1)
    A0_42 = -im * be * (rho .* (D*lam) + (1 .+ lam .* rho)).*D
    A0_43 = im * R * rho .* (be * G .- omega) .* I(N_cheb + 1) + Ro * rho.^2 .*  ((D*H) .* I(N_cheb + 1) + H .* D) - rho .* (2 .+ lam .* rho) .* D2 + be^2 * T .* I(N_cheb + 1)
    A0_44 = (gamma*Ma^2)^(-1) * R * rho .* ((D*T) .* I(N_cheb + 1) + T .* D)
    A0_45 = -im * rho .* (be .* (D*G)) .* I(N_cheb + 1) + (gamma*Ma^2)^(-1) * R * rho .* (D*rho.* I(N_cheb + 1) + rho .* D)

    A0_51 = -2 * (gamma - 1) * Ma^2 * rho .* (D*F) .* D
    A0_52 = -2 * (gamma - 1) * Ma^2 * rho .* (D*G) .* D
    A0_53 = -2 * im * (gamma - 1) * Ma^2 * (be * (D*G))  .* I(N_cheb + 1 ) + R * rho.^2 .* (D*T) .* I(N_cheb + 1)
    A0_54 = rho .* H .* (D*T) .* I(N_cheb + 1) - ((gamma - 1) *  (gamma)^(-1) * ( (im * R * (be * G .- omega) .* T .* I(N_cheb + 1) )+ rho .* H .* D*T .* I(N_cheb + 1) + rho .* H .* T .* D))
    A0_55 = im * R * rho .* (be * G .- omega) .* I(N_cheb + 1)+ be^2 * kappa .* I(N_cheb + 1) + rho.^2 .* (H .* D - kappa .* D2) + (1/sigma) * (-rho .* ((D*rho) .* (D*T) .* I(N_cheb + 1) + rho .* (D2 * T) .* I(N_cheb + 1) - rho .* (D*T) .* D))
        + (-(gamma - 1) * Ma^2 * rho.^2 .* ((D*F).^2 + (D*G).^2) .* I(N_cheb + 1)) - ((gamma - 1) * (gamma)^(-1) * ((im * R * (be * G .- omega) .* rho .* I(N_cheb + 1)) + rho .* H .* D*rho .* I(N_cheb + 1) + rho .* H .* rho .* D))

    A1_11 = im * R * rho .* I(N_cheb + 1)
    A1_12 = zeros(N_cheb + 1, N_cheb + 1)
    A1_13 = zeros(N_cheb + 1, N_cheb + 1)
    A1_14 = im * R * F .* I(N_cheb + 1)
    A1_15 = zeros(N_cheb + 1, N_cheb + 1)

    A1_21 = im * R * rho .* F .* I(N_cheb + 1)
    A1_22 = be * (lam .+ T) .* I(N_cheb + 1)
    A1_23 = -im * ((rho .* (D*T) + (1 .+ lam .* rho)).* D)
    A1_24 = im * R .* I(N_cheb + 1) * (gamma*Ma^2)^(-1) * T .* I(N_cheb + 1)
    A1_25 = im * R .* I(N_cheb + 1) * (gamma*Ma^2)^(-1) * rho .* I(N_cheb + 1)

    A1_31 = be * (lam + T) .* I(N_cheb + 1)
    A1_32 = im * R * rho .* F .* I(N_cheb + 1)
    A1_33 = zeros(N_cheb + 1, N_cheb + 1)
    A1_34 = zeros(N_cheb + 1, N_cheb + 1)
    A1_35 = zeros(N_cheb + 1, N_cheb + 1)

    A1_41 = -im * (rho .* (D*lam) .* I(N_cheb + 1) + (1 .+ lam .* rho) .* D)
    A1_42 = zeros(N_cheb + 1, N_cheb + 1)
    A1_43 = im * R * rho .* F .* I(N_cheb + 1)
    A1_44 = zeros(N_cheb + 1, N_cheb + 1)
    A1_45 = -im * rho .* ( (D*F)) .* I(N_cheb + 1)

    A1_51 = zeros(N_cheb + 1, N_cheb + 1)
    A1_52 = zeros(N_cheb + 1, N_cheb + 1)
    A1_53 = -2 * im * (gamma - 1) * Ma^2 * D * F .* I(N_cheb + 1)
    A1_54 = -(gamma - 1) * (gamma)^(-1) *  im * R * F  .* T .* I(N_cheb + 1)
    A1_55 = im * R * rho .* F .* I(N_cheb + 1) - (gamma - 1) * (gamma)^(-1) * im * R * F .* rho .* I(N_cheb + 1)

    A2_11 = zeros(N_cheb + 1, N_cheb + 1)
    A2_12 = zeros(N_cheb + 1, N_cheb + 1)
    A2_13 = zeros(N_cheb + 1, N_cheb + 1)
    A2_14 = zeros(N_cheb + 1, N_cheb + 1)
    A2_15 = zeros(N_cheb + 1, N_cheb + 1)

    A2_21 = (lam + 2 * T) .* I(N_cheb + 1)
    A2_22 = zeros(N_cheb + 1, N_cheb + 1)
    A2_23 = zeros(N_cheb + 1, N_cheb + 1)
    A2_24 = zeros(N_cheb + 1, N_cheb + 1)
    A2_25 = zeros(N_cheb + 1, N_cheb + 1)

    A2_31 = zeros(N_cheb + 1, N_cheb + 1)
    A2_32 = T .* I(N_cheb + 1)
    A2_33 = zeros(N_cheb + 1, N_cheb + 1)
    A2_34 = zeros(N_cheb + 1, N_cheb + 1)
    A2_35 = zeros(N_cheb + 1, N_cheb + 1)

    A2_41 = zeros(N_cheb + 1, N_cheb + 1)
    A2_42 = zeros(N_cheb + 1, N_cheb + 1)
    A2_43 = T .* I(N_cheb + 1)
    A2_44 = zeros(N_cheb + 1, N_cheb + 1)
    A2_45 = zeros(N_cheb + 1, N_cheb + 1)

    A2_51 = zeros(N_cheb + 1, N_cheb + 1)
    A2_52 = zeros(N_cheb + 1, N_cheb + 1)
    A2_53 = zeros(N_cheb + 1, N_cheb + 1)
    A2_54 = zeros(N_cheb + 1, N_cheb + 1)
    A2_55 = kappa .* I(N_cheb + 1)

    A0 = [A0_11 A0_12 A0_13 A0_14 A0_15 ; A0_21 A0_22 A0_23 A0_24 A0_25 ; A0_31 A0_32 A0_33 A0_34 A0_35 ; A0_41 A0_42 A0_43 A0_44 A0_45 ; A0_51 A0_52 A0_53 A0_54 A0_55  ]
    A1 = [A1_11 A1_12 A1_13 A1_14 A1_15 ; A1_21 A1_22 A1_23 A1_24 A1_25 ; A1_31 A1_32 A1_33 A1_34 A1_35 ; A1_41 A1_42 A1_43 A1_44 A1_45 ; A1_51 A1_52 A1_53 A1_54 A1_55  ]
    A2 = [A2_11 A2_12 A2_13 A2_14 A2_15 ; A2_21 A2_22 A2_23 A2_24 A2_25 ; A2_31 A2_32 A2_33 A2_34 A2_35 ; A2_41 A2_42 A2_43 A2_44 A2_45 ; A2_51 A2_52 A2_53 A2_54 A2_55  ]

    A0 = A0[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5))]
    A1 = A1[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5))]
    A2 = A2[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5))]

    return A0,A1,A2
  end
function Spatial_mode_BEK(F,G,H,rho,lam,kappa,T,sigma,gamma,R,Ma,omega,be,N_cheb,Ro,Co,D,D2)

    # (A0 +A1*alpha +A2*alpha^2)ϕ=0 
    # phi_hat = [u,v,w,rho,t]
    A0_11 = rho .* I(N_cheb + 1)
    A0_12 = im * be * R * rho .* I(N_cheb + 1)
    A0_13 = R * rho .* (D*rho .* I(N_cheb + 1) + rho .* D)
    A0_14 = im * R * (be * G .- omega ) .* I(N_cheb + 1)  + 2 * F .* I(N_cheb + 1)  + rho .* (D * H) .* I(N_cheb + 1) + rho .* H .* D 
    A0_15 = zeros(N_cheb + 1, N_cheb + 1)

    A0_21 = im * R * rho .* (be * G .- omega ) .* I(N_cheb + 1) + Ro * rho .*  F .* I(N_cheb + 1) + be^2 * T .* I(N_cheb + 1) + Ro * rho.^2 .* H .* D - rho .* D2
    A0_22 = -1 * rho .* (2*Ro * G .+ Co) .* I(N_cheb + 1)
    A0_23 = R * rho.^2 .* D*F .* I(N_cheb + 1)
    A0_24 = rho .* (D2 * F) .* I(N_cheb + 1)
    A0_25 = -rho .* (D * rho .* D * F + rho .* (D2 * F)) .* I(N_cheb + 1) - rho.^2 .* (D*F) .* D

    A0_31 =  rho .* (2 * Ro * G .+ Co) .* I(N_cheb + 1)
    A0_32 = im * R * rho .* (be * G .- omega) .* I(N_cheb + 1) + Ro * rho .* F .* I(N_cheb + 1) + be^2 * (lam + 2 * T) .* I(N_cheb + 1) + Ro * rho.^2 .* H .* D - rho .* D2
    A0_33 = R * rho.^2  .* (D*G) .* I(N_cheb + 1) - im * be * (rho .* (D*T) .* I(N_cheb + 1) + (1 .+ lam .* rho) .* D)
    A0_34 = F .* (2 * Ro * G .+ Co) .* I(N_cheb + 1) + Ro * rho .* H .* (D*G) .* I(N_cheb + 1) + im * be * R * (gamma*Ma^2)^(-1) * T .* I(N_cheb + 1)
    A0_35 = -rho .* (D * rho .* D * G + rho .* (D2 * G)) .* I(N_cheb + 1) - rho.^2 .* (D*G) .* D  + im * be * R * (gamma*Ma^2)^(-1) * rho .* I(N_cheb + 1)

    A0_41 = zeros(N_cheb + 1, N_cheb + 1)
    A0_42 = -im * be * (rho .* (D*lam) + (1 .+ lam .* rho)).*D
    A0_43 = im * R * rho .* (be * G .- omega) .* I(N_cheb + 1) + Ro * rho.^2 .*  ((D*H) .* I(N_cheb + 1) + H .* D) - rho .* (2 .+ lam .* rho) .* D2 + be^2 * T .* I(N_cheb + 1)
    A0_44 = (gamma*Ma^2)^(-1) * R * rho .* ((D*T) .* I(N_cheb + 1) + T .* D)
    A0_45 = -im * rho .* (be .* (D*G)) .* I(N_cheb + 1) + (gamma*Ma^2)^(-1) * R * rho .* (D*rho.* I(N_cheb + 1) + rho .* D)

    A0_51 = -2 * (gamma - 1) * Ma^2 * rho .* (D*F) .* D
    A0_52 = -2 * (gamma - 1) * Ma^2 * rho .* (D*G) .* D
    A0_53 = -2 * im * (gamma - 1) * Ma^2 * (be * (D*G))  .* I(N_cheb + 1 ) + R * rho.^2 .* (D*T) .* I(N_cheb + 1)
    A0_54 = rho .* H .* (D*T) .* I(N_cheb + 1) - ((gamma - 1) *  (gamma)^(-1) * ( (im * R * (be * G .- omega) .* T .* I(N_cheb + 1) )+ rho .* H .* D*T .* I(N_cheb + 1) + rho .* H .* T .* D))
    A0_55 = im * R * rho .* (be * G .- omega) .* I(N_cheb + 1)+ be^2 * kappa .* I(N_cheb + 1) + rho.^2 .* (H .* D - kappa .* D2) + (1/sigma) * (-rho .* ((D*rho) .* (D*T) .* I(N_cheb + 1) + rho .* (D2 * T) .* I(N_cheb + 1) - rho .* (D*T) .* D))+ (-(gamma - 1) * Ma^2 * rho.^2 .* ((D*F).^2 + (D*G).^2) .* I(N_cheb + 1)) - ((gamma - 1) * (gamma)^(-1) * ((im * R * (be * G .- omega) .* rho .* I(N_cheb + 1)) + rho .* H .* D*rho .* I(N_cheb + 1) + rho .* H .* rho .* D))

    A1_11 = im * R * rho .* I(N_cheb + 1)
    A1_12 = zeros(N_cheb + 1, N_cheb + 1)
    A1_13 = zeros(N_cheb + 1, N_cheb + 1)
    A1_14 = im * R * F .* I(N_cheb + 1)
    A1_15 = zeros(N_cheb + 1, N_cheb + 1)

    A1_21 = im * R * rho .* F .* I(N_cheb + 1)
    A1_22 = be * (lam .+ T) .* I(N_cheb + 1)
    A1_23 = -im * ((rho .* (D*T) + (1 .+ lam .* rho)).* D)
    A1_24 = im * R .* I(N_cheb + 1) * (gamma*Ma^2)^(-1) * T .* I(N_cheb + 1)
    A1_25 = im * R .* I(N_cheb + 1) * (gamma*Ma^2)^(-1) * rho .* I(N_cheb + 1)

    A1_31 = be * (lam + T) .* I(N_cheb + 1)
    A1_32 = im * R * rho .* F .* I(N_cheb + 1)
    A1_33 = zeros(N_cheb + 1, N_cheb + 1)
    A1_34 = zeros(N_cheb + 1, N_cheb + 1)
    A1_35 = zeros(N_cheb + 1, N_cheb + 1)

    A1_41 = -im * (rho .* (D*lam) .* I(N_cheb + 1) + (1 .+ lam .* rho) .* D)
    A1_42 = zeros(N_cheb + 1, N_cheb + 1)
    A1_43 = im * R * rho .* F .* I(N_cheb + 1)
    A1_44 = zeros(N_cheb + 1, N_cheb + 1)
    A1_45 = -im * rho .* ( (D*F)) .* I(N_cheb + 1)

    A1_51 = zeros(N_cheb + 1, N_cheb + 1)
    A1_52 = zeros(N_cheb + 1, N_cheb + 1)
    A1_53 = -2 * im * (gamma - 1) * Ma^2 * D * F .* I(N_cheb + 1)
    A1_54 = -(gamma - 1) * (gamma)^(-1) *  im * R * F  .* T .* I(N_cheb + 1)
    A1_55 = im * R * rho .* F .* I(N_cheb + 1) - (gamma - 1) * (gamma)^(-1) * im * R * F .* rho .* I(N_cheb + 1)

    A2_11 = zeros(N_cheb + 1, N_cheb + 1)
    A2_12 = zeros(N_cheb + 1, N_cheb + 1)
    A2_13 = zeros(N_cheb + 1, N_cheb + 1)
    A2_14 = zeros(N_cheb + 1, N_cheb + 1)
    A2_15 = zeros(N_cheb + 1, N_cheb + 1)

    A2_21 = (lam + 2 * T) .* I(N_cheb + 1)
    A2_22 = zeros(N_cheb + 1, N_cheb + 1)
    A2_23 = zeros(N_cheb + 1, N_cheb + 1)
    A2_24 = zeros(N_cheb + 1, N_cheb + 1)
    A2_25 = zeros(N_cheb + 1, N_cheb + 1)

    A2_31 = zeros(N_cheb + 1, N_cheb + 1)
    A2_32 = T .* I(N_cheb + 1)
    A2_33 = zeros(N_cheb + 1, N_cheb + 1)
    A2_34 = zeros(N_cheb + 1, N_cheb + 1)
    A2_35 = zeros(N_cheb + 1, N_cheb + 1)

    A2_41 = zeros(N_cheb + 1, N_cheb + 1)
    A2_42 = zeros(N_cheb + 1, N_cheb + 1)
    A2_43 = T .* I(N_cheb + 1)
    A2_44 = zeros(N_cheb + 1, N_cheb + 1)
    A2_45 = zeros(N_cheb + 1, N_cheb + 1)

    A2_51 = zeros(N_cheb + 1, N_cheb + 1)
    A2_52 = zeros(N_cheb + 1, N_cheb + 1)
    A2_53 = zeros(N_cheb + 1, N_cheb + 1)
    A2_54 = zeros(N_cheb + 1, N_cheb + 1)
    A2_55 = kappa .* I(N_cheb + 1)

    A0 = [A0_11 A0_12 A0_13 A0_14 A0_15 ; A0_21 A0_22 A0_23 A0_24 A0_25 ; A0_31 A0_32 A0_33 A0_34 A0_35 ; A0_41 A0_42 A0_43 A0_44 A0_45 ; A0_51 A0_52 A0_53 A0_54 A0_55  ]
    A1 = [A1_11 A1_12 A1_13 A1_14 A1_15 ; A1_21 A1_22 A1_23 A1_24 A1_25 ; A1_31 A1_32 A1_33 A1_34 A1_35 ; A1_41 A1_42 A1_43 A1_44 A1_45 ; A1_51 A1_52 A1_53 A1_54 A1_55  ]
    A2 = [A2_11 A2_12 A2_13 A2_14 A2_15 ; A2_21 A2_22 A2_23 A2_24 A2_25 ; A2_31 A2_32 A2_33 A2_34 A2_35 ; A2_41 A2_42 A2_43 A2_44 A2_45 ; A2_51 A2_52 A2_53 A2_54 A2_55  ]

    A0 = A0[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5))]
    A1 = A1[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5))]
    A2 = A2[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5))]

    return A0,A1,A2
 end
function DevelopingSpatialMode()
    c = -2/3
    eye = I(N_cheb + 1)
    Zero = zeros(N_cheb + 1, N_cheb + 1)
    A0_11 = (1/R) * rho .* eye
    A1_11 = im * rho .* eye
    A2_11 = Zero

    A0_12 = im * be * rho .* eye
    A1_12 = Zero
    A2_12 = Zero

    A0_13 = rho .* D + D * rho .* eye
    A1_13 = Zero
    A2_13 = Zero

    A0_14 = im * be * G .* eye .- im * omega + (1/R) * H .* D + (2*F)/R .* eye + (D*H)/R .* eye
    A1_14 = im * F .* eye
    A2_14 = Zero

    A0_15 = Zero
    A1_15 = Zero
    A2_15 = Zero

    A0_21 = im * be * rho .* G .* eye - im * omega .* rho .* eye + (rho.*F)./R .* eye + ((rho.*H)./R) .* D + be^2 / R .* T .* eye - D * T .* D - T .* D2
    A1_21 = im * rho .* F .* eye
    A2_21 = ((2+c)/R) * T .* eye

    A0_22 = -(2/R) * rho .* G .* eye - 2 * rho .* eye
    A1_22 = (c/R) * be * T .* eye + (be/R) * T .* eye
    A2_22 = Zero

    A0_23 = rho .* D*F .* eye 
    A1_23 = -(1/R) * (im * c * T .* D + im * D * T .* eye + im * T .* D)
    A2_23 = Zero

    A0_24 = F.^2 .* eye + (1/R) * H .* D * F .* eye - (1/R) * G.^2 .* eye - 2 * G .* eye - R .* eye
    A1_24 = im * (gamma * Ma^2)^(-1) * T .* eye
    A2_24 = Zero

    A0_25 = -(D2 * F .* eye + D*F .* D)
    A1_25 = im * (gamma * Ma^2)^(-1) * rho .* eye
    A2_25 = Zero

    A0_31 = (1/R) * (2 * rho .* G .* eye + 2 * rho .* eye - c * T * im * be /R .* eye)  
    A1_31 = (c/R) * be * T .* eye
    A2_31 = Zero

    A0_32 = im * be * rho .* G .* eye - im * omega .* rho .* eye + (1/R) * ( rho .* H .* D + rho .* F .* eye + c * be^2 *T .* eye + im * be * (c/R) * T .* eye + 2 * be^2 * T .* eye) - D*T .* D - T .* D2
    A1_32 = im * rho .* F .* eye
    A2_32 = Zero

    A0_33 = rho .* (D*G) .* eye - (c/R^2) * im * be * T .* D - (im*be/R) * D * T .* eye - im * be / R * T .* D
    A1_33 = Zero
    A2_33 = Zero

    A0_34 = (1/R) * (2 * F .* G .* eye + H .* D * G .* eye + 2 * F .* eye ) + im * be * (gamma * Ma^2)^(-1) * T .* eye
    A1_34 = Zero
    A2_34 = Zero

    A0_35 = -(D2 * G .* eye + D*G .* D) + im * be *(gamma * Ma^2)^(-1) * rho .* eye
    A1_35 = Zero
    A2_35 = Zero

    A0_41 = - (c/R^2) * T .* D
    A1_41 = -(im * c/R) * T .* D - im * T .* D
    A2_41 = Zero

    A0_42 = (-1/R) * (im * be * c * T .* D + im * be * T .* D)
    A1_42 = Zero
    A2_42 = Zero

    A0_43 = im * be * rho .* G .* eye - im * omega .* rho .* eye + (1/R) * (rho .* H .* D + rho .* D * H .* eye - c * T .* D2 + be^2 * T .* eye) - 2 * D * T .* D - 2 * T .* D2
    A1_43 = im * rho .* F .* eye
    A2_43 = (1/R) * T .* eye

    A0_44 = (H .* D * H .* eye) / (R^2) + (gamma * Ma^2)^(-1) * (D*T .* eye + T .* D)
    A1_44 = Zero
    A2_44 = Zero

    A0_45 = - (1/R) * ((c/R) * D * F .* eye + (c/R) * D2 * H .* eye + im * be * D * G .* eye + 2 * D2 * H .* eye + 2 * D * H .* D) + (gamma * Ma^2)^(-1) * (D * rho .* eye + rho .* D)
    A1_45 = - im * D * F .* eye
    A2_45 = Zero

    A0_51 = -(1/R^2) * (2 * Ma^2 * (gamma - 1)) * im * be * T .* G .* eye - (2 * Ma^2 * (gamma - 1)/R) * R * (T .* D * F .* D) - (c * Ma^2 * (gamma - 1) / R^2) * (4 * T .* F .* eye + D * H .* T .* eye)
    A1_51 = -(4 * im * Ma^2 * (gamma - 1) / R^2) * T .* F .* eye - (c * Ma^2 * (gamma - 1) / R^2) * (4 * im * R  * T .* F .* eye + im * R * D * H .* F.* eye )
    A2_51 = Zero

    A0_52 = - (Ma^2 * (gamma - 1) / R) * (4 * im * be / R * T .* F .* eye + 2 * R * D * G .* T .* D + (c/R) * (4 * im * be * R * T .* F .* eye + 2 * im * be * R * D * H .* eye))
    A1_52 = -(2 * im * Ma^2 * (gamma - 1) / R^2) * T .* G .* eye
    A2_52 = Zero

    A0_53 = D * T .* rho .* eye - (Ma^2*(gamma - 1) / R) * ((4/R) * T .* D * H .* eye + 2 * im * be * T .* D * G .* eye + (c/R) * (4 * R * T .* F .* D + 2 * R * D * H .* D))
    A1_53 = - (Ma^2*(gamma - 1) / R) * (2 * im *  T .* D * F) .* eye
    A2_53 = Zero

    A0_54 = (1/R) * H .* D * T .* eye + (gamma-1)/(gamma) * ((im * be * G .* eye - im * omega .* eye) .* T + H .* D * T .* eye + H .* T .* D) 
    A1_54 = (gamma-1)/(gamma) * (im * F .* T .* eye)
    A2_54 = Zero

    A0_55 = im * be * rho .* G .* eye -im * omega * rho .* eye + (1/R) * (rho .* H .* D ) + (1/sigma) * (be^2 .* eye - D^2 ) - (Ma^2 * (gamma - 1) / R) * (4 * F.^2 / R^2 .* eye + (2 * D * H + G.^2) / R^2 .* eye 
    + R * (D*F).^2 .* eye + R * (D*G).^2 .* eye + (c/R) * (4 * F.^2 .* eye + 4 * F .* D * H .* eye + (D*H).^2 .* eye)) + (gamma-1)/(gamma) * (im * be * G .* eye -im * omega .* eye) .* rho + (gamma-1)/(gamma) * (H .* D * rho .* eye + H .* rho .* D )
    A1_55 = im * rho .* F .* eye - (im/(sigma * R^2) .* eye) + (gamma-1)/gamma * im * F .* rho .* eye
    A2_55 = 1/(sigma * R) .* eye

    A0 = [A0_11 A0_12 A0_13 A0_14 A0_15 ; A0_21 A0_22 A0_23 A0_24 A0_25 ; A0_31 A0_32 A0_33 A0_34 A0_35 ; A0_41 A0_42 A0_43 A0_44 A0_45 ; A0_51 A0_52 A0_53 A0_54 A0_55  ]
    A1 = [A1_11 A1_12 A1_13 A1_14 A1_15 ; A1_21 A1_22 A1_23 A1_24 A1_25 ; A1_31 A1_32 A1_33 A1_34 A1_35 ; A1_41 A1_42 A1_43 A1_44 A1_45 ; A1_51 A1_52 A1_53 A1_54 A1_55  ]
    A2 = [A2_11 A2_12 A2_13 A2_14 A2_15 ; A2_21 A2_22 A2_23 A2_24 A2_25 ; A2_31 A2_32 A2_33 A2_34 A2_35 ; A2_41 A2_42 A2_43 A2_44 A2_45 ; A2_51 A2_52 A2_53 A2_54 A2_55  ]

    A0 = A0[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5))]
    A1 = A1[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5))]
    A2 = A2[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5))]

 end