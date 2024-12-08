" Designed by hj_zhang on Dec 13,2024,TJU "

module CRD_BF

    export In_Su,sol_baseflowODE,velocity,T_start,Cheb,y_DiffMat

    using LinearAlgebra
    using BoundaryValueDiffEq
    using DifferentialEquations
    using BSplineKit
    using IterativeSolvers

function In_Su(tspan,Num)
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
    F_phi = interpolate(t, PHI , BSplineOrder(4))
    return U,dU,V,dV,W,F_U,F_dU,F_dV,F_W,F_phi

 end

function Cheb(u,v,w,T,t,N)
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
        for i=1:N+1
            D[i,:]=D[i,:].*((2*x[i]^3-x[i]^2+3*x[i]-4)^2/(20*(6*x[i]^2-2*x[i]+3)))
        end
        for i=1:N+1
            x[i]=(4*x[i]^3-2*x[i]^2+6*x[i]+12)/(-2*x[i]^3+x[i]^2-3*x[i]+4)
            if x[i]>20
                x[i]=20
            end
        end
        itpw = itp = interpolate(t, w , BSplineOrder(4))
        itpu = itp = interpolate(t, u , BSplineOrder(4))
        itpv = itp = interpolate(t, v , BSplineOrder(4))
        itpT = itp = interpolate(t, T , BSplineOrder(4))
        u = zeros(N+1,1)
        v = zeros(N+1,1)
        w = zeros(N+1,1)
        T = zeros(N+1,1)
        for i=1:N+1
            u[i,1] = itpu(x[i])
            v[i,1] = itpv(x[i])
            w[i,1] = itpw(x[i])
            T[i,1] = itpT(x[i])
        end
        return D,x,u,v,w,T
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
function f_q(sigma,gamma,Tw,F_du,F_dv,F_u,F_phi)
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
    prob = BVProblem(ODE_f!, BC_f!, [0.0, 0.0], tspan)
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
    prob = BVProblem(ODE_q!, BC_q!, [1,0], tspan)
    sol1 = solve(prob,MIRK4(), dt=0.01)
    q = sol1(t)
    f = f[1,:]
    q = q[1,:]
    return f,q
 end
function T_var(Mx,f,q,Tw)
    T = 1 .- ( (gamma-1)/2 )*Mx^2 * f + (Tw - 1) * q
    return T
 end
using.CRD_BF
function baseflow_var(R,Ma,Tw,N_cheb)
    N = 1001
    tspan = (0,20)
    t = range(0,20,N)
    sigma = 0.7
    gamma = 1.4
    Mx = Ma * R
    u,z = CRD_BF.sol_baseflowODE(tspan,N)
    PHI = phi_var(u,t.step.hi,N)
    u0,du0,v0,dv0,w0,F_u,F_du,F_dv,F_w = CRD_BF.velocity(u,t,PHI)
    f,q = f_q(sigma,gamma,Tw,F_du,F_dv,F_u,F_phi)
    T = T_var(Mx,f,q,Tw)    
    D,x,F,G,W,T = Cheb(u0,v0,w0,T,t,N_cheb)
    RHO = 1 ./ T
    H = W ./ T
    return RHO,F,G,H,T,D,x
 end