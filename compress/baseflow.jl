
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
    function velocity_CP(t,u)

        U = u[1 , :]
        dU = u[2 , :]
        V = u[3 , :]
        dV = u[4 , :]
        W = u[5 , :]
        F_U = itp = interpolate(t, U , BSplineOrder(4))
        F_dU = itp = interpolate(t, dU , BSplineOrder(4))
        F_dV = itp = interpolate(t, dV , BSplineOrder(4))
        F_W= itp = interpolate(t, W , BSplineOrder(4))
        
        return U,dU,V,dV,W,F_U,F_dU,F_dV,F_W
     end
    function velocity_SL(t,u)
        U = u[1 , :]
        V = u[2 , :]
        W = u[3 , :]
        T = u[4 , :]
        dU = u[5 , :]
        dV = u[6 , :]
        dT = u[7 , :]
        return U,V,W,T,dU,dV,dT
     end
    function T_start(F_dU,F_dV,F_W,sigma,gamma,Tw,tspan,Num)
        gamma = gamma
        Tw = Tw 
        sigma = sigma
        function ODE_T!(du,u,p,t)
            T = u[1]
            dT = u[2]
            du[1] = dT
            du[2] = (F_W(t) * u[2] )*sigma
        end
        function bc2!(residual, u, p, t)
            residual[1] = u[begin][1] - Tw
            residual[2] = u[end][1] - 1
        end
            prob = BVProblem(ODE_T! , bc2! , [0,0] , tspan)
            sol = solve(prob , Shooting(Vern7()) , dt = 0.001)
            t=range(0.0, 20, Num)
            T=sol(t)
            dT = T[2,:]
            T = T[1,:]
        return T , dT

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
    function y_DiffMat(y_range,N)
    y0 = minimum(y_range)
    y1 = maximum(y_range)
    dely = (y1 - y0) / (N-1)
    D = zeros(N,N)
    D2 = zeros(N,N)
    for i = 1:N
        if i == 1
            D[i,i+2] = -1 / (2*dely)
            D[i,i+1] = 4 / (2*dely)
            D[i,i] = -3 / (2*dely)
            D2[i,i] = 2 / (dely^2)
            D2[i,i+1] = -5 / (dely^2)
            D2[i,i+2] = 4 / (dely^2)
            D2[i,i+3] = -1 / (dely^2)
        elseif i == N 
            D[i,i] = 3 / (2*dely)
            D[i,i-1] = -4 / (2*dely)
            D[i,i-2] = 1 / (2*dely)
            D2[i,i] = 2 / (dely^2)
            D2[i,i-1] = -5 / (dely^2)
            D2[i,i-2] = 4 / (dely^2)
            D2[i,i-3] = -1 /(dely^2)
        elseif i != 1;2;N-1;N
            D[i,i-1] = -1 / (2*dely)
            D[i,i+1] = 1 / (2*dely)
            D2[i,i] = -2 / (dely^2)
            D2[i,i+1] = 1 / (dely^2)
            D2[i,i-1] = 1 / (dely^2)
        end
    end

    return D,D2

    end

    function Physical_Interpretation(T,delt,Num)
        z = zeros(Num,1)
        integral_T = delt * T
        for i = 1 : 1 : Num
            z[i,1] = sum(integral_T[1:i])
        end
        return z
 end

using.CRD_BF
using DelimitedFiles

global N = 1001
global x_step = 10^-2
global y_step = 2 * 10^-2
global tspan = (0,20)
global sigma = 0.7
global Tw = 0.5
global gamma = 1.4
global r = x_step/(y_step^2)
global q = x_step / (2*y_step)
global data = empty!

println("type in viscouslaw,x_max,Ma,Tw")
parameters = readline(stdin)
c = split(parameters,',')
viscouslaw = c[1,1]
Ma = parse(Int64,c[3,1])
x_max = parse(Int64,c[2,1])
it_num = Int(round(x_max/x_step))
A = zeros(N,N)
t = range(0,20,N)
m_left = zeros(N,N)
m_right = zeros(N,N)
D,D2 = y_DiffMat(t,N)
Title = " Variables= \"x\" \"z\" \"eta\" \"u0\" \"v0\" \"w0\" \"W\"  \"Temp\" \"mu\" 
        \n Zone T=1 I=1001 J=$(it_num) "
if viscouslaw == "CP"
    local u,z = CRD_BF.sol_baseflowODE(tspan,N)
    local u0,du0,v0,dv0,w0,F_u,F_du,F_dv,F_w = CRD_BF.velocity_CP(t,u)
    local T0,dT0 = CRD_BF.T_start(F_du,F_dv,F_w,sigma,gamma,Tw,tspan,N)
    local data = empty!
    for j = 1 : 1 : it_num
            for i = 1 : 1 : N
            if i == 1
                    m_left[i,i] = 1
            elseif i == N
                    m_left[N,N] = 1
            else
                    m_left[i,i-1] = -(r/sigma + q * w0[i,1])
                    m_left[i,i] = j*x_step*u0[i,1] + 2 * (r/sigma) 
                    m_left[i,i+1] = (q * w0[i,1]-r/sigma)
            end
            end
            T_temp = T0
            local A = Ma^2 * (gamma-1) .*(du0 .* du0 + dv0 .* dv0)
            local m_right = x_step * (j*x_step)^2 * A + j*x_step * u0 .* T_temp
            m_right[1,1] = Tw
            m_right[end,1] = 1
            T0 = m_left\m_right
            x = x_step * j * ones(N,1)
            eta = Physical_Interpretation(T0,y_step,N)
            W = w0 ./ T0
            data_temp = [x t eta u0 v0 w0 W T0 H]
            data = [data;data_temp]  
    end
    data = data[2:end,:]
    data_full = [Title empty empty empty empty empty empty empty empty; data]
    writedlm("Ma = $(Ma)_Tw = $(Tw)_1.dat",data_full,'\t')

elseif viscouslaw == "SL"
    local u,z = CRD_BF.In_Su(tspan,N)
    local u0,v0,w0,T0,du0,dv0,dT0 = CRD_BF.velocity_SL(t,u)
    local data = empty!
    for j = 1 : 1 : it_num
            G = (T0.^(1/2).*(273+114))./(T0.*273 .+ 114)       
            dG = D * G
            for i = 1 : 1 : N
            if i == 1
                    m_left[i,i] = 1
            elseif i == N
                    m_left[N,N] = 1
            else
                    m_left[i,i-1] = -((r/sigma).*G[i,1] - (q/sigma)*dG[i,1] + q * w0[i,1])
                    m_left[i,i] = j*x_step*u0[i,1] + 2 * (r/sigma) *G[i,1]
                    m_left[i,i+1] = (q * w0[i,1]-(r/sigma).*G[i,1]-(q/sigma).*dG[i,1])
            end
            end
            T_temp = T0
            local A = Ma^2 * (gamma-1) .*(du0 .* du0 + dv0 .* dv0)
            local m_right = x_step * (j*x_step)^2 * A .*G + j*x_step * u0 .* T_temp
            m_right[1,1] = Tw
            m_right[end,1] = 1
            T0 = m_left\m_right
            x = x_step * j * ones(N,1)
            H = (T0.^(3/2).*(273+114))./(T0.*273 .+ 114)
            eta = Physical_Interpretation(T0,y_step,N)
            W = w0 ./ T0
            data_temp = [x t eta u0 v0 w0 W T0 H]
            data = [data;data_temp]  
    end
    data = data[2:end,:]
    data_full = [Title empty empty empty empty empty empty empty empty; data]
    DelimitedFiles.writedlm("Ma = $(Ma)_Tw = $(Tw)_1.dat",data_full,'\t')
else
            println("Undifine law,this module only support Sutherland Law(SL) and Chapman Law(CP)")
 end
 end