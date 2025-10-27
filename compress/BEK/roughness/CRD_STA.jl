" 
    ∊ Designed by hj_zhang on Nov 16,2024,TJU 
    ∊ Add Stability Analysis on Jan 16,2025,TJU
                                                "
module CRD_BF

    export sol_baseflowODE,velocity,Cheb,phi_var,f_q,T_var,Physical_Interpretation

    using LinearAlgebra
    using DifferentialEquations
    using BSplineKit
    using IterativeSolvers
    using PyCall
    
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

 function Cheb(N)
        θ = range(0,length=N+1,stop=pi)
        x = reshape(-cos.(θ), N+1, 1)
        c = [2; ones(N-1, 1) ; 2] .* (-1) .^ (0:N)
        X = repeat(x, 1, N+1);
        dX = X - X';
        D = (c * (1 ./ c)') ./ (dX .+ I(N+1));
        D = D - diagm(vec(sum(D, dims=2))); 
        a = 1
        b = 0.6
        c = 0.5
        for i=1:N+1
            D[i,:]=D[i,:].* (1-b*x[i]-(1-b)*(x[i]^3+c*(1-x[i]^2)))^2/(2a*(b .+ 3 * (1-b)*x[i]^2 - 2 * c * (1-b) * x[i]))
        end
        for i=1:N+1
            x[i] = a * (1+b*x[i]+(1-b)*(x[i]^3+c*(1-x[i]^2)))/(1-b*x[i]-(1-b)*(x[i]^3+c*(1-x[i]^2)))
            if x[i]>30
                x[i]=30
            end
        end
        
        D2 = D^2;

        return D,D2,x

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
using DifferentialEquations
using BSplineKit
using LinearAlgebra
struct COF
        Ta :: Matrix{ComplexF64}
        A :: Matrix{ComplexF64}
        B :: Matrix{ComplexF64}
        C :: Matrix{ComplexF64}
        D1 :: Matrix{ComplexF64}
        Vxx :: Matrix{ComplexF64}
        Vyy :: Matrix{ComplexF64}
        Vzz :: Matrix{ComplexF64}
        Vxy :: Matrix{ComplexF64}
        Vxz :: Matrix{ComplexF64}
        Vyz :: Matrix{ComplexF64}
end
function baseflow_var(N_cheb,Ro,Co)
    N = 20000
    tspan = (0,30)
    t = range(0,30,N)
    sigma = 0.72
    u0,v0,w0,du0,dv0,x = CRD_BF.sol_baseflowODE(Ro)
    PHI = CRD_BF.phi_var(u0,t.step.hi,N)
    u0,du0,v0,dv0,w0,F_u,F_du,F_dv,F_w,F_phi = CRD_BF.velocity(u0,v0,w0,du0,dv0,t,PHI)
    D,D2,x = CRD_BF.Cheb(N_cheb)
    f,q = CRD_BF.f_q(sigma,F_du,F_dv,F_u,F_phi,tspan,t)
    if Ro >= 0
        u0 = -1 * u0
        v0 = -1 * v0
        w0 = -1 * w0
    end
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
function Spatial_mode_BEK(F,G,H,rho,lam,kappa,T,sigma,gamma,R,Ma,N_cheb,Ro,Co,D,D2)
    eye = I(N_cheb + 1)
    Zero = zeros(N_cheb + 1,N_cheb + 1)
    if Ro == -1 
        Ro = 1
    end
    Ta_11 = Ta_12 = Ta_13 = Ta_15 = Ta_22 = Ta_23 = Ta_24 = Ta_25 = Ta_31 = Ta_33 = Ta_34 = Ta_35 = Ta_41 = Ta_42 = Ta_44 = Ta_45 = Ta_51 = Ta_52 = Ta_53 = Zero
    Ta_14 = R * eye
    Ta_21 = R .* rho .* eye
    Ta_32 = R .* rho .* eye
    Ta_43 = R .* rho .* eye
    Ta_54 = (1-gamma)/(gamma) * R .*T .* eye
    Ta_55 = 1/(gamma) * R .*rho .* eye

    A_12 = A_13 = A_15 = A_22 = A_31 = A_33 = A_34 = A_35 = A_42 = A_44 = A_51 = A_52 = Zero
    A_11 = R * rho .* eye
    A_14 = R * F .* eye
    A_21 = R * rho .* F .* eye
    A_23 = - rho .* (D*T) .* eye
    A_24 = (gamma*Ma^2)^(-1) * R .* T .* eye
    A_25 = (gamma*Ma^2)^(-1) * R .* rho .* eye
    A_32 = R * rho .* F .* eye
    A_41 = -rho .* (D*lam) .* eye
    A_43 = R * rho .* F .* eye
    A_45 = - rho .* (D*F) .* eye
    A_53 = -2 * (gamma-1) * Ma^2 * (D*F) .* eye 
    A_54 = -(gamma-1)/(gamma) * R .* T .* F .* eye
    A_55 = 1/gamma * R * rho .* F .* eye

    B_11 = B_13 = B_15 = B_22 = B_23 = B_24 = B_25 = B_31 = B_41 = B_44 = B_51 = B_52 = Zero
    B_12 = R* rho .* eye
    B_14 = R * G .* eye
    B_21 = R * rho .* G .* eye
    B_32 = R * rho .* G .* eye
    B_33 = -rho .* (D*T) .* eye
    B_34 = (gamma * Ma^2)^(-1)  * R .* T .* eye
    B_35 = (gamma * Ma^2)^(-1)  * R .* rho .* eye
    B_42 = -rho .* (D*lam) .* eye
    B_43 = R * rho .* G .* eye
    B_45 = -rho .* (D*G) .* eye
    B_53 = -2 * (gamma-1) * Ma^2 * (D*G) .* eye
    B_54 = -(gamma-1)/(gamma) * R * T .* G .* eye
    B_55 = 1/gamma * R .* rho .* G .* eye

    C_11 = C_12 = C_15 = C_22 = C_23 = C_24 = C_31 = C_33 = C_34 = C_41 = C_42 = C_53 = Zero
    C_13 = R * rho.^2 .* eye
    C_14 = rho .* H .* eye
    C_21 = Ro * rho.^2  .* H .* eye
    C_25 = -rho.^2 .* (D*F) .* eye
    C_32 = Ro * rho.^2 .* H .* eye
    C_35 = -rho.^2 .* (D*G) .* eye
    C_43 = Ro * rho.^2 .* H .* eye
    C_44 = R * (gamma * Ma^2)^(-1) * rho .* T .* eye
    C_45 = R * (gamma * Ma^2)^(-1) * rho .* rho .* eye
    C_51 = -2 * (gamma-1) * Ma^2 * rho .* (D*F) .* eye
    C_52 = -2 * (gamma-1) * Ma^2 * rho .* (D*G) .* eye
    C_54 = -(gamma-1)/gamma * rho .* H .* T .* eye
    C_55 = 1/gamma * rho.^2 .* H .* eye + 1/sigma * rho.^2 .* (D*T) .* eye

    D_12 = D_15 = D_41 = D_42 = D_51 = D_52 = Zero
    D_11 = rho .* eye
    D_13 = R * rho .* (D*rho) .* eye
    D_14 = 2 * F .* eye + rho .* (D*H) .* eye
    D_21 = Ro * rho .* F .* eye
    D_22 = -rho .* (2* Ro * G .+ Co) .* eye
    D_23 = R * rho.^2 .* (D*F) .* eye
    D_24 = D*rho .* (D2*F) .* eye
    D_25 = -rho .* D*(rho.*(D*F)) .* eye
    D_31 = rho .* (2* Ro * G .+ Co) .* eye
    D_32 = Ro * rho .* F .* eye
    D_33 = R * rho.^2 .* (D*G) .* eye
    D_34 = F .* (2* Ro * G .+ Co) .* eye + Ro * rho .* H .* (D*G) .* eye
    D_35 = -rho.* (D * (rho .* (D*G))) .* eye
    D_43 = Ro * rho.^2 .* (D*H) .* eye
    D_44 = R * (gamma * Ma^2)^(-1) * rho .* (D*T) .* eye
    D_45 = R * (gamma * Ma^2)^(-1) * rho .* (D*rho) .* eye
    D_53 = R * rho.^2 .* (D*T) .* eye
    D_54 = 1/gamma * rho .* H .* (D*T) .* eye
    D_55 = -(gamma-1)/gamma * rho .* H .* (D*rho) .* eye - (1/sigma) * (rho .* (D*rho) .* (D*T) .+ rho.^2 .* (D2 * T)) .* eye - (gamma-1) * Ma^2 * rho.^2 .* ((D*F).^2 + (D*G).^2) .* eye

    Vxx_11 = Vxx_12 = Vxx_13 = Vxx_14 = Vxx_15 = Vxx_22 = Vxx_23 = Vxx_24 = Vxx_25 = Vxx_31 = Vxx_33 = Vxx_34 = Vxx_35 = Vxx_41 = Vxx_42 = Vxx_44 = Vxx_45 = Vxx_51 = Vxx_52 = Vxx_53 = Vxx_54 = Zero
    Vxx_21 = -(lam + 2*T) .* eye
    Vxx_32 = -T .* eye
    Vxx_43 = -T .* eye
    Vxx_55 = -kappa .* eye

    Vyy_11 = Vyy_12 = Vyy_13 = Vyy_14 = Vyy_15 = Vyy_22 = Vyy_23 = Vyy_24 = Vyy_25 = Vyy_31 = Vyy_33 = Vyy_34 = Vyy_35 = Vyy_41 = Vyy_42 = Vyy_44 = Vyy_45 = Vyy_51 = Vyy_52 = Vyy_53 = Vyy_54 = Zero
    Vyy_21 = -T .* eye
    Vyy_32 = -(lam .+ 2*T) .* eye
    Vyy_43 = -T .* eye
    Vyy_55 = -kappa .* eye

    Vzz_11 = Vzz_12 = Vzz_13 = Vzz_14 = Vzz_15 = Vzz_22 = Vzz_23 = Vzz_24 = Vzz_25 = Vzz_31 = Vzz_33 = Vzz_34 = Vzz_35 = Vzz_41 = Vzz_42 = Vzz_44 = Vzz_45 = Vzz_51 = Vzz_52 = Vzz_53 = Vzz_54 = Zero
    Vzz_21 = -rho .* eye
    Vzz_32 = -rho .* eye
    Vzz_43 = -rho .* (2 .+ lam .* rho) .* eye
    Vzz_55 = -rho.^2 .* kappa .* eye

    Vxy_11 = Vxy_12 = Vxy_13 = Vxy_14 = Vxy_15 = Vxy_21 = Vxy_23 = Vxy_24 = Vxy_25 = Vxy_32 = Vxy_33 = Vxy_34 = Vxy_35 = Vxy_41 = Vxy_42  = Vxy_43 = Vxy_44 = Vxy_45 = Vxy_51 = Vxy_52 = Vxy_53 = Vxy_54 = Vxy_55 =  Zero
    Vxy_22 = -(lam .+ T) .* eye
    Vxy_31 = -(lam .+ T) .* eye

    Vxz_11 = Vxz_12 = Vxz_13 = Vxz_14 = Vxz_15 = Vxz_21 = Vxz_22 = Vxz_24 = Vxz_25 = Vxz_31 = Vxz_32 = Vxz_33 = Vxz_34 = Vxz_35 = Vxz_42 = Vxz_43 = Vxz_44 = Vxz_45 = Vxz_51 = Vxz_52 = Vxz_53 = Vxz_54 = Vxz_55 =  Zero
    Vxz_23 = - (1 .+ rho.*lam) .* eye
    Vxz_41 = - (1 .+ rho.*lam) .* eye

    Vyz_11 = Vyz_12 = Vyz_13 = Vyz_14 = Vyz_15 = Vyz_21 = Vyz_22 = Vyz_23 = Vyz_24 = Vyz_25 = Vyz_31 = Vyz_32  = Vyz_34 = Vyz_35 = Vyz_41 = Vyz_43 = Vyz_44 = Vyz_45 = Vyz_51 = Vyz_52 = Vyz_53 = Vyz_54 = Vyz_55 =  Zero
    Vyz_33 = - (1 .+ rho.*lam) .* eye
    Vyz_42 = - (1 .+ rho.*lam) .* eye

    Ta = [Ta_11 Ta_12 Ta_13 Ta_14 Ta_15;Ta_21 Ta_22 Ta_23 Ta_24 Ta_25;Ta_31 Ta_32 Ta_33 Ta_34 Ta_35;Ta_41 Ta_42 Ta_43 Ta_44 Ta_45;Ta_51 Ta_52 Ta_53 Ta_54 Ta_55]

    A = [A_11 A_12 A_13 A_14 A_15;A_21 A_22 A_23 A_24 A_25;A_31 A_32 A_33 A_34 A_35;A_41 A_42 A_43 A_44 A_45;A_51 A_52 A_53 A_54 A_55]

    B = [B_11 B_12 B_13 B_14 B_15;B_21 B_22 B_23 B_24 B_25;B_31 B_32 B_33 B_34 B_35;B_41 B_42 B_43 B_44 B_45;B_51 B_52 B_53 B_54 B_55]

    C = [C_11 C_12 C_13 C_14 C_15;C_21 C_22 C_23 C_24 C_25;C_31 C_32 C_33 C_34 C_35;C_41 C_42 C_43 C_44 C_45;C_51 C_52 C_53 C_54 C_55]

    D1 = [D_11 D_12 D_13 D_14 D_15;D_21 D_22 D_23 D_24 D_25;D_31 D_32 D_33 D_34 D_35;D_41 D_42 D_43 D_44 D_45;D_51 D_52 D_53 D_54 D_55]

    Vxx = [Vxx_11 Vxx_12 Vxx_13 Vxx_14 Vxx_15;Vxx_21 Vxx_22 Vxx_23 Vxx_24 Vxx_25;Vxx_31 Vxx_32 Vxx_33 Vxx_34 Vxx_35;Vxx_41 Vxx_42 Vxx_43 Vxx_44 Vxx_45;Vxx_51 Vxx_52 Vxx_53 Vxx_54 Vxx_55]

    Vyy = [Vyy_11 Vyy_12 Vyy_13 Vyy_14 Vyy_15;Vyy_21 Vyy_22 Vyy_23 Vyy_24 Vyy_25;Vyy_31 Vyy_32 Vyy_33 Vyy_34 Vyy_35;Vyy_41 Vyy_42 Vyy_43 Vyy_44 Vyy_45;Vyy_51 Vyy_52 Vyy_53 Vyy_54 Vyy_55]

    Vzz = [Vzz_11 Vzz_12 Vzz_13 Vzz_14 Vzz_15;Vzz_21 Vzz_22 Vzz_23 Vzz_24 Vzz_25;Vzz_31 Vzz_32 Vzz_33 Vzz_34 Vzz_35;Vzz_41 Vzz_42 Vzz_43 Vzz_44 Vzz_45;Vzz_51 Vzz_52 Vzz_53 Vzz_54 Vzz_55]

    Vxy = [Vxy_11 Vxy_12 Vxy_13 Vxy_14 Vxy_15;Vxy_21 Vxy_22 Vxy_23 Vxy_24 Vxy_25;Vxy_31 Vxy_32 Vxy_33 Vxy_34 Vxy_35;Vxy_41 Vxy_42 Vxy_43 Vxy_44 Vxy_45;Vxy_51 Vxy_52 Vxy_53 Vxy_54 Vxy_55]

    Vxz = [Vxz_11 Vxz_12 Vxz_13 Vxz_14 Vxz_15;Vxz_21 Vxz_22 Vxz_23 Vxz_24 Vxz_25;Vxz_31 Vxz_32 Vxz_33 Vxz_34 Vxz_35;Vxz_41 Vxz_42 Vxz_43 Vxz_44 Vxz_45;Vxz_51 Vxz_52 Vxz_53 Vxz_54 Vxz_55]

    Vyz = [Vyz_11 Vyz_12 Vyz_13 Vyz_14 Vyz_15;Vyz_21 Vyz_22 Vyz_23 Vyz_24 Vyz_25;Vyz_31 Vyz_32 Vyz_33 Vyz_34 Vyz_35;Vyz_41 Vyz_42 Vyz_43 Vyz_44 Vyz_45;Vyz_51 Vyz_52 Vyz_53 Vyz_54 Vyz_55]
    return COF(Ta,A,B,C,D1,Vxx,Vyy,Vzz,Vxy,Vxz,Vyz)
end
function assemble_mat(cof :: COF,D,D2,be,omega)
    L0 = cof.D1  + im * be * cof.B - im * omega * cof.Ta - be^2 * cof.Vyy + (cof.C .+ im * be * cof.Vyz) * kron(I(5), D)  + (cof.Vzz) * kron(I(5),D2) 
    L1 = im * cof.A - be * cof.Vxy + im *  cof.Vxz * kron(I(5),D)
    L2 = -cof.Vxx 
    L0 = L0[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5))]
    L1 = L1[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5))]
    L2 = L2[setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5)),setdiff(1:end , (1,N_cheb + 1,N_cheb + 2,2N_cheb + 2,2N_cheb + 3,3N_cheb + 3,3N_cheb + 4,4N_cheb + 4,4N_cheb + 5,5N_cheb + 5))]
    return L0,L1,L2
end
