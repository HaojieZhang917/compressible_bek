module CRC_BF
    using LinearAlgebra
    using DifferentialEquations
    using BSplineKit
    using IterativeSolvers
    using PyCall
    function BaseFlow(Re_s)
        py"""
        import numpy as np
        from scipy.integrate import solve_bvp
        import matplotlib.pyplot as plt
        import math
        from math import sqrt
        Re_s = $Re_s
        def oneDiskODE(z,y):
                # Y0 = H, Y1 = F,Y2 = F', Y3 = F'', Y4 = G, Y5 = G'
                dH = -2 * sqrt(Re_s) * y[1]
                dydz = np.zeros((6, len(z)))
                dydz = np.array([dH , y[2] , y[3] , Re_s * ( (1/sqrt(Re_s)) * ( y[3] * y[0] + y[2] * dH) - 2 * y[4] * y[5] + 2 * y[1] * y[2]) , y[5] , Re_s * ( (1/sqrt(Re_s)) * y[5] * y[0] + 2 * y[1] * y[4])])
                return dydz 

        def oneDiskBC(ya, yb):
                resa = np.array([ya[0],ya[1], ya[4] - 0])
                
                resb = np.array([yb[0],yb[1], yb[4] - 1])
                
                return np.concatenate((resa, resb))


        z = np.linspace(0, 1, 20000)
        y = np.zeros((6, len(z)))
        y_guess = np.zeros((6, z.size))
        y_guess[0] = 1
        y_guess[1] = 0
        y_guess[2] = 0
        y_guess[3] = 0
        y_guess[4] = 1
        y_guess[5] = 0
        solution = solve_bvp(oneDiskODE, oneDiskBC, z, y_guess,max_nodes=5000000)

        x_plot = np.linspace(0, 1, 40000)


        y1_plot = solution.sol(x_plot)[0]
        y2_plot = solution.sol(x_plot)[1]
        y3_plot = solution.sol(x_plot)[4]
        y4_plot = solution.sol(x_plot)[2]
        y5_plot = solution.sol(x_plot)[5]
        """
        w0 = py"y1_plot"
        u0 = py"y2_plot"
        v0 = py"y3_plot"
        du0 = py"y4_plot"
        dv0 = py"y5_plot"
        x = py"x_plot"
        return u0,v0,w0,du0,dv0,x
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
        D = 2 * D
        x = 0.5 * (x .+1) 
        D2 = D^2;
        return D,D2,x
    end
end