module CRC_BF
    using LinearAlgebra
    using BSplineKit
    using IterativeSolvers
    using PyCall
 function BaseFlow(Re_s,Ro,Ts,mode)
    if mode == 1 
        py"""
        import numpy as np
        from scipy.integrate import solve_bvp
        import matplotlib.pyplot as plt
        import math
        from math import sqrt
        Re_s = $Re_s
        Ts = $Ts
        def oneDiskODE(z,y):
                # Y0 = H, Y1 = F,Y2 = F', Y3 = F'', Y4 = G, Y5 = G'
                dH = -2 * sqrt(Re_s) * y[1]
                dydz = np.zeros((6, len(z)))
                dydz = np.array([dH , y[2] , y[3] , Re_s * ( (1/sqrt(Re_s)) * ( y[3] * y[0] + y[2] * dH) - 2 * y[4] * y[5] + 2 * y[1] * y[2]) , y[5] , Re_s * ( (1/sqrt(Re_s)) * y[5] * y[0] + 2 * y[1] * y[4])])
                return dydz 

        def oneDiskBC(ya, yb):
                resa = np.array([ya[0] + Ts,ya[1], ya[4] - 1])
                
                resb = np.array([yb[0],yb[1], yb[4] - 0])
                
                return np.concatenate((resa, resb))


        z = np.linspace(0, 1, 1000)
        y = np.zeros((6, len(z)))
        y_guess = np.zeros((6, z.size))
        y_guess[0] = 1
        y_guess[1] = 0
        y_guess[2] = 0
        y_guess[3] = 0
        y_guess[4] = 1
        y_guess[5] = 0
        solution = solve_bvp(oneDiskODE, oneDiskBC, z, y_guess,max_nodes=5000000)

        x_plot = np.linspace(0, 1, 2000)


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
    elseif mode == 2
                py"""
        import numpy as np
        from scipy.integrate import solve_bvp
        import matplotlib.pyplot as plt
        kappa = $Ro
        Ts = $Ts
        def oneDiskODE(z, y):
        
                # Y0 = H, Y1 = F', Y2 = F, Y3 = G', Y4 = G
                dydz = np.zeros((5, len(z)))
                dydz = np.array([-2*y[2], y[2] * y[2] - y[4] * y[4] + y[0] * y[1], y[1], 2 * y[2] * y[4] + y[3] * y[0], y[3]])
                return dydz 
        
        def oneDiskBC(ya, yb):
                resa = np.array([ya[0]+Ts,
                                ya[2],
                                ya[4]-1.0])
                
                resb = np.array([yb[2],
                                yb[4]])
                
                return np.concatenate((resa, resb))
        
        
        z = np.linspace(0, 30, 2000)
        y = np.zeros((5, len(z)))
        y_guess = np.zeros((5, z.size))
        y_guess[0] = 1.2
        y_guess[1] = 0
        y_guess[2] = 0
        y_guess[3] = 0
        y_guess[4] = 0
        solution = solve_bvp(oneDiskODE, oneDiskBC, z, y_guess,tol=1e-10,max_nodes=5000000)
        
        x_plot = np.linspace(0, 30, 2000)
        
        
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
    elseif mode == 3 
    py"""
        import numpy as np
        from scipy.integrate import solve_bvp
        kappa = $Ro
        Ts = $Ts
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
                resa = np.array([ya[0]+Ts,
                                ya[2],
                                ya[4]])
                
                resb = np.array([yb[2],
                                yb[4] - 1.0])
                
                return np.concatenate((resa, resb))
        
        
        z = np.linspace(0, 30, 2000)
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
        
        x_plot = np.linspace(0, 30, 2000)
        
        
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
    end
    return u0,v0,w0,du0,dv0,x
 end
 function Cheb(N,mode)
    if mode == 1
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
    else
        θ = range(0,length=N+1,stop=pi)
        x = reshape(-cos.(θ), N+1, 1)
        c = [2; ones(N-1, 1) ; 2] .* (-1) .^ (0:N)
        X = repeat(x, 1, N+1);
        dX = X - X';
        D = (c * (1 ./ c)') ./ (dX .+ I(N+1));
        D = D - diagm(vec(sum(D, dims=2))); 
        a = 2
        b = 0.6
        c = 0.5
        for i=1:N+1
            D[i,:]=D[i,:].* (1-b*x[i]-(1-b)*(x[i]^3+c*(1-x[i]^2)))^2/(2a*(b .+ 3 * (1-b)*x[i]^2 - 2 * c * (1-b) * x[i]))
        end
        for i=1:N+1
            x[i] = a * (1+b*x[i]+(1-b)*(x[i]^3+c*(1-x[i]^2)))/(1-b*x[i]-(1-b)*(x[i]^3+c*(1-x[i]^2)))
            if x[i] > 30
                x[i] = 30
            end
        end
        D2 = D^2;
    end
    return D,D2,x
 end
 function interp(u0,v0,w0,x,N,mode)
        F = Base.zeros(N+1,1)
        G = Base.zeros(N+1,1)
        H = Base.zeros(N+1,1)
        T = Base.zeros(N+1,1)
        if mode ==1 
            z = range(0,1,2000)
            itu = BSplineKit.interpolate(z, u0 , BSplineOrder(4))
            itv = BSplineKit.interpolate(z, v0 , BSplineOrder(4))
            itw = BSplineKit.interpolate(z, w0 , BSplineOrder(4))
            for i = 1 : N + 1
                F[i,1] = itu(x[i])
                G[i,1] = itv(x[i])
                H[i,1] = itw(x[i])
            end
        elseif mode == 2
            z = range(0,30,2000)
            itu = BSplineKit.interpolate(z, u0 , BSplineOrder(4))
            itv = BSplineKit.interpolate(z, v0 , BSplineOrder(4))
            itw = BSplineKit.interpolate(z, w0 , BSplineOrder(4))
            for i = 1 : N + 1
                F[i,1] = itu(x[i])
                G[i,1] = itv(x[i])
                H[i,1] = itw(x[i])
            end
        elseif mode == 3
            z = range(0,30,2000)
            itu = BSplineKit.interpolate(z, u0 , BSplineOrder(4))
            itv = BSplineKit.interpolate(z, v0 , BSplineOrder(4))
            itw = BSplineKit.interpolate(z, w0 , BSplineOrder(4))
            for i = 1 : N + 1
                F[i,1] = itu(x[i])
                G[i,1] = itv(x[i])
                H[i,1] = itw(x[i])
            end
        end
    return F,G,H
 end
end
println("mode = 1:cavity; mode = 2:stationary; mode = 3:rotation;")