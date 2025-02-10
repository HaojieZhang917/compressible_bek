include("CRD_STA.jl")
include("LST_BEK.jl")
using Plots
using LinearAlgebra
using NonlinearEigenproblems
using ProgressMeter
using DelimitedFiles
N_cheb = 99
Mr = 0.3
Tw = 1
gamma = 1.4
sigma = 0.72
omega = 0
global netural_cur = [-1 -1 -1]
global netural_cur_real = [-1 -1 -1]
global netural_cur_imag = [-1 -1 -1]
@showprogress for R = 280 : 1 : 290
    theta = 0.3
    Ma = Mr/R
    u0,v0,w0,f,q,D,D2,x = baseflow_var(N_cheb)
    H,T = T_ca(Mr,f,q,w0,gamma,Tw)
    F,G,H,T,rho,z = interp(u0,v0,H,T,x,N_cheb,"sim")
    lam = - (2/3) * T
    kappa = (1/sigma) * T
    num :: Int64 = 3
    for be = 0.02 : 0.001 : 0.15
        A0,A1,A2 = Spatial_mode(F,G,H,rho,lam,kappa,T,sigma,gamma,R,Ma,omega,be)
        nep = PEP([A0,A1,A2]); 
        eigval,eigvec = iar(nep,σ = theta , neigs = num ,maxit = 500,tol=1e-10)
        point = filter(x -> - 0.002 < imag(x) < 0 && real(x)>0.05, eigval)
        if point == []

            eig = - 1

        elseif netural_cur[end,3] == -1 && imag(point[1]) < -0.0005 || imag(point[1]) > 0.0005

            eig = -1

        else
            
            eig = point[findmin(x-> abs(imag(x)), point)[2]]
            theta = eig
            num = 2

        end
        global netural_cur = [netural_cur; [R be eig]]
        global netural_cur_real = [netural_cur_real; [R be real(eig)]]
        global netural_cur_imag = [netural_cur_imag; [R be imag(eig)]]
    end
    println("R = ",R)
end
writedlm("CRD_STA_test.dat",netural_cur)
writedlm("CRD_STA_test_imag.dat",netural_cur_imag)
writedlm("CRD_STA_test_real.dat",netural_cur_real)