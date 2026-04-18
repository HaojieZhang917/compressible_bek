using LinearAlgebra
using BSplineKit
using PyCall
using DelimitedFiles
using Plots
using NonlinearEigenproblems
include("BaseFlow_cavity.jl")
include("Stability_Cavity.jl")

LinearAlgebra.BLAS.set_num_threads(1)
function EigenCore(cof,D,D2,be,alpha,R,N_cheb)
    H0,H1 = CRC_STA.assemble_time_mat(cof,D,D2,be,alpha,R,N_cheb)
    C = eigen(H0,H1)
    val = C.values
    vec = C.vectors
    map_index0 = map(x-> abs(real(x)) < 0.3 && abs(imag(x)) < 0.5, val)
    val_filter0 = val[map_index0]
    vec_filter0 = vec[:,map_index0]
    indictor = sum(real(val_filter0))/length(val_filter0)
    map_index = map(x-> (real(x)) > indictor + 0.005 , val_filter0)
    val_filter = val_filter0[map_index]
    vec_filter = vec_filter0[:,map_index]
    val_target = val_filter[findmax(imag.(val_filter))[2]]
    vec_target = vec_filter[:,findmax(imag.(val_filter))[2]];
    return val_target,vec_target
end

function interation(R_ini, R_end, alpha_ini, alpha_end, be_up, be_down, N_cheb)
    u0,v0,w0,du0,dv0,x = CRC_BF.BaseFlow(1000,-1,0.0,1)
    D,D2,z = CRC_BF.Cheb(N_cheb,1)
    F,G,H = CRC_BF.interp(u0,v0,w0,z,N_cheb,1)
    be_pos_range = 0.005 : 0.005 : be_up
    be_neg_range = -0.005 : -0.005 : be_down
    N_total_steps = 1 + length(be_pos_range) + length(be_neg_range)

    idx_root = length(be_neg_range) + 1

    alpha_range = collect(alpha_ini : 0.025 : alpha_end)

    io_lock = ReentrantLock()

    for R = R_ini : -5 : R_end
        println("正在计算雷诺数 = $R ...")
        open("R = $R.dat ", "w") do io
            cof = CRC_STA.Spatial_mode_BEK1((F),(G.-1),(H),R,N_cheb,D,D2, 1000)

            Threads.@threads for alpha in alpha_range
                
                results_mat = zeros(Float64, N_total_steps, 5)

                val_root,vec_root = EigenCore(cof,D,D2,0.0,alpha,R,N_cheb)
                results_mat[idx_root, :] .= (R, alpha, 0.0, real(val_root), imag(val_root))                
                
                val, vec = val_root, vec_root
                for (i, be) in enumerate(be_pos_range)
                    H0,H1 = CRC_STA.assemble_time_mat(cof,D,D2,be,alpha,R,N_cheb)
                    val,vec = rayleigh_quotient_iteration(H0,H1,val, vec)
                    results_mat[idx_root + i, :] .= (R, alpha, be, real(val), imag(val))              
                end
                val, vec = val_root, vec_root
                for (i, be) in enumerate(be_neg_range)
                    H0,H1 = CRC_STA.assemble_time_mat(cof,D,D2,be,alpha,R,N_cheb)
                    val,vec = rayleigh_quotient_iteration(H0,H1,val, vec)
                    results_mat[idx_root - i, :] .= (R, alpha, be, real(val), imag(val))
                end

                lock(io_lock)
                try
                    writedlm(io, results_mat)
                    flush(io)
                finally
                    unlock(io_lock)
                end
                
            end
        end
    end
end
function rayleigh_quotient_iteration(A, B, sigma, q0=rand(size(A, 1), 1))

    flg = true
    i=1
    while flg
        i=i+1
        sigma0 = real(sigma[1]) + abs(imag(sigma[1]))im + 0.0e0im
        q = (A - sigma*B) \ (B*q0)
        q0 = q/maximum(abs.(q))
        sigma = ((q0'*(A*q0))/(q0'*(B*q0)))[1]
        if abs(sigma-sigma0)<=eps(1.0f0)
            flg = false
        end
        if i==20
            flg=false
        end
    end

      return sigma, q0
end
R_ini= 495
R_end = 90
be_up = 0.3
be_down = -0.3
alpha_ini = 0.001
alpha_end = 0.71
N_cheb = 129
interation(R_ini, R_end, alpha_ini, alpha_end, be_up, be_down, N_cheb)