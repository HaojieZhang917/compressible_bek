using LinearAlgebra
using BSplineKit
using PyCall
using DelimitedFiles
using Plots
using NonlinearEigenproblems
include("BaseFlow_cavity.jl")
include("Stability_Cavity.jl")
R_ini= 450
R_end = 85
be_up = 0.3
be_down = -0.3
alpha_ini = 0.005
alpha_end = 0.705
N_cheb = 129
LinearAlgebra.BLAS.set_num_threads(1)
function EigenCore(cof,D,D2,be,alpha,R,N_cheb)
    H0,H1 = CRC_STA.assemble_time_mat(cof,D,D2,be,alpha,R,N_cheb)
    C = eigen(H0,H1)
    val = C.values
    vec = C.vectors
    map_index0 = map(x-> abs(real(x)) < 0.3 && abs(imag(x)) < 0.02, val)
    val_filter0 = val[map_index0]
    vec_filter0 = vec[:,map_index0]
    indictor = sum(real(val_filter0))/length(val_filter0)
    if alpha < 0.031
        addition = 0.003
    else
        addition = 0.007
    end
    map_index = map(x-> (real(x)) > abs(indictor) + addition , val_filter0)
    val_filter1 = val_filter0[map_index]
    vec_filter1 = vec_filter0[:,map_index]
    val_filter2 = val_filter1[findmax(imag.(val_filter1))[2]]
    vec_filter2 = vec_filter1[:,findmax(imag.(val_filter1))[2]];
    if 0.2 < alpha < 0.4
        idx = partialsortperm(real(val_filter1), 1:3, by = x -> abs(x - real(val_filter2)))
        val_filter = val_filter1[idx]
        vec_filter = vec_filter1[:,idx]
    else
        val_filter = val_filter2
        vec_filter = vec_filter2
    end
    val_target = val_filter
    vec_target = vec_filter
    return val_target,vec_target
end

function interation(R_ini, R_end, alpha_ini, alpha_end, be_up, be_down, N_cheb)
    u0,v0,w0,du0,dv0,x = CRC_BF.BaseFlow(1000,-1,0.0,1)
    D,D2,z = CRC_BF.Cheb(N_cheb,1)
    F,G,H = CRC_BF.interp(u0,v0,w0,z,N_cheb,1)
    be_pos_range = 0.005 : 0.005 : be_up
    be_neg_range = -0.005 : -0.005 : be_down
    be_pos_scan_range = 0.005 : 0.01 : be_up
    be_neg_scan_range = -0.005 : -0.01 : be_down
    N_total_steps = 1 + length(be_pos_range) + length(be_neg_range)

    idx_root = length(be_neg_range) + 1
    alpha_range = collect(alpha_ini : 0.01 : alpha_end)

    for R = R_ini : -5 : R_end
        t_start = time()
        println("正在计算雷诺数 = $R ...")
        cof = CRC_STA.Spatial_mode_BEK1((F),(G.-1),(H),R,N_cheb,D,D2, 1000)

        # ====================================================
        # 时间模式矩阵算子预处理 (全局极速组装准备)
        # ====================================================
        N_sys = size(cof.D1, 1)
        idx_bc = setdiff(1:N_sys, (1, N_cheb + 1, N_cheb + 2, 2N_cheb + 2, 2N_cheb + 3, 3N_cheb + 3))
        KD  = kron(I(4), D)
        KD2 = kron(I(4), D2)

        M_const   = cof.D1 .+ cof.C * KD .+ cof.Vzz * KD2
        M_alpha   = im .* cof.A .+ im .* cof.Vxz * KD
        M_alpha2  = -cof.Vxx
        M_be      = im * R .* cof.B .+ im * R .* cof.Vyz * KD
        M_be2     = -R^2 .* cof.Vyy
        M_alphabe = -R .* cof.Vxy

        sM_const   = M_const[idx_bc, idx_bc]
        sM_alpha   = M_alpha[idx_bc, idx_bc]
        sM_alpha2  = M_alpha2[idx_bc, idx_bc]
        sM_be      = M_be[idx_bc, idx_bc]
        sM_be2     = M_be2[idx_bc, idx_bc]
        sM_alphabe = M_alphabe[idx_bc, idx_bc]

        sH1 = (im .* cof.Ta)[idx_bc, idx_bc]

        # ====================================================
        # 阶段一：动态大步长探测 (Probing Phase)
        # ====================================================
        println(">>> [阶段一] 开始动态探测失稳边界...")
        probe_start = 0.3
        probe_step = 0.1
        alpha_cutoff = alpha_end + 1.0
        
        for alpha_probe in probe_start : probe_step : alpha_end - probe_step
            is_globally_stable = true
            
            # 【优化】：预估阶段也使用极速多项式组装
            H0_base_probe = sM_const .+ alpha_probe .* sM_alpha .+ (alpha_probe^2) .* sM_alpha2
            H0_be_linear_probe = sM_be .+ alpha_probe .* sM_alphabe
            H0_work_probe = zeros(ComplexF64, size(sH1))
            H1_work_probe = copy(sH1)

            vals_root, vecs_root = EigenCore(cof, D, D2, 0.0, alpha_probe, R, N_cheb)
            
            # 只要任意一个模态不稳定，就标记
            if any(imag.(vals_root) .> 0)
                is_globally_stable = false
            end
            
            if is_globally_stable
                # 遍历所有危险模态进行探测
                for m in 1:length(vals_root)
                    val, vec = vals_root[m], vecs_root[:, m]
                    for be in be_pos_scan_range
                        H0_work_probe .= H0_base_probe .+ be .* H0_be_linear_probe .+ (be^2) .* sM_be2
                        val, vec = rayleigh_quotient_iteration(H0_work_probe, H1_work_probe, val, vec)
                        if imag(val) > 0
                            is_globally_stable = false
                            break 
                        end
                    end
                    if !is_globally_stable; break; end
                end
            end

            if is_globally_stable
                for m in 1:length(vals_root)
                    val, vec = vals_root[m], vecs_root[:, m]
                    for be in be_neg_scan_range
                        H0_work_probe .= H0_base_probe .+ be .* H0_be_linear_probe .+ (be^2) .* sM_be2
                        val, vec = rayleigh_quotient_iteration(H0_work_probe, H1_work_probe, val, vec)
                        if imag(val) > 0
                            is_globally_stable = false
                            break
                        end
                    end
                    if !is_globally_stable; break; end
                end
            end
            
            if is_globally_stable
                alpha_cutoff = alpha_probe
                println("  [!] 探测结果：alpha = $alpha_cutoff 及其后续所有波均已稳定！截断点已锁定。")
                break
            else
                println("  [-] 探测结果：alpha = $alpha_probe 仍存在不稳定波，继续探测...")
            end
        end

        # ====================================================
        # 阶段二：并行主计算程序 (Main Parallel Execution)
        # ====================================================
        println(">>> [阶段二] 启动高精度并行扫描，截断阈值 alpha_cutoff = $alpha_cutoff ...")
        tasks = Vector{Task}(undef, length(alpha_range))

        for i in 1 : length(alpha_range)
            alpha = alpha_range[i]

            tasks[i] = Threads.@spawn begin
                # 【多模态核心】：初始化包络线矩阵
                local_mat_env = zeros(Float64, N_total_steps, 5)
                for r_idx in 1:N_total_steps
                    b_val = (r_idx == idx_root) ? 0.0 : 
                            (r_idx > idx_root ? be_pos_range[r_idx - idx_root] : be_neg_range[idx_root - r_idx])
                    local_mat_env[r_idx, :] .= (R, alpha, b_val, -1.0, -1.0)
                end

                if alpha >= alpha_cutoff
                    return local_mat_env # 剪枝直接返回包络线
                end

                H0_base = sM_const .+ alpha .* sM_alpha .+ (alpha^2) .* sM_alpha2
                H0_be_linear = sM_be .+ alpha .* sM_alphabe
                H0_work = zeros(ComplexF64, size(sH1))
                H1_work = copy(sH1)

                vals_root, vecs_root = EigenCore(cof,D,D2,0.0,alpha,R,N_cheb)
                
                # ====================================================
                # 分别追踪提取出的前几个最危险模态
                # ====================================================
                for m in 1:length(vals_root)
                    val_root = vals_root[m]
                    vec_root = vecs_root[:, m]
                    
                    # 临时矩阵用于记录当前单一模态的轨迹
                    local_mat = copy(local_mat_env)
                    local_mat[idx_root, :] .= (R, alpha, 0.0, real(val_root), imag(val_root))
                    
                    val, vec = val_root, vec_root
                    crossed_neutral = false 
                    post_steps = 0
                    
                    # --- 正向追踪 ---
                    for (ib, be) in enumerate(be_pos_range)
                        H0_work .= H0_base .+ be .* H0_be_linear .+ (be^2) .* sM_be2
                        if  300 < R < 400 && 0.2 < alpha < 0.4
                        val, vec = rayleigh_quotient_iteration(H0_work, H1_work, val + 0.03 * alpha * im, vec)
                        else
                        val, vec = rayleigh_quotient_iteration(H0_work, H1_work, val, vec)
                        end
                        curr_idx = idx_root + ib
                        local_mat[curr_idx, :] .= (R, alpha, be, real(val), imag(val))

                        if !crossed_neutral && local_mat[curr_idx, 5] < -1e-4 && local_mat[curr_idx - 1, 5] > -1e-4 && abs(be) > 0.1
                            crossed_neutral = true
                        end

                        if crossed_neutral
                            if post_steps >= 10
                                if curr_idx < N_total_steps
                                    for remain_idx in (curr_idx + 1) : N_total_steps
                                        local_mat[remain_idx, :] .= (R, alpha, be_pos_range[remain_idx - idx_root], local_mat[curr_idx, 4], local_mat[curr_idx, 5])
                                    end
                                end
                                break 
                            else
                                post_steps += 1 
                            end
                        end
                    end
                    
                    # --- 负向追踪 ---
                    val, vec = val_root, vec_root
                    crossed_neutral = false 
                    post_steps = 0
                    for (ib, be) in enumerate(be_neg_range)
                        H0_work .= H0_base .+ be .* H0_be_linear .+ (be^2) .* sM_be2
                        val, vec = rayleigh_quotient_iteration(H0_work, H1_work, val + val + 0.03 * alpha * im, vec)
                        
                        curr_idx = idx_root - ib
                        local_mat[curr_idx, :] .= (R, alpha, be, real(val), imag(val))

                        if !crossed_neutral && local_mat[curr_idx, 5] < -1e-4 && local_mat[curr_idx + 1, 5] > -1e-4 && abs(be) > 0.1
                            crossed_neutral = true
                        end
                        
                        if crossed_neutral
                            if post_steps >= 10
                                if curr_idx > 1
                                    for remain_idx in 1 : (curr_idx - 1)
                                        local_mat[remain_idx, :] .= (R, alpha, be_neg_range[idx_root - remain_idx], local_mat[curr_idx, 4], local_mat[curr_idx, 5])                        
                                    end
                                end
                                break 
                            else
                                post_steps += 1
                            end
                        end
                    end

                    # ====================================================
                    # 【包络线合并】：将当前模态合并入总环境矩阵
                    # ====================================================
                    for r_idx in 1:N_total_steps
                        if local_mat[r_idx, 5] != -1.0
                            # 如果该点还没数据，或者当前模态的虚部(增长率)更大，就覆盖它
                            if local_mat_env[r_idx, 5] == -1.0 || local_mat[r_idx, 5] > local_mat_env[r_idx, 5]
                                local_mat_env[r_idx, 4] = local_mat[r_idx, 4]
                                local_mat_env[r_idx, 5] = local_mat[r_idx, 5]
                            end
                        end
                    end
                end # 单一模态追踪结束

                return local_mat_env
            end
        end
        
        results_mat = [fetch(t) for t in tasks]
        full_R_data = vcat(results_mat...)
        full_R_data = sortslices(full_R_data, dims = 1, by = x->(x[2], x[3]))
        full_R_data = filter_boundary_instability!(full_R_data)
        
        filename = "R=$R.dat"
        open(filename, "w") do io
            println(io, "TITLE = \"Stability Analysis\"")
            println(io, "VARIABLES = \"R\", \"alpha\", \"beta\", \"omega_r\", \"omega_i\"")
            println(io, "ZONE T=\"R=$R\", I=$N_total_steps, J=$(length(alpha_range)), F=POINT")
            writedlm(io, full_R_data)
        end
        t_end = time()
            println(">>> R = $R 计算完成，耗时: $(round(t_end - t_start, digits=2)) 秒")    
        end
end
function rayleigh_quotient_iteration(A, B, sigma, q0=rand(ComplexF64, size(A, 1)))
    tol = 1e-8
    sigma_current = ComplexF64(sigma[1]) 
    q = q0 / norm(q0) # 必须使用欧几里得范数
    
    for i in 1:20
        sigma_old = sigma_current
        
        v = (A - sigma_current * B) \ (B * q)
        q = v / norm(v)
        
        num = dot(q, A * q)
        den = dot(q, B * q)
        sigma_current = num / den
        
        if abs(sigma_current - sigma_old) < tol
            return sigma_current, q
        end
    end
    @warn "RQI failed to converge"
    return sigma_current, q
end
function filter_boundary_instability!(data::Matrix{Float64})
    # 获取所有的 alpha 唯一值 (也就是提取所有不同的波数)
    alphas = unique(data[:, 2])
    
    for alpha in alphas

        idx = findall(x -> x == alpha, data[:, 2])
        
        if isempty(idx)
            continue
        end
        
        # ---------------------------------------------------------
        # 1. 检查上限边界 (即 beta = 0.3，对应这段索引的最后一行)
        # ---------------------------------------------------------
        idx_up = idx[end]
        if data[idx_up, 5] > 0.0 # 如果上边界点是不稳定的
            # 从边界开始，倒推向内部扫描
            for i in length(idx):-1:1
                curr_row = idx[i]
                if data[curr_row, 5] > 0.0
                    # 只要还是不稳定的，就抹除（填充为 -1.0）
                    data[curr_row, 4] = -1.0
                    data[curr_row, 5] = -1.0
                else
                    # 【核心剪枝】：一旦碰到稳定的点，说明这个碰到边界的不稳定分支结束了。
                    # 立刻 break 跳出，保护更内部的其他不稳定模态！
                    break 
                end
            end
        end
        
        # ---------------------------------------------------------
        # 2. 检查下限边界 (即 beta = -0.3，对应这段索引的第一行)
        # ---------------------------------------------------------
        idx_down = idx[1]
        if data[idx_down, 5] > 0.0 # 如果下边界点是不稳定的
            # 从边界开始，正推向内部扫描
            for i in 1:length(idx)
                curr_row = idx[i]
                if data[curr_row, 5] > 0.0
                    data[curr_row, 4] = -1.0
                    data[curr_row, 5] = -1.0
                else
                    # 碰到稳定点，立刻跳出，保护内部模态
                    break 
                end
            end
        end
    end
    
    return data
end
interation(R_ini, R_end, alpha_ini, alpha_end, be_up, be_down, N_cheb)
