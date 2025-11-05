using DelimitedFiles
using LinearAlgebra
function eigvector(eigval, values, vec)
    index = findall(x-> abs(x - eigval) < 1e-8, values)  # 使用容差比较避免浮点误差
    if isempty(index)
        # 如果找不到精确匹配，返回最近的特征向量
        idx = argmin(abs.(values .- eigval))
        return vec[:, idx]
    else
        return vec[:, index[1]]  # 取第一个匹配的
    end
end
function RQI(A, B, sigma; q0=rand(size(A, 1), 1))

    flg = true
    while flg
        sigma0 = sigma[1]+ 0.0e0im
        q = (A - sigma*B) \ (B*q0)
        q0 = q/maximum(abs.(q))
        sigma = ((q0'*(A*q0))/(q0'*(B*q0)))[1]
        if abs(sigma-sigma0)<=eps(1.0f0)
            flg = false
        end

    end

      return sigma, q0
end


function diff1(y, h)
    """
    使用六阶中心差分法计算向量y的导数。
    
    参数:
    y : 需要求导的函数值向量
    h : 差分步长
    
    返回:
    dy_dx : 与y等长的导数向量
    """
    n = length(y)
    dy_dx = zeros(n)
    
    # 内部点使用六阶中心差分公式
    for i in 4:n-3
        dy_dx[i] = (-1*y[i+3] + 9*y[i+2] - 45*y[i+1] + 45*y[i-1] - 9*y[i-2] + y[i-3]) / (60*h)
    end
    
    # 边界点使用较低阶的中心差分或向前/向后差分
    # 第1点使用向前差分
    dy_dx[1] = (-25*y[1] + 48*y[2] - 36*y[3] + 16*y[4] - 3*y[5]) / (12*h)
    # 第2点使用向前差分
    dy_dx[2] = (-3*y[1] + 27*y[2] - 27*y[3] + 9*y[4] - y[5]) / (12*h)
    # 第3点使用二阶中心差分
    dy_dx[3] = (y[5] - y[1]) / (4*h)
    # 倒数第3点使用二阶中心差分
    dy_dx[n-2] = (y[n] - y[n-4]) / (4*h)
    # 倒数第2点使用向前差分
    dy_dx[n-1] = (y[n] - 27*y[n-1] + 27*y[n-2] - 9*y[n-3] + y[n-4]) / (12*h)
    # 最后一点使用向后差分
    dy_dx[n] = (25*y[n] - 48*y[n-1] + 36*y[n-2] - 16*y[n-3] + 3*y[n-4]) / (12*h)
    
    return dy_dx
end
function extline(para,dat)
    data = readdlm(dat)[3:end,:]
    c = data[:,1] .== para
    if norm(c) == 0
        data_positive = (data[data[:,2] .== para,:])
    else
        data_positive = (data[data[:,1] .== para,:])
    end
    data = data_positive[data_positive[:,4].>0,:]
    return data
end