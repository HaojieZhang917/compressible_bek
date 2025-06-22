function eigvector(eigval,values,vec)
    index = findall(x-> x== eigval,values)
    vector = vec[:,index] 
    return vector
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
    使用四阶中心差分法计算向量y的导数。
    
    参数:
    y : 需要求导的函数值向量
    h : 差分步长
    
    返回:
    dy_dx : 与y等长的导数向量
    """
    n = length(y)
    dy_dx = zeros(n)
    
    # 内部点使用四阶中心差分公式
    for i in 3:n-2
        dy_dx[i] = (-y[i+2] + 8*y[i+1] - 8*y[i-1] + y[i-2]) / (12*h)
    end
    
    # 边界点使用二阶中心差分或向前/向后差分
    # 第1点使用向前差分
    dy_dx[1] = (-3*y[1] + 4*y[2] - y[3]) / (2*h)
    # 第2点使用二阶中心差分
    dy_dx[2] = (y[3] - y[1]) / (2*h)
    # 倒数第2点使用二阶中心差分
    dy_dx[n-1] = (y[n] - y[n-2]) / (2*h)
    # 最后一点使用向后差分
    dy_dx[n] = (3*y[n] - 4*y[n-1] + y[n-2]) / (2*h)
    
    return dy_dx
end