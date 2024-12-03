using Distributed
addprocs(4)  # 添加 4 个工作进程

@everywhere function power(x, n)
    return x^n
end

results = pmap(x -> power(x, 3), 1:10)
println(results)

