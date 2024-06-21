# Module2.jl

module Module2
using .Module1  # 导入 Module1 模块

function call_function_from_module1()
    println("Calling a function from Module1:")
    Module1.function_from_module1()  # 调用 Module1 中的函数
end

end