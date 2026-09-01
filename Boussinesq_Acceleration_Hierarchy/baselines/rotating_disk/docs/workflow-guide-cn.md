# 基本流与中性曲线计算使用说明

本文说明如何在 `Vonkarmen_bone` 工作区中生成基本流、中性曲线及其汇总数据。
所有命令默认从以下目录运行：

```bash
cd /home/zhj/Rotating-Flow-ToolKit/compress/BEK/Vonkarmen_bone
```

核心 Julia 实现在 `/home/zhj/Rotating-Flow-ToolKit/RotatingDiskFlow` 中，工作区
顶层的 `LopezBaseflow.jl`、`NeutralCurveRunner.jl` 和
`SutherlandMarching.jl` 是兼容入口。

## 1. 环境准备

首次使用工作区环境时执行：

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
python3 -m pip install numpy scipy matplotlib pandas
```

检查 Julia 和 Python：

```bash
julia --version
python3 --version
```

Lopez 基本流通过 `PyCall` 调用 SciPy。如果 Julia 无法导入 SciPy，在 Julia 中将
`PyCall` 指向包含 SciPy 的 Python，然后重新构建：

```julia
using Pkg
ENV["PYTHON"] = Sys.which("python3")
Pkg.build("PyCall")
```

重新构建后应重启 Julia 或 Jupyter 内核。

## 2. 生成基本流

### 2.1 单个 Lopez 广义 Boussinesq 基本流

在 `Bounssinesq.ipynb` 单元格或 Julia REPL 中运行：

```julia
include("LopezBaseflow.jl")
using .LopezBaseflow
using Plots

Tw = 1.10
profile = LopezBaseflow.solve_lopez_baseflow(
    Tw;
    zmax=40.0,
    tol=1.0e-8,
    nout=2001,
    save=false,
)

plot(profile.F, profile.z; label="F", xlabel="velocity", ylabel="z")
plot!(profile.G, profile.z; label="G")
plot!(profile.H, profile.z; label="H")
```

返回的主要变量为：

```text
profile.z       轴向坐标
profile.F       径向速度相似函数
profile.G       周向速度亏损，G(0)=0，G(infinity)=1
profile.H       轴向速度相似函数
profile.T       温度
profile.Fp      dF/dz
profile.Gp      dG/dz
profile.Tp      dT/dz
profile.metrics 峰值、壁面导数、边界残差等信息
```

需要保存时：

```julia
LopezBaseflow.write_profile(
    "lopez_boussinesq_comparison/lopez_Tw1p10.csv",
    profile,
)
```

也可以在求解时设置 `save=true` 和 `output_path=...`。

### 2.2 批量比较传统与 Lopez Boussinesq 基本流

```bash
python3 scripts/compare_lopez_boussinesq.py
```

脚本从 `Tw=1` 的经典 von Karman 解开始进行壁温延拓，输出到：

```text
lopez_boussinesq_comparison/
```

主要文件包括：

```text
profiles_traditional.csv       传统 Boussinesq 基本流
profiles_lopez.csv             Lopez 基本流
baseflow_metrics.csv           峰值、壁面导数和 H_infinity
profile_differences.csv        两种模型的剖面误差
lopez_domain_convergence.csv   不同 z_max 的远场收敛结果
```

该脚本中的壁温序列目前在 `main()` 内定义。改变批量温度范围时，应修改
`common` 和 `generalized` 两个数组。

### 2.3 生成 Boussinesq 与完全可压缩基本流对比数据

Linux/WSL 终端示例：

```bash
TW_VALUES=1.0:0.1:1.4 \
MR=0.3 \
RO=-1.0 \
N_CHEB=199 \
OUT_DIR=baseflow_comparison_data \
python3 scripts/generate_boussinesq_compressible_report.py
```

`TW_VALUES` 支持两种形式：

```text
1.0:0.1:1.4        起点:步长:终点
1.0,1.05,1.1,1.2   显式温度列表
```

PowerShell 中可使用：

```powershell
$env:TW_VALUES = "1.0:0.1:1.4"
$env:MR = "0.3"
$env:OUT_DIR = "baseflow_comparison_data"
python scripts/generate_boussinesq_compressible_report.py
```

默认输出目录为 `baseflow_comparison_data/`，其中包含物理坐标原始剖面、共同网格
插值、Chebyshev 网格插值、导数、特征位置和误差汇总。

### 2.4 指定 `Tw`、`Mr` 和 `R` 的 Sutherland 径向推进基本流

```bash
julia --project=. SutherlandMarching.jl \
  --Tw=1.5 \
  --Mr=0.3 \
  --R=20 \
  --r0=1 \
  --Nr=81 \
  --Nz=121 \
  --zmax=25 \
  --out_dir=sutherland_Tw1.5_Mr0.3_R20
```

这里 `Mr` 是目标位置 `R` 处的局部马赫数，程序内部使用

```math
M_\infty=\frac{M_r(R)}{R},\qquad M_r(r)=M_\infty r.
```

输出包括：

```text
sutherland_profile_at_R.csv       指定 R 位置的基本流
sutherland_marching_profiles.csv  所有径向站位的基本流
sutherland_marching_summary.csv   每个站位的速度峰值和位置
metadata.txt                      参数、速度定义和离散方式
```

速度定义为：

```math
u=rU,\qquad v_{rot}=rV,\qquad v_{inertial}=r(V+1),\qquad w=W.
```

增加 `Nr` 用于检查径向 BDF2 收敛，增加 `Nz` 和 `zmax` 用于检查轴向网格与远场
边界。三者应分别做网格无关性验证。

## 3. 生成中性曲线

### 3.1 单个 Lopez 中性曲线

在 notebook 或 Julia 脚本中运行：

```julia
include("NeutralCurveRunner.jl")
using .NeutralCurveRunner

config = CurveConfig(
    Tw=1.10,
    omega=0.0,
    R_initial=500.0,
    beta_initial=0.04,
    alpha_target=0.10,
    num_modes=2,
    model=:lopez,
    Mr=0.3,
    Ro=-1.0,
    N_cheb=69,
    beta_step=8.0e-4,
    neutral_tol=1.0e-7,
    output_dir=joinpath(pwd(), "neutral_curve_batch"),
    keep_logs=false,
)

result = NeutralCurveRunner.run_case_with_retries(config)
println(result.path)
```

`run_case_with_retries` 会搜索初始中性点、用特征向量重叠追踪同一模态、自动减小
`beta` 步长，并在主分支失败时尝试 Type-I 搜索。研究新分支结构时，可改用只进行
一次延拓的低层接口：

```julia
result = compute_neutral_curve(config)
```

### 3.2 单个完全可压缩中性曲线

只需将配置改为：

```julia
config = CurveConfig(
    Tw=1.10,
    omega=0.0,
    R_initial=500.0,
    beta_initial=0.04,
    alpha_target=0.10,
    num_modes=2,
    model=:compressible,
    Mr=0.3,
    Ro=-1.0,
    N_cheb=69,
    beta_step=8.0e-4,
    neutral_tol=1.0e-7,
    property_perturbations=true,
    base_property_variation=true,
    output_dir=joinpath(pwd(), "neutral_curve_batch"),
    keep_logs=false,
)

result = NeutralCurveRunner.run_case_with_retries(config)
```

物性开关对应关系为：

| 工况 | `property_perturbations` | `base_property_variation` |
|---|---:|---:|
| 完全可压缩 | `true` | `true` |
| 关闭扰动物性 | `false` | `true` |
| 关闭扰动物性并冻结基本流物性 | `false` | `false` |

这两个开关只影响 `model=:compressible`。

### 3.3 标准温度序列批量计算

标准批次计算：

```math
T_w=1.06,1.08,\ldots,1.20,
```

每个温度包含 Lopez 和三种可压缩物性配置，共 32 个工况。

只计算 Lopez，串行执行：

```bash
NEUTRAL_PARALLEL=false \
NEUTRAL_WORKER=lopez \
NEUTRAL_OUTPUT_DIR=neutral_curve_batch \
NEUTRAL_N_CHEB=69 \
NEUTRAL_RESUME=true \
julia --project=. NeutralCurveRunner.jl
```

计算全部四类工况并使用四个进程：

```bash
NEUTRAL_PARALLEL=true \
NEUTRAL_OUTPUT_DIR=neutral_curve_batch \
NEUTRAL_N_CHEB=69 \
NEUTRAL_RESUME=true \
julia --project=. NeutralCurveRunner.jl
```

`NEUTRAL_WORKER` 可取：

```text
lopez
compressible_full
compressible_no_perturb
compressible_frozen
all
```

建议保持 `NEUTRAL_RESUME=true`。程序会重新验证已有结果，只跳过完整工况。

### 3.4 基准温度与特殊分支脚本

重新计算 `Tw=1` 的 Lopez 和可压缩基准曲线：

```bash
julia --project=. scripts/RecomputeBenchmarkCurves.jl
```

只检查 Malik Type-I 和 Type-II 基准点处的特征值：

```bash
julia --project=. scripts/CheckMalikBenchmarks.jl
```

Lopez Type-I 专用脚本目前内置了 `Tw=1.08` 的种子：

```bash
julia --project=. scripts/ComputeLopezTypeI.jl 1.08
```

Type-II 专用脚本目前内置 `Tw=1.08、1.18、1.20` 的种子：

```bash
julia --project=. scripts/ComputeLopezTypeII.jl 1.18 1.20
```

如果要计算其他壁温，应先通过模态扫描获得合适的 `beta` 和 `alpha`，再向脚本的
`DEFAULT_TYPEI_SEEDS` 或 `DEFAULT_TYPEII_SEEDS` 添加种子。不要直接复用相距很远的
壁温种子。

## 4. 合并、验证和分析中性曲线

### 4.1 合并 Lopez Type-I 与 Type-II 分支

当主曲线和 `_branch=typeII.dat` 都存在时：

```bash
julia --project=. scripts/MergeLopezNeutralBranches.jl 1.08 1.18 1.20
```

输出仍位于 `neutral_curve_batch/`，文件名带 `_allbranches.dat`。

### 4.2 生成所有模型的统一 Tecplot 数据

```bash
julia --project=. scripts/IntegrateNeutralCurves.jl
```

脚本读取 `neutral_curve_batch/`，严格检查分支完整性和中性残差，并写入：

```text
neutral_curve_integrated/neutral_curves_all.dat
neutral_curve_integrated/neutral_curves_complete.dat
neutral_curve_integrated/manifest.tsv
neutral_curve_integrated/completeness_manifest.tsv
```

如果任何工况缺少 Type-I 或 Type-II，脚本会拒绝发布汇总文件并报告不完整工况。

### 4.3 提取临界点与模型误差

```bash
julia --project=. scripts/AnalyzeNeutralCriticalErrors.jl
```

输出：

```text
neutral_curve_integrated/neutral_critical_point_errors.tsv
neutral_curve_integrated/neutral_critical_point_errors.dat
neutral_curve_integrated/neutral_critical_point_errors.md
```

进行 Lopez、完全可压缩及物性开关的曲线区间比较：

```bash
julia --project=. scripts/AnalyzeBoussinesqRange.jl
```

## 5. 中性曲线输出格式与绘图

每个 `.dat` 文件包含七列：

```text
omega  R  beta  alpha_r_1  alpha_i_1  alpha_r_2  alpha_i_2
```

空间中性条件为：

```math
\alpha_i=0.
```

Julia 绘图示例：

```julia
using Plots

data = result.data
plot(
    data[:, 2], data[:, 3];
    xlabel="R",
    ylabel="beta",
    linewidth=2,
    label="Tw=1.10",
)

residual = min.(abs.(data[:, 5]), abs.(data[:, 7]))
plot(
    data[:, 3], residual;
    yscale=:log10,
    xlabel="beta",
    ylabel="neutral residual",
    label=false,
)
```

正式使用结果前至少检查：

1. `R-beta` 曲线连续，没有异常长跳跃。
2. 被追踪特征值的 `|alpha_i|` 满足 `neutral_tol`。
3. Type-I 与 Type-II 的临界点均为曲线内部极小值，而不是文件端点。
4. 关键工况在多个 `N_cheb` 下收敛。
5. 基本流和特征函数在 `zmax` 附近已经进入远场状态。

网格验证脚本为：

```bash
julia --project=. scripts/LopezGridIndependence.jl
```

## 6. 推荐计算顺序

1. 使用 `LopezBaseflow.solve_lopez_baseflow` 或 `SutherlandMarching.jl` 检查单个基本流。
2. 检查壁面边界条件、远场条件、速度峰值和温度范围。
3. 运行 `CheckMalikBenchmarks.jl` 验证稳定性入口。
4. 用 `CurveConfig` 计算一个壁温下的一条中性曲线并检查残差。
5. 确认模态追踪正确后再启动标准批次。
6. 补全独立 Type-II 分支，然后运行 `IntegrateNeutralCurves.jl`。
7. 最后提取临界点、计算模型误差并导入 Tecplot。

## 7. 常见问题

### Jupyter 中出现旧的 `eigsol` 方法或关键字不匹配

通常是同一内核多次 `include` 后保留了旧方法。重启内核，并且每个兼容入口只加载
一次。

### 初始位置附近找不到中性点

依次检查 `R_initial`、`beta_initial` 和 `alpha_target`。其中 `alpha_target` 用于从少量
候选特征值中选择初始模态，不是最终中性曲线上的固定波数。

### Type-II 分支没有出现在主曲线中

Type-II 可能与 Type-I 不连通，或者需要低 `beta` 种子。使用
`ComputeLopezTypeII.jl` 单独计算后再合并，不能仅凭主文件缺失就判断模态不存在。

### 配置点越多结果反而变差

高阶微分矩阵的条件数随 `N_cheb` 增长。常规生产计算使用 `N_cheb=69`，通过
`49、59、69、79、89` 等网格序列验证收敛，而不是盲目增加配置点。

### 批量计算被中断

保持 `NEUTRAL_RESUME=true` 后重新运行相同命令。程序会验证已有文件并从未完成工况
继续，不需要删除完整结果。
