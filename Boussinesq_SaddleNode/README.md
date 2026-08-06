# Boussinesq Saddle-Node Project

本项目汇总传统离心 Boussinesq 闭合在两类旋转流相似性模型中产生的鞍结问题：

- 半无限域 von Kármán 旋转盘边界层；
- 有限间隙 rotor--stator 双盘系统；
- 两种模型的零模态、时间稳定性、平方根律和动力系统比较。

这是从原研究目录复制并重新组织得到的独立工作区。原目录
`compress/BEK/Vonkarmen_bone` 和 `TwoDisk` 没有被移动或覆盖；本项目中的脚本路径已经改为只读取本项目内的数据和参考文件。

## 已确认的核心结果

### von Kármán 半无限域

有理 Chebyshev 真无穷远映射确认至少两个鞍结：

$$
\begin{aligned}
H_{\infty,c1}&\approx-0.53276166,
&T_{w,c1}&\approx1.048021731,\\
H_{\infty,c2}&\approx-0.11331,
&T_{w,c2}&\approx1.03014466.
\end{aligned}
$$

第二鞍结后的分支仍有正实时间特征值，所以它是数学上存在但不吸引的相似性平衡态。靠近 $H_\infty=0$ 还观察到第三个候选折返，但该区域处于长热尾奇异极限，目前没有把它列为与前两个同等级的确认结论。

### rotor--stator，$\mathrm{Re}_h=1000$

传统模型存在两个鞍结：

$$
T_{w,c}^{\mathrm{lower}}\approx1.155517549,
\qquad
T_{w,c}^{\mathrm{principal}}\approx1.167676484.
$$

两个外支在轴对称相似性子空间内稳定，中支有一个实不稳定方向，因此模型内存在双稳态和迟滞窗口。

## 目录结构

```text
Boussinesq_SaddleNode/
├── von_karman/
│   ├── scripts/       # 有限域、真无穷远映射和谱分析
│   ├── notebooks/     # 原 Boussinesq 基本流笔记本
│   ├── reference/     # 历史基本流参考数据
│   └── data/
│       ├── finite_fold_analysis/
│       ├── finite_domain_branches/
│       └── infinite_mapping/
├── rotor_stator/
│   ├── scripts/       # 双盘延拓和三解动力学
│   ├── notebooks/     # 原基本流与稳定性笔记本
│   ├── reference/     # 等温参考解与 Julia 方程实现
│   └── data/          # Re 扫描、三解和时间谱
├── cross_model/
│   ├── scripts/       # 两模型鞍结统一比较
│   └── data/          # 零模态、平方根律和比较报告
├── src/               # 共用 Chebyshev、延拓、稳定性和 I/O 模块
├── docs/
│   ├── DATA_CATALOG.md
│   ├── SOURCE_MANIFEST.md
│   └── ROTOR_STATOR_SINGULARITY_REPORT.md
├── Project.toml
├── Manifest.toml
├── AGENTS.md           # 后续新程序统一使用 Julia
├── check_project.py    # 只读项目检查辅助工具
└── requirements.txt    # 既有 Python 数据处理工具依赖
```

## 推荐阅读顺序

1. [`von_karman/data/infinite_mapping/INFINITE_MAPPING_UPPER_BRANCH_REPORT.md`](von_karman/data/infinite_mapping/INFINITE_MAPPING_UPPER_BRANCH_REPORT.md)
2. [`docs/ROTOR_STATOR_SINGULARITY_REPORT.md`](docs/ROTOR_STATOR_SINGULARITY_REPORT.md)
3. [`cross_model/data/dynamical_singularity_comparison/BEGINNER_DYNAMICAL_SYSTEMS_COMPARISON_REPORT.md`](cross_model/data/dynamical_singularity_comparison/BEGINNER_DYNAMICAL_SYSTEMS_COMPARISON_REPORT.md)
4. [`docs/DATA_CATALOG.md`](docs/DATA_CATALOG.md)

## 计算环境

主要数值求解、延拓和稳定性程序已迁移为原生 Julia，不调用 Python、SciPy 或 PythonCall。建议使用 Julia 1.10 或更高版本，在项目根目录执行：

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

已有的 Tecplot 转换、数据校验、结果汇总和项目检查仍保留为 Python 辅助程序；需要运行这些工具时安装：

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

只读完整性检查仍使用原 Python 工具：

```bash
python3 check_project.py
```

## 主要复现入口

### von Kármán 有限域分支

```bash
julia --project=. von_karman/scripts/compute_domain_branches.jl
```

### von Kármán 真无穷远映射和时间谱

```bash
julia --project=. von_karman/scripts/analyze_vonkarman_infinite_mapping.jl \
  --degree 110 --scale 8 --h-stop -0.02 \
  --tolerance 3e-10 --temporal-degrees 50,70,90 \
  --output-dir von_karman/data/infinite_mapping/stability_N110_L8

python3 von_karman/scripts/summarize_vonkarman_infinite_mapping.py
```

### rotor--stator 分支和三解动力学

```bash
julia --project=. rotor_stator/scripts/two_disk_boussinesq_singularity.jl

julia --project=. rotor_stator/scripts/analyze_three_solution_dynamics.jl
```

### 两模型动力系统比较

```bash
julia --project=. cross_model/scripts/compare_saddle_node_dynamics.jl
```

### rotor--stator 的 Tecplot 数据

现有 rotor--stator 数据已经转换为 Tecplot ASCII DAT，统一位于
`rotor_stator/tecplot/`。转换过程不修改原始 CSV、JSON 或 NPZ 文件。

重新生成并校验全部 DAT：

```bash
python3 rotor_stator/scripts/convert_all_to_tecplot.py
python3 rotor_stator/scripts/verify_tecplot_data.py
```

生成等温相连主支在 $T_w=1.00,1.04,1.08,1.12,1.16$ 的基本流剖面：

```bash
julia --project=. rotor_stator/scripts/generate_principal_baseflow_tecplot.jl
```

各文件的变量、分区数、数据行数及源文件映射见
`rotor_stator/tecplot/README.md` 和 `conversion_manifest.json`。

这些命令可能覆盖对应的派生输出。需要保留当前快照时，应通过脚本的输出参数写入单独的试算目录。

## 数据解释原则

- `finite_*` 保存原有限截断域结果，用于说明第二折返点的域敏感性。
- `infinite_mapping/convergence_full` 是映射尺度和谱阶数的主扫描。
- `convergence_strict_second_fold` 与 `convergence_tight_high_order` 用于第二鞍结严格残差交叉验证。
- `prototype`、`preview_upper_limit`、空的严格收敛目录属于探索性试算，不应优先于汇总表。
- 推荐引用数值见 `confirmed_fold_summary.csv` 和各模型的 JSON 汇总，而不是从图上读取。
- 当前“稳定”均只指轴对称相似性子空间，不能替代完整三维稳定性分析。

## 版本控制说明

仓库根目录默认忽略 `*.csv`、`*.dat` 和 `*.lpk`。本项目的 `.gitignore` 对整理后的科学数据作了局部反向声明，使这些文件可以作为项目资产被版本控制；本次整理没有执行暂存或提交。
