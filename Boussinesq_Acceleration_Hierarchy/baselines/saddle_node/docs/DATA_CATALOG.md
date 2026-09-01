# Data catalog

## 1. von Kármán

### `finite_fold_analysis/`

传统有限域主分支、第一鞍结、多解剖面及早期相似性时间谱。关键文件包括：

- `branch_by_Hinf.csv`：以 $H_\infty$ 参数化的分支；
- `turning_points.csv`：有限域折返点；
- `domain_sensitivity.csv`：主折叠的域敏感性；
- `upper_branch_domain_stability.csv`：有限域上支时间谱；
- `report.md`：早期分析报告。

### `finite_domain_branches/`

$\eta_{\max}=15$--$80$ 的系统扫描。`turning_points_by_zmax.csv` 清楚显示第一折返点快速收敛，而第二有限域折返点明显漂移。

### `infinite_mapping/`

有理 Chebyshev 映射

$$
\eta=L\frac{1-x}{1+x}
$$

得到的半无限域结果。

推荐入口：

- `INFINITE_MAPPING_UPPER_BRANCH_REPORT.md`：详细报告；
- `confirmed_fold_summary.csv`：推荐鞍结数值；
- `upper_branch_stability_summary.csv`：关键分支位置的时间谱；
- `upper_branch_infinite_mapping_summary.png`：汇总图；
- `stability_N110_L8/mapped_temporal_samples.csv`：不同时间谱阶数的原始特征值。

子目录状态：

| 目录 | 用途 | 状态 |
|---|---|---|
| `convergence_full/` | $N=80,110,140$ 和 $L=4,8,12$ 主映射扫描 | 主数据 |
| `convergence_strict_second_fold/` | 较低条件数、严格残差的第二鞍结验证 | 主数据 |
| `convergence_tight_high_order/` | 较高阶严格残差补充 | 部分完成，但有效文件保留 |
| `stability_N110_L8/` | 上支时间谱 | 主数据 |
| `near_zero_N110_L8/` | $H_\infty\to0^-$ 探索 | 探索性 |
| `prototype/`, `preview_upper_limit/` | 早期试算 | 探索性 |
| `convergence_strict/`, `convergence_N140_L8/` | 中止试算留下的空目录 | 不用于结论 |

## 2. rotor--stator

### `boussinesq_singularity_convergence/`

$\mathrm{Re}_h=1000$ 的精细分支和收敛临界值，是双盘临界参数的推荐数据源：

$$
T_{w,c}^{\mathrm{lower}}=1.155517549054,
\qquad
T_{w,c}^{\mathrm{principal}}=1.167676484047.
$$

### `boussinesq_singularity_results/`

较早但完整的 $\mathrm{Re}_h=1000$ 分支，同时保留 Soong rotating-force 替代模型检查。与 convergence 目录数值非常接近，但文件并非逐字相同，因此两套均保留。

### `boussinesq_singularity_re_scan/`

- $\mathrm{Re}_h=400$：当前扫描区间没有检测到折叠；
- $\mathrm{Re}_h=800$：检测到两个折叠；
- 每个 Reynolds 数均保存分支、验证 JSON、三解样本和图。

`reynolds_scan_summary.csv` 是跨 Reynolds 数的紧凑汇总。

### `three_solution_dynamics/`

$T_w=1.160$ 三组双盘相似性解的：

- 中面空间 Jacobian；
- 时间特征值阶数收敛；
- 主导时间模态；
- 分支稳定性交换；
- 初学者动力系统说明。

### `rotor_stator/tecplot/`

`rotor_stator/data/` 内全部 CSV、JSON 数据和
`rotor_stator/reference/baseflow_Res1000.npz` 的 Tecplot ASCII DAT 派生版本。
目录结构与源数据对应；含 `README.md` 文件索引和
`conversion_manifest.json` 机器可读转换清单。多分支数据按 branch 拆分为多个
ZONE，文本标签保存在 ZONE AUXDATA 中。

## 3. cross-model

`dynamical_singularity_comparison/` 统一比较：

- von Kármán 第一鞍结；
- rotor--stator 两个鞍结；
- 时间零模态和稳态零奇异向量；
- 非退化系数；
- 平方根律；
- 两模型局部普适类与不同反馈机制。

推荐先读 `BEGINNER_DYNAMICAL_SYSTEMS_COMPARISON_REPORT.md`，再查看相应 CSV 和 JSON。

## 4. 使用限制

1. “稳定”只表示轴对称相似性扰动稳定。
2. `finite_domain_branches` 的远上支不能替代真无穷远映射结果。
3. 第三个 von Kármán 折返目前仍为候选结构。
4. 传统 Boussinesq 方程内存在稳定解，不等于该闭合在大温差下定量准确。
