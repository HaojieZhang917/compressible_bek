# Rotor--stator Tecplot data

本目录由 `convert_all_to_tecplot.py` 从原始 CSV、JSON 和 NPZ 数据只读生成。
所有文件采用 Tecplot ASCII、`DATAPACKING=POINT`。包含 `branch` 的剖面或谱数据按 branch 分为多个 `ZONE`；文本字段写入 zone `AUXDATA`。

| DAT 文件 | 来源 | 行数 | ZONE 数 | 数值变量数 |
|---|---|---:|---:|---:|
| `data/boussinesq_singularity_convergence/branch_Re1000_traditional_centrifugal.dat` | `boussinesq_singularity_convergence/branch_Re1000_traditional_centrifugal.csv` | 1048 | 1 | 10 |
| `data/boussinesq_singularity_convergence/folds_Re1000_traditional_centrifugal.dat` | `boussinesq_singularity_convergence/folds_Re1000_traditional_centrifugal.json` | 2 | 2 | 8 |
| `data/boussinesq_singularity_convergence/three_solutions_Tw1.160.dat` | `boussinesq_singularity_convergence/three_solutions_Tw1.160.csv` | 3003 | 3 | 7 |
| `data/boussinesq_singularity_convergence/validation_Re1000_traditional_centrifugal.dat` | `boussinesq_singularity_convergence/validation_Re1000_traditional_centrifugal.json` | 1 | 1 | 11 |
| `data/boussinesq_singularity_re_scan/Re400/branch_Re400_traditional_centrifugal.dat` | `boussinesq_singularity_re_scan/Re400/branch_Re400_traditional_centrifugal.csv` | 384 | 1 | 10 |
| `data/boussinesq_singularity_re_scan/Re400/folds_Re400_traditional_centrifugal.dat` | `boussinesq_singularity_re_scan/Re400/folds_Re400_traditional_centrifugal.json` | 1 | 1 | 1 |
| `data/boussinesq_singularity_re_scan/Re400/three_solutions_Tw1.160.dat` | `boussinesq_singularity_re_scan/Re400/three_solutions_Tw1.160.csv` | 1001 | 1 | 7 |
| `data/boussinesq_singularity_re_scan/Re400/validation_Re400_traditional_centrifugal.dat` | `boussinesq_singularity_re_scan/Re400/validation_Re400_traditional_centrifugal.json` | 1 | 1 | 5 |
| `data/boussinesq_singularity_re_scan/Re800/branch_Re800_traditional_centrifugal.dat` | `boussinesq_singularity_re_scan/Re800/branch_Re800_traditional_centrifugal.csv` | 387 | 1 | 10 |
| `data/boussinesq_singularity_re_scan/Re800/folds_Re800_traditional_centrifugal.dat` | `boussinesq_singularity_re_scan/Re800/folds_Re800_traditional_centrifugal.json` | 2 | 2 | 8 |
| `data/boussinesq_singularity_re_scan/Re800/three_solutions_Tw1.160.dat` | `boussinesq_singularity_re_scan/Re800/three_solutions_Tw1.160.csv` | 3003 | 3 | 7 |
| `data/boussinesq_singularity_re_scan/Re800/validation_Re800_traditional_centrifugal.dat` | `boussinesq_singularity_re_scan/Re800/validation_Re800_traditional_centrifugal.json` | 1 | 1 | 5 |
| `data/boussinesq_singularity_results/branch_Re1000_traditional_centrifugal.dat` | `boussinesq_singularity_results/branch_Re1000_traditional_centrifugal.csv` | 609 | 1 | 10 |
| `data/boussinesq_singularity_results/folds_Re1000_traditional_centrifugal.dat` | `boussinesq_singularity_results/folds_Re1000_traditional_centrifugal.json` | 2 | 2 | 8 |
| `data/boussinesq_singularity_results/soong_model/branch_Re1000_soong_rotating_forces.dat` | `boussinesq_singularity_results/soong_model/branch_Re1000_soong_rotating_forces.csv` | 101 | 1 | 10 |
| `data/boussinesq_singularity_results/soong_model/folds_Re1000_soong_rotating_forces.dat` | `boussinesq_singularity_results/soong_model/folds_Re1000_soong_rotating_forces.json` | 1 | 1 | 1 |
| `data/boussinesq_singularity_results/soong_model/validation_Re1000_soong_rotating_forces.dat` | `boussinesq_singularity_results/soong_model/validation_Re1000_soong_rotating_forces.json` | 1 | 1 | 11 |
| `data/boussinesq_singularity_results/three_solutions_Tw1.160.dat` | `boussinesq_singularity_results/three_solutions_Tw1.160.csv` | 3003 | 3 | 7 |
| `data/boussinesq_singularity_results/validation_Re1000_traditional_centrifugal.dat` | `boussinesq_singularity_results/validation_Re1000_traditional_centrifugal.json` | 1 | 1 | 11 |
| `data/reynolds_scan_summary.dat` | `reynolds_scan_summary.csv` | 3 | 3 | 4 |
| `data/three_solution_dynamics/branch_diagnostics.dat` | `three_solution_dynamics/branch_diagnostics.csv` | 3 | 3 | 32 |
| `data/three_solution_dynamics/branch_temporal_growth.dat` | `three_solution_dynamics/branch_temporal_growth.csv` | 74 | 1 | 6 |
| `data/three_solution_dynamics/midplane_spatial_eigenvalues.dat` | `three_solution_dynamics/midplane_spatial_eigenvalues.csv` | 21 | 3 | 5 |
| `data/three_solution_dynamics/temporal_convergence.dat` | `three_solution_dynamics/temporal_convergence.json` | 12 | 3 | 5 |
| `data/three_solution_dynamics/temporal_eigenvalues_convergence.dat` | `three_solution_dynamics/temporal_eigenvalues_convergence.csv` | 240 | 3 | 7 |
| `reference/baseflow_Res1000.dat` | `baseflow_Res1000.npz` | 20000 | 1 | 6 |

## 直接计算的基本流剖面

`data/baseflow_profiles/principal_branch_Tw1.00_1.16_step0.04.dat` 包含
等温相连主支在 $T_w=1.00,1.04,1.08,1.12,1.16$ 的五个基本流；每个壁温
对应一个 2001 点 ZONE。变量为 `z,H,F,F_z,G,G_z,T,T_z,pressure_gradient,Tw`。
对应数值校验见同目录的 `_validation.json` 文件，复现脚本为
`rotor_stator/scripts/generate_principal_baseflow_tecplot.jl`。

## Tecplot 导入

在 Tecplot 360 中选择 `File > Load Data File(s)`，数据类型选择 `Tecplot Data Loader`，再选择一个或多个 `.dat` 文件。
缺失数值写成 `NaN`。空的折返点 JSON 会转换成一个 `record_count=0` 的占位 zone。
