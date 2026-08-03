# Source manifest

整理日期：2026-08-02。

本项目采用复制方式建立，原始文件保持不变。下表记录目录级来源；项目中的 Python 脚本只修改了路径常量和默认输入/输出位置，没有改变控制方程、离散方法或物理参数。

| 类别 | 原始路径（相对 `Rotating-Flow-ToolKit`） | 新项目路径 | 处理 |
|---|---|---|---|
| von Kármán 主折叠程序 | `compress/BEK/Vonkarmen_bone/scripts/investigate_boussinesq_fold.py` | `von_karman/scripts/` | 复制；输出路径改为项目内 |
| von Kármán 有限域扫描 | `compress/BEK/Vonkarmen_bone/scripts/compute_domain_branches.py` | `von_karman/scripts/` | 复制；输出路径改为项目内 |
| von Kármán 有限域数据 | `compress/BEK/Vonkarmen_bone/boussinesq_fold_analysis/` | `von_karman/data/finite_fold_analysis/` | 完整复制 |
| von Kármán 域敏感性 | `compress/BEK/Vonkarmen_bone/boussinesq_domain_branches/` | `von_karman/data/finite_domain_branches/` | 完整复制 |
| von Kármán 无穷远映射程序 | `TwoDisk/analyze_vonkarman_infinite_mapping.py` | `von_karman/scripts/` | 复制；改为导入同目录求解器 |
| von Kármán 映射汇总程序 | `TwoDisk/summarize_vonkarman_infinite_mapping.py` | `von_karman/scripts/` | 复制；改为读取项目内数据 |
| von Kármán 映射数据 | `TwoDisk/vonkarman_infinite_mapping/` | `von_karman/data/infinite_mapping/` | 完整复制，包括探索性试算 |
| von Kármán 原始笔记本 | `compress/BEK/Vonkarmen_bone/Bounssinesq.ipynb` | `von_karman/notebooks/` | 复制 |
| 双盘鞍结求解器 | `TwoDisk/two_disk_boussinesq_singularity.py` | `rotor_stator/scripts/` | 复制；参考解和输出路径改为项目内 |
| 双盘三解动力学 | `TwoDisk/analyze_three_solution_dynamics.py` | `rotor_stator/scripts/` | 复制；默认输入/输出改为项目内 |
| 双盘等温参考解 | `TwoDisk/baseflow_Res1000.npz` | `rotor_stator/reference/` | 复制 |
| 双盘原方程实现 | `TwoDisk/BaseFlow_cavity.jl`, `TwoDisk/Stability_Cavity.jl` | `rotor_stator/reference/` | 复制 |
| 双盘原始笔记本 | `TwoDisk/BaseFlow.ipynb`, `TwoDisk/stability.ipynb` | `rotor_stator/notebooks/` | 复制 |
| 双盘主结果 | `TwoDisk/boussinesq_singularity_results/` | `rotor_stator/data/boussinesq_singularity_results/` | 完整复制 |
| 双盘收敛结果 | `TwoDisk/boussinesq_singularity_convergence/` | `rotor_stator/data/boussinesq_singularity_convergence/` | 完整复制 |
| 双盘 Reynolds 数扫描 | `TwoDisk/boussinesq_singularity_re_scan/` | `rotor_stator/data/boussinesq_singularity_re_scan/` | 完整复制 |
| 双盘三解谱数据 | `TwoDisk/three_solution_dynamics/` | `rotor_stator/data/three_solution_dynamics/` | 完整复制 |
| 跨模型比较程序 | `TwoDisk/compare_saddle_node_dynamics.py` | `cross_model/scripts/` | 复制；所有输入和输出改为项目内 |
| 跨模型比较数据 | `TwoDisk/dynamical_singularity_comparison/` | `cross_model/data/dynamical_singularity_comparison/` | 完整复制 |
| 双盘报告 | `TwoDisk/BOUSSINESQ_TWO_DISK_SINGULARITY_REPORT.md` | `docs/ROTOR_STATOR_SINGULARITY_REPORT.md` | 复制 |

没有从 `Boussinesq_approximation/` 再次复制文件，因为该目录本身是早期 von Kármán 工作区的整理副本。新项目以用户指定的原始 `Vonkarmen_bone` 和当前 `TwoDisk` 为来源，避免重复和来源歧义。

