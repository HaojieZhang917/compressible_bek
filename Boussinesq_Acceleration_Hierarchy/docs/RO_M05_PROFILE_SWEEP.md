# Ro=-0.5 加热基本流剖面扫描

## 目的

在传统 Boussinesq BEK 相似性模型内，计算 `Ro=-0.5`、
`Tw=1.0:0.1:1.6` 的基本流，观察接近热尾失效上限时速度场与温度场的连续变化。
本扫描属于完全自洽基本流计算，不是冻结基本流或温度项开关诊断。

## 模型与数值参数

- 方程和热强迫参数化：`work/src/BEKTraditionalForcing.jl`；
- `Pr=0.72`；
- Pier 型无穷映射：`N=120, a=2, b=0.6, c=0.5`；
- 正则延拓参数：`B=Lambda_cf(Ro)*(Tw-1)`；
- Newton 残差容限：`1e-10`；
- 壁温：`Tw=[1.0,1.1,1.2,1.3,1.4,1.5,1.6]`；
- 每个壁温点使用前一温度点作为初值，按固定 `B` 求解。

## 输出

所有新结果写入
`work/results/ro_m05_profile_evolution/N120_a2.0_b0.6_c0.5/`，不修改 `baselines/`。

- `summary.csv`：`Tw, B, Hinf, ell_T, residual, iterations`；
- `profiles.csv`：全部映射节点上的 `H,F,G,Theta,T`；
- `ro_m05_profile_evolution.svg`：速度、温度及热尾演化图。

## 解释边界

这里的“失效”指传统相似性基本流的热尾极限。它不等价于完整三维稳定性失效，
也不能单独证明有限温差下物理 Boussinesq 近似的定量精度。
