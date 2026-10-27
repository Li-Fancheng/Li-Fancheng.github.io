# 随机过程与随机动力学

QBIO 中以理论联系组织的三篇学习笔记，阅读顺序见[主题目录](stochastic_dynamics_note.md)。网站入口为 `/notes/qbio_course/stochastic_dynamics/`，标题、阅读次序与篇间导航由 `data/qbio_series.json` 管理。

三篇依次讨论状态与粗粒化、生成元与连续极限、涨落与首达及非平衡。生物模型用来检验状态选择、事件结构和边界条件。主要推导随正文展开，技术计算归入连续展示的附录 A、B、C；各篇末尾提供原始论文链接。12 张图片均为模型计算或示意，没有使用实验数据或转载论文图片。

## 复现图与核对计算

安装本目录 `requirements.txt` 后，在本目录运行：

```bash
python generate_figures.py
```

脚本使用 QBIO 共用的 `../plot_style.py`，按 `data/qbio_figures.json` 中的展示宽度导出 SVG 与 300 dpi PNG，随机种子固定。SVG 用于网站，PNG 便于独立查看。图号按各篇阅读次序编号，文件名维持稳定。

| 图号 | 文件名 | 内容 |
|---|---|---|
| 1.1 | `08_marginals_and_memory` | 同一高斯边缘分布，不同路径和相关 |
| 1.2 | `09_persistence_and_diffusion` | 持续运动的短时弹道与长时扩散 |
| 1.3 | `01_walk_and_diffusion` | 二项分布到扩散极限 |
| 1.4 | `02_hidden_state_memory` | 隐藏构象、混合等待与危险率 |
| 2.1 | `03_birth_death_dynamics` | 路径、瞬态分布和均值 |
| 2.2 | `10_discrete_and_gaussian` | 连续近似、低计数和边界 |
| 2.3 | `04_same_distribution_different_clocks` | 同一稳态下的不同关联和谱 |
| 3.1 | `05_noise_filtering` | 有色输入经线性响应核的传递 |
| 3.2 | `06_barrier_and_first_passage` | 首达积分与 Kramers 渐近式 |
| 3.3 | `11_quasistationary_extinction` | 生存概率衰减与准稳态形状 |
| 3.4 | `07_same_occupancy_different_currents` | 同一占据和驻留时间下的不同环流 |
| 3.5 | `12_gaussian_density_and_circulation` | 同一二维高斯密度下的零流与旋转 |

脚本输出 `numerical_checks.json`。计算核对包括泊松稳态的概率流、下游方差的频谱积分与矩方程、首达双积分与独立后向边值解、灭绝递推与线性方程组、持续运动的协方差积分、批量生成的精确稳态与矩、准稳态特征值与平均寿命、多变量 OU 的 Lyapunov 方程与时间关联。

图 3.1 的输入采用精确 OU 更新，下游采用左端点分段恒定输入推进，步长 `0.005`，略去最初 8 个时间单位；方差曲线来自解析结果。图 3.2 的左边界 `-2` 反射、目标 `0.8` 吸收、起点 `-1`，对应越过势垒后的首达事件。图 3.3 用完整有限状态矩阵指数演化，条件分布除以尚未吸收的总概率。图 3.5 的两个状态坐标在时间反演下均按偶变量处理。
