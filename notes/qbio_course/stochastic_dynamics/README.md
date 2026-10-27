# 生命系统中的随机动力学

QBIO 学习笔记，阅读顺序见[主题目录](stochastic_dynamics_note.md)。网站入口为 `/notes/qbio_course/stochastic_dynamics/`，标题、阅读次序与篇间导航由 `data/qbio_series.json` 管理。

三篇从单细胞表达差异出发，依次讨论随机性与隐藏状态、从事件建立动力学、时间关联与稀有事件的后果。用同一个生灭模型连接主方程、随机模拟、Langevin、Fokker–Planck 与 OU，再加入启动子切换、首达边界和最小循环。主要推导随问题展开，技术计算归入连续展示的附录 A、B、C，各篇末尾提供原始论文链接。

正文与附录共使用 10 张模型图，没有使用实验数据或转载论文图片。生成脚本另保留两张独立高斯／相关序列和二维旋转 OU 的对照图，作为源码中的额外算例，不进入三篇阅读主线。

## 复现图与核对计算

安装本目录 `requirements.txt` 后，在本目录运行：

```bash
python generate_figures.py
```

脚本也可通过其文件路径从其他目录调用。它使用 QBIO 共用的 `../plot_style.py`，按 `data/qbio_figures.json` 中的展示宽度导出 SVG 与 300 dpi PNG，随机种子固定。SVG 用于网站，PNG 便于独立查看。

| 图号 | 文件名 | 内容 |
|---|---|---|
| 1.1 | `01_walk_and_diffusion` | 二项分布到扩散极限 |
| 1.2 | `02_hidden_state_memory` | 隐藏构象、混合等待与危险率 |
| 1.3 | `09_persistence_and_diffusion` | 附录中的持续运动与长时扩散 |
| 2.1 | `03_birth_death_dynamics` | 生灭轨迹、瞬态分布和均值 |
| 2.2 | `10_discrete_and_gaussian` | 连续近似、低计数和边界 |
| 2.3 | `04_same_distribution_different_clocks` | 相同稳态下的不同关联和谱 |
| 3.1 | `05_noise_filtering` | 表达输入经线性下游的传递 |
| 3.2 | `06_barrier_and_first_passage` | 首达积分与 Kramers 近似 |
| 3.3 | `11_quasistationary_extinction` | 生存概率衰减与准稳态形状 |
| 3.4 | `07_same_occupancy_different_currents` | 相同占据下的不同环流 |

脚本输出 `numerical_checks.json`。计算核对包括泊松稳态的概率流、下游方差的频谱积分与矩方程、首达双积分与独立后向边值解、灭绝递推与线性方程组、持续运动的协方差积分、双态启动子的联合主方程与稳态矩，以及准稳态寿命。

图 3.1 的输入采用精确 OU 更新，下游采用左端点分段恒定输入推进，步长 `0.005`，略去最初 8 个时间单位；方差曲线来自解析结果。图 3.2 的左边界 `-2` 反射、目标 `0.8` 吸收、起点 `-1`，对应越过势垒后的首达事件。图 3.3 用完整有限状态主方程演化，条件分布除以尚未吸收的总概率。
