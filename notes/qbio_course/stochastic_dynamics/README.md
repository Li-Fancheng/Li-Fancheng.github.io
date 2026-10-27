# 生命系统中的随机动力学

QBIO 学习笔记，阅读顺序见 [主题目录](stochastic_dynamics_note.md)。三篇依次讨论随机性与隐藏状态、从事件建立动力学、时间关联与罕见事件的后果。网站入口为 `/notes/qbio_course/stochastic_dynamics/`，目录与篇间导航由仓库的 `data/qbio_series.json` 统一生成。

六个主要推导是随机游走到扩散、隐藏构象的等待时间与危险率、生灭主方程及泊松稳态、跳跃过程到 Langevin 与 Fokker–Planck、相关输入通过线性下游后的方差、首达时间与弱噪声逃逸。主要计算随物理问题在正文展开，补充计算在篇末折叠。引用的原始论文在各篇末尾；7 张图片均为模型计算或示意，没有使用实验数据或转载论文图片。

## 复现图与核对计算

在本目录安装 `requirements.txt` 后运行：

```bash
python generate_figures.py
```

脚本可从任意工作目录运行。它使用 QBIO 共用的 `../plot_style.py`，按 `data/qbio_figures.json` 中的展示宽度导出 SVG 与 300 dpi PNG，随机种子固定。SVG 用于网站，PNG 便于独立查看。

| 图号 | 文件名 | 内容 |
|---|---|---|
| 1 | `01_walk_and_diffusion` | 二项分布的扩散极限与均方位移 |
| 2 | `02_hidden_state_memory` | 隐藏构象的混合等待时间与危险率 |
| 3 | `03_birth_death_dynamics` | Gillespie 轨迹、平均恢复与瞬态分布 |
| 4 | `04_same_distribution_different_clocks` | 相同泊松稳态下的关联和谱 |
| 5 | `05_noise_filtering` | 快慢 OU 输入通过同一个线性下游 |
| 6 | `06_barrier_and_first_passage` | 双阱势、首达积分与 Kramers 近似 |
| 7 | `07_same_occupancy_different_currents` | 相同占据和驻留时间下的不同环流 |

脚本还输出 `numerical_checks.json`：核对泊松分布的概率收支；比较频谱积分与下游矩方程；独立求解后向边值问题，和首达时间双积分比较；比较离散灭绝时间的递推与线性方程组；检查三态环的稳态、驻留时间与概率流。

图 5 通过精确 OU 更新生成输入，采用左端点分段恒定输入推进下游，步长为 `0.005`，略去下游最初 8 个时间单位。它用于显示轨迹差别；方差曲线和正文数值来自解析计算。图 6 的边界为左侧 `-2` 反射、右侧 `0.8` 吸收，起点为 `-1`，不把“到达势垒顶”与“越过势垒”混用。
