复杂神经网络扩大数据、模型和训练计算时，为什么经常表现出稳定的幂律？这些曲线可能反映有限样本的统计波动、待学习结构的分布，也可能反映训练中不断变化的有效学习尺度。多种过程都能产生幂律，曲线形状本身不足以识别机制。

从经验发现进入具体模型，再综合它们解释的误差、所取的资源极限和相互联系，可以看清已有理论怎样组织起来，以及各项结论的适用范围。

## 1. 八篇笔记与各自的任务 {#contents}

| 笔记 | 理论角色 | 核心内容 |
| --- | --- | --- |
| [01 · 从学习曲线到 Neural Scaling Laws](/notes/scaling_laws/empirical_origins/) | 经验入口 | 经典学习曲线，Hestness 的组合曲线，Rosenfeld 的二维误差面，Kaplan／Chinchilla 的计算前沿 |
| [02 · 统计估计与几何分辨率](/notes/scaling_laws/statistical_geometry/) | 静态机制 | 统计波动、局部插值、teacher–student 回归、投影平台与噪声复杂度 |
| [03 · 谱与学习分辨率](/notes/scaling_laws/spectral_learning/) | 静态机制 | 目标加权速率谱，时间松弛、随机采样、容量截止与目标对齐 |
| [04 · 离散技能与长尾结构](/notes/scaling_laws/discrete_structure/) | 静态机制及其联系 | 技能代价、有限覆盖、组合共享，以及谱松弛与随机覆盖的指数混合结构 |
| [05 · 特征学习与动态加速](/notes/scaling_laws/feature_learning/) | 动态机制 | 时变学习度量、任务难度、反馈时钟与有限模型瓶颈 |
| [06 · 随机图与语言标度](/notes/scaling_laws/graph_language/) | 受控机制检验 | 已知生成器、熵率和计数基线；检验输入幂律是否为必要条件 |
| [07 · 计算配置与经验外推](/notes/scaling_laws/compute_allocation/) | 预测与配置 | 已知曲面上的优化、未知曲面的辨认、梯度噪声与 batch 效率 |
| [08 · 理论脉络与综合认识](/notes/scaling_laws/mechanism_synthesis/) | 理论综述 | 不同解释层次、静态与动态表示的联系、已有认识与适用范围 |

## 2. 快速阅读与完整推导 {#reading-paths}

**机制概览。** 先读约 29 分钟的 [01 篇](/notes/scaling_laws/empirical_origins/)，了解经验上究竟需要解释什么；再按下面的结果导航选择问题，最后读 [08 篇的理论综述](/notes/scaling_laws/mechanism_synthesis/#comparison)。各篇开头的学习问题和结尾的边界讨论可以串成一条较短路径。

**完整理论推导。** 从 02–04 的统计、谱与技能出发，进入 05 的动态表示；用 06 的受控任务检验必要条件，再到 07 的风险面与资源配置。01 和 08 分别提供问题入口与综合。

| 想先弄清的问题 | 核心结果或比较入口 |
| --- | --- |
| 更多数据是在测准旧参数，还是分辨新结构？ | [四种受限区域](/notes/scaling_laws/statistical_geometry/#four-regimes)与[噪声下的最优分辨率](/notes/scaling_laws/statistical_geometry/#bias-variance) |
| 单个模态按指数恢复，总体为何是幂律？ | [小速率端与时间指数](/notes/scaling_laws/spectral_learning/#小速率端决定长期指数) |
| 谱学习与技能覆盖有什么共同数学结构？ | [累计测度、Gamma 常数与物理区别](/notes/scaling_laws/discrete_structure/#mixture-bridge) |
| 表示学习何时只改常数，何时改变指数？ | [自洽学习时钟](/notes/scaling_laws/feature_learning/#self-consistent-clock) |
| 输入没有预置幂律，网络还能有经验幂律吗？ | [ER 实验排除的必要条件](/notes/scaling_laws/graph_language/#discussion-necessity) |
| 预测风险面与求预算最优，差别在哪里？ | [风险面拟合](/notes/scaling_laws/compute_allocation/#fitting-surface)与[边际最优条件](/notes/scaling_laws/compute_allocation/#最优条件的一般边际解释) |
| 已有理论怎样联系，已经解释了哪些问题？ | [各理论的解释任务与相互联系](/notes/scaling_laws/mechanism_synthesis/#comparison) |

## 3. 几个贯穿专题的区别 {#summary}

**总损失与超额损失。** 有不可约随机性或表示平台时，应说明减去的是哪个极限。交叉熵的非零损失底会改变未经校正的双对数斜率。

**数据量与训练时间。** 新的独立数据带来信息，重复读取同一数据主要改变优化过程。在线一遍训练能把步数和已见样本数联系起来，但这种映射依赖 batch 与训练协议。

**参数数、宽度与有效秩。** 换横轴可能改变指数。把网络参数换成可解析模态、技能数或每步成本，需要明确映射。

**充分机制与起源解释。** 给定幂律谱或长尾技能，能够推出学习曲线；真实网络为什么具有这样的加权结构，还需要独立解释。输入、初始核、训练后权重与技能频率也不是同一张谱。

**预测与机制辨识。** 回归器可以准确预测预算内的风险面，却未识别微观来源；机制模型可以解释某种效应，却未覆盖全部实验区间。留出尺度与受控干预分别约束这两类结论。

## 4. 文献入口 {#literature}

| 机制与问题 | 代表性原始文献 |
| --- | --- |
| 经典学习曲线 | [Statistical mechanics of learning from examples](https://doi.org/10.1103/PhysRevA.45.6056)，1992；[A Universal Theorem on Learning Curves](https://doi.org/10.1016/0893-6080%2893%2990013-M)，1993 |
| 现代经验入口 | [Deep Learning Scaling is Predictable, Empirically](https://arxiv.org/abs/1712.00409)，2017 |
| 二维误差面 | [A Constructive Prediction of the Generalization Error Across Scales](https://arxiv.org/abs/1909.12673)，2019 / ICLR 2020 |
| 统计学习与四种区域 | [Explaining Neural Scaling Laws](https://arxiv.org/abs/2102.06701)，PNAS 2024 |
| 核谱与目标对齐 | [Spectrum Dependent Learning Curves in Kernel Regression and Wide Neural Networks](https://proceedings.mlr.press/v119/bordelon20a.html)，ICML 2020 |
| 联合谱与随机特征模型 | [A Solvable Model of Neural Scaling Laws](https://arxiv.org/abs/2210.16859)，2022 |
| 离散技能与共享 | [The Quantization Model of Neural Scaling](https://arxiv.org/abs/2303.13506)，NeurIPS 2023 |
| 动态特征与指数改善 | [How Feature Learning Can Improve Neural Scaling Laws](https://arxiv.org/abs/2409.17858)，ICLR 2025 |
| 受控随机图与语言 | [On the origin of neural scaling laws: from random graphs to natural language](https://arxiv.org/abs/2601.10684)，2026 |
| 经验资源配置 | [Scaling Laws for Neural Language Models](https://arxiv.org/abs/2001.08361)，2020；[Training Compute-Optimal Large Language Models](https://arxiv.org/abs/2203.15556)，2022 |
| 梯度噪声与并行效率 | [An Empirical Model of Large-Batch Training](https://arxiv.org/abs/1812.06162)，2018 |
| 配置复现与拟合可靠性 | [Resolving Discrepancies in Compute-Optimal Scaling of Language Models](https://arxiv.org/abs/2406.19146)，2024；[(Mis)Fitting: A Survey of Scaling Laws](https://arxiv.org/abs/2502.18969)，2025 |

[神经网络理论与可解模型](/notes/neural_network_theory/)提供高维泛化、核方法与训练动力学背景；[层次结构数据与学习理论](/notes/hierarchical_data/)进一步讨论组合生成、长程语言关联和共享表示。

## 5. 文献与补充材料 {#resources}

补充材料包含八篇正文、十三张计算图、文献来源总表、计算记录与公式核对结果。数值图来自可解模型与有限图统计，独立算例的假设和范围随图说明。

[文献来源总表](sources.json) · [数值记录与补充材料](reproducibility.zip)
