训练数据、模型和计算增加时，测试损失经常在一段范围内呈现稳定的幂律下降。但相似的曲线可以来自不同过程：估计波动减小、目标结构被逐渐分辨、稀有技能被覆盖，或者训练中的表示变化让学习加速。这个专题整理这些机制的模型、推导、直观例子与文献联系。

我关心的不只是指数是多少，还包括指数对应什么观测量、哪些资源被固定、何种结构正在被学会，以及换一个训练协议以后预测是否仍然成立。各篇保留前言、计算的中间步骤、跨模型比较和结语，使单篇可以连续阅读，也能沿不同机制往返。

## 六篇笔记 {#contents}

| 笔记 | 核心问题 | 主要内容 |
| --- | --- | --- |
| [01 · 统计估计与分辨率](/notes/scaling_laws/statistical_geometry/) | 更多数据是在反复测量，还是在发现新结构？ | 波动与分辨率、teacher–student 回归、随机特征投影、噪声与最优复杂度 |
| [02 · 谱与学习分辨率](/notes/scaling_laws/spectral_learning/) | 目标能量与学习速率怎样共同决定曲线？ | 统一谱风险、时间松弛、数据学习、容量截止与目标对齐 |
| [03 · 离散技能与长尾结构](/notes/scaling_laws/discrete_structure/) | 平滑损失为何能伴随局部能力转变？ | 技能尾部、随机计数、容量代价、组合共享与数据选择 |
| [04 · 特征学习与动态加速](/notes/scaling_laws/feature_learning/) | 表示变化何时改变指数？ | 自适应度量、初始核下的任务难度、自洽时钟与有限模型瓶颈 |
| [05 · 随机图与语言标度](/notes/scaling_laws/graph_language/) | 数据没有幂律，学习还能有幂律吗？ | ER 与 BA 图、Markov 预测、熵率、计数基线与简化语言 |
| [06 · 计算配置与经验外推](/notes/scaling_laws/compute_allocation/) | 各种机制怎样进入同一份预算？ | 条件最优解、数据回用、二维风险拟合、损失底与协议敏感性 |

## 从什么问题进入 {#reading-paths}

想先理解幂律的不同来源，可以从第一篇的均值估计和地图分辨率开始，再进入第二篇的连续模态与第三篇的离散技能。它们分别强调统计波动、可恢复方向和随机覆盖，也提供了彼此可以比较的简单模型。

关心神经网络与固定核的差异，可以直接阅读第四篇。它把隐藏层更新写成时变的学习度量，解释为什么某些目标只有常数改善，而另一些目标能够获得新的时间指数。第五篇进一步把数据生成过程控制起来：用简单的图游走检验哪些输入结构是必要的，哪些网络问题仍未解决。

实际的资源选择汇入第六篇。计算最优公式依赖风险曲面与成本定义；即使若干一维曲线近似幂律，也不能直接认定二维曲面可加。损失底、优化程度、参数计数和数据复用都参与最终的指数解释。

## 几个贯穿专题的区别 {#summary}

**总损失与超额损失。** 有不可约随机性或表示平台时，应说明减去的是哪个极限。交叉熵的非零损失底会改变未经校正的双对数斜率。

**数据量与训练时间。** 新的独立数据带来信息，重复读取同一数据主要改变优化过程。在线一遍训练能把步数和已见样本数联系起来，但这种映射依赖 batch 与训练协议。

**参数数、宽度与有效秩。** 一个指数换了横轴就可能改变。把网络参数换成可解析模态、技能数或每步成本，需要明确的映射。

**充分机制与必要条件。** 幂律谱能产生幂律学习曲线，长尾技能也能；它们各自的可解推导不能证明所有网络都必须沿同一机制工作。ER 随机游走提供了很直接的对照。

**预测与解释。** 一个回归器可以准确拟合预算内的风险面，却未识别微观来源；一个机制模型可以解释某种效应，却未覆盖全部实验区间。二者可以互相约束。

## 文献入口 {#literature}

| 机制与问题 | 代表性原始文献 |
| --- | --- |
| 统计学习与四种区域 | [Explaining Neural Scaling Laws](https://arxiv.org/abs/2102.06701)，初稿 2021，PNAS 2024 |
| 核谱与目标对齐 | [Spectrum Dependent Learning Curves in Kernel Regression and Wide Neural Networks](https://proceedings.mlr.press/v119/bordelon20a.html)，ICML 2020 |
| 联合谱与随机特征模型 | [A Solvable Model of Neural Scaling Laws](https://arxiv.org/abs/2210.16859)，2022 |
| 离散知识与技能 | [The Quantization Model of Neural Scaling](https://arxiv.org/abs/2303.13506)，NeurIPS 2023 |
| 数据选择 | [Beyond neural scaling laws: beating power law scaling via data pruning](https://arxiv.org/abs/2206.14486)，NeurIPS 2022 |
| 动态特征与指数改善 | [How Feature Learning Can Improve Neural Scaling Laws](https://arxiv.org/abs/2409.17858)，ICLR 2025 |
| 受控随机图与语言 | [On the origin of neural scaling laws: from random graphs to natural language](https://arxiv.org/abs/2601.10684)，ICML 2026 |
| 经验资源配置 | [Scaling Laws for Neural Language Models](https://arxiv.org/abs/2001.08361)，2020；[Training Compute-Optimal Large Language Models](https://arxiv.org/abs/2203.15556)，2022 |
| 数据回用与协议比较 | [Scaling Data-Constrained Language Models](https://arxiv.org/abs/2305.16264)，2023；[Resolving Discrepancies in Compute-Optimal Scaling of Language Models](https://arxiv.org/abs/2406.19146)，2024 |

[神经网络理论与可解模型](/notes/neural_network_theory/)提供高维泛化、核方法与训练动力学背景；[层次结构数据与学习理论](/notes/hierarchical_data/)进一步讨论组合生成、长程语言关联和共享表示。三个专题可以沿“数据结构怎样变成可学习的表示”这条线连接起来。

专题中的数值图包括可解模型的直接计算与有限随机图统计，各图注明观测量、假设和计算范围。原论文的网络训练结果通过原始文献讨论，机制示意与独立计算不会被当作训练复现。

[文献来源表](sources.json) · [数值记录与补充材料](reproducibility.zip)
