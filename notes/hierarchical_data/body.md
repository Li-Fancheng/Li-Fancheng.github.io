一张图像或一句话可以有指数多个表面实现，却由有限、可复用的局部规律生成。学习器怎样发现这些规律，把它们组织成表示，并在不同语言或模态之间保留共同关系？这个专题从层次结构数据出发，连接统计学习、生成推断和表示几何。

五篇笔记按问题展开。每篇先用直观例子建立对象，再保留推导中的中间步骤，区分精确计算、均场近似与神经网络实验。这里整理的是个人学习中的文献综述和理解，重点放在模型假设、机制之间的联系，以及结论的适用范围。

## 五篇笔记 {#contents}

| Note | 主要问题 | 关键计算 |
| --- | --- | --- |
| [01 · 层次生成与表征学习](/notes/hierarchical_data/generative_learning/) | 不同表面实现怎样被认成同一个部件？ | 小语法逐例展开、标签条件统计、梯度与样本阈值 |
| [02 · 语言统计与学习曲线](/notes/hierarchical_data/language_statistics/) | 更远的上下文何时真正有用？ | 关联平台、条件熵分解、有限样本分辨率与长尾积分 |
| [03 · 层级推断与扩散生成](/notes/hierarchical_data/inference_generation/) | 加噪以后哪些层次仍能恢复？ | 树上后验、联合采样、均场固定点与动力学关联 |
| [04 · 语法解析与隐变量学习](/notes/hierarchical_data/parsing_latents/) | 组合边界未知时怎样学习？ | 多重解析、区间动态规划、隐变量预测与可辨识性 |
| [05 · 跨语言表示与共享几何](/notes/hierarchical_data/shared_representations/) | 翻译句为何在中间层更相似？ | 共享语法、几何与坐标、邻域指标、仿射残差与概念谱 |

## 几条贯穿的联系 {#connections}

**规则与统计信号。** 生成树提供可复用的对象，但不保证它们能从有限数据辨认。第一篇建立同义关系，第二篇追踪信号怎样随距离衰减，第四篇再讨论边界和目标也未知时的困难。

**推断与学习。** 第三篇先假设字典已知，计算真正的后验；其他各篇讨论如何从数据接近这个参照。最优后验仍有歧义，与模型没有学会正确后验，是两种不同限制。

**共享内容与表面实现。** 第五篇将前面的单语言模型扩展到共享高层规则的语言。重点纳入 [A theory of platonic representations in language models](https://arxiv.org/abs/2610.07168)，并结合概念层次、统计对称性与相似性校准，区分共享几何、坐标对齐和实际计算。

所有笔记使用叶层为 $0$、根层为 $L$ 的记号；原论文的层编号不同，在比较标度公式时重新对齐。数据树的层次、网络计算层和训练时间分别定义。

## 与其他笔记的联系 {#perspective}

[神经网络理论与可解模型](/notes/neural_network_theory/)提供核方法、特征学习和深度动力学背景；其中的[层级组合学习](/notes/neural_network_theory/hierarchical_learning/)用奇偶机制模型分析弱信号和链式梯度，与这里典型随机语法的统计信号形成对照。[Scaling Laws 机制与学习曲线](/notes/scaling_laws/)进一步比较统计估计、谱结构、离散技能和动态特征的标度机制，[扩散模型的动力学](/notes/diffusion_dynamics/)则从生成过程的时间演化进入。

## 文献与研究脉络 {#references}

各篇在相关推导与结论旁链接原文。这里保留总表，便于沿同一问题比较不同模型；预印本与正式发表版本按同一项工作处理。

| 问题 | 原始文献 | 关注的联系 |
| --- | --- | --- |
| 层级生成与可学习性的早期形式化 | [Deep Learning and Hierarchal Generative Models](https://arxiv.org/abs/1612.09057)，2016 | 区分层级信息和受限低阶统计所能实现的分类 |
| 逐层学习算法 | [A Provably Correct Algorithm for Deep Learning that Actually Works](https://arxiv.org/abs/1803.09522)，2018 | 局部梯度更新、聚类与层级恢复的算法路线 |
| 随机语法与监督分类 | [How Deep Neural Networks Learn Compositional Data: The Random Hierarchy Model](https://arxiv.org/abs/2307.02129)，PRX 2024 | 随机层级模型、同义不变性和样本阈值 |
| 稀疏性与空间稳定性 | [How Deep Networks Learn Sparse and Hierarchical Data](https://arxiv.org/abs/2404.10727)，ICML 2024 | 空白位置、离散形变与权重共享 |
| 语言关联与有效上下文 | [Towards a Theory of How the Structure of Language Is Acquired by Deep Neural Networks](https://arxiv.org/abs/2406.00048)，NeurIPS 2024 | 共同祖先、关联平台与有限样本分辨率 |
| 推断与卷积架构 | [U-Nets as Belief Propagation](https://arxiv.org/abs/2404.18444)，2024 年初稿 | 层级生成模型中的推断近似、分类与去噪 |
| Transformer 的层级计算 | [How Transformers Learn Structured Data](https://proceedings.mlr.press/v267/garnier-brun25a.html)，ICML 2025 | 用层级过滤检验不同关联尺度的恢复 |
| 最优去噪中的重建转变 | [A Phase Transition in Diffusion Models Reveals the Hierarchical Nature of Data](https://arxiv.org/abs/2402.16991)，PNAS 2025 | 前向—后向实验与不同层次的重合率 |
| 用协同变化探测隐结构 | [Probing the Latent Hierarchical Structure of Data via Diffusion Models](https://arxiv.org/abs/2410.13770)，ICLR 2025 | 动力学关联、易感性和临界长度 |
| 长尾规则与学习曲线 | [Learning Curves Theory for Hierarchically Compositional Data with Power-Law Distributed Features](https://proceedings.mlr.press/v267/cagnetta25a.html)，ICML 2025 | 规则稀有性、层级前因子与任务依赖 |
| 架构与语言标度 | [Scaling Laws and Representation Learning in Simple Hierarchical Languages](https://arxiv.org/abs/2505.07070)，PRE 2025 | 局部性、共享与高层规则的统计效率 |
| 有限数据下的生成学习 | [How Compositional Generalization and Creativity Improve as Diffusion Models Are Trained](https://arxiv.org/abs/2502.12089)，ICML 2025 | 从局部合法到全局连贯的逐尺度生成 |
| 可变树与解析歧义 | [Deep Networks Learn to Parse Uniform-Depth Context-Free Languages from Local Statistics](https://arxiv.org/abs/2602.06065)，2026 年预印本 | 不固定组合边界时的统计学习 |
| 隐变量预测的效率 | [Learn from Your Own Latents and Not from Tokens](https://arxiv.org/abs/2605.27734)，2026 年预印本 | 将预测目标提升到隐层次及其适用条件 |
| 跨语言共享表示 | [A theory of platonic representations in language models](https://arxiv.org/abs/2610.07168)，2026 年 10 月预印本 | 共享高层语法、相似性与对齐、逐层仿射残差 |
| 表示收敛假说 | [Position: The Platonic Representation Hypothesis](https://proceedings.mlr.press/v235/huh24a.html)，ICML 2024 | 跨模型与跨模态表示收敛的经验假说 |
| 表示相似性的统计校准 | [Revisiting the Platonic Representation Hypothesis: An Aristotelian View](https://arxiv.org/abs/2602.14486)，ICML 2026 | 置换空模型、维度与选层效应、局部关系 |
| 概念层次与共现谱 | [Hierarchical Concept Geometry in Language Models Emerges from Word Co-occurrence](https://arxiv.org/abs/2605.23821)，2026 | 概念分类树中的由粗到细谱结构 |
| 连续统计对称性与几何 | [Symmetry in language statistics shapes the geometry of model representations](https://arxiv.org/abs/2602.15029)，2026 | 共现对称性、Fourier 模式与连续几何 |
