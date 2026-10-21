Grokking 指训练拟合之后明显滞后的泛化改善。它把几件通常同时发生的事情分开：训练样本被拟合、预测间隔扩大、共享表示形成，以及原有函数成分之间的相对权重改变。这个专题围绕这些时间尺度，整理可解模型、表示机制、实验记录与文献联系。

我更关心怎样区分解释，而不是给所有延迟寻找一个共同的故事。固定特征的慢方向可以产生延迟，正则化可以改变已有解的选择，特征学习可以改变网络实际可用的方向；同一条准确率曲线也可能同时包含这些因素。各篇从直观问题进入，保留中间推导和条件，再把结论接回实际观测。

## 五篇笔记 {#contents}

| 笔记 | 核心问题 | 主要内容 |
| --- | --- | --- |
| [01 · 现象与可解基线](/notes/grokking/phenomenology/) | 很晚的准确率改善究竟说明什么？ | 延迟定义、阈值评分、固定模态、随机设计谱与噪声 |
| [02 · 正则化与解的竞争](/notes/grokking/solution_competition/) | 为什么快找到的解会被慢解取代？ | 精确竞争模型、乘法平台、电路效率、范数几何与优化器 |
| [03 · 表示形成与算法电路](/notes/grokking/representation_circuits/) | 网络如何形成可复用的计算？ | Lazy-to-rich、二次教师、Fourier 算法、对称投影与消融 |
| [04 · 有效理论与学习相图](/notes/grokking/effective_theory/) | 数据何时足以约束一种共享表示？ | 平行四边形关系、零模、尺度保持、谱隙与有限时间相图 |
| [05 · 实验与机制诊断](/notes/grokking/experiments_diagnosis/) | 哪些观测能区分这些机制？ | 真实模加法记录、失败对照、末态函数干预与预训练类比 |

## 几条可以连接起来的主线 {#connections}

**慢时间与评分。** 第一篇算清固定表示下的延迟，说明训练与测试的不同权重、容差和有限谱隙怎样改变表观时间。它提供一个参照，使表示解释能够说明自己额外解释了什么。

**长期偏好与可达路径。** 第二篇将快拟合、慢选择和有效表示代价分开。弱正则化可以在很长时间后产生有限影响，深度参数化则同时改变启动速度和偏好的解；范数变化本身并不足以辨认这两种作用。

**形成结构与使用结构。** 第三、四篇从两个方向进入表示：一个解析网络怎样实现周期算法，另一个解析数据怎样约束几何。存在紧凑解、训练形成了它、预测依赖它，以及它何时改变泛化，需要不同证据。

**从终点到时间因果。** 第五篇保留同一组训练权重与曲线，并补充可重算的函数干预。原模型、对称投影和单频函数可以具有同样的准确率，却有很不同的交叉熵。这些对照帮助区分性能、结构和机制。

## 原始文献入口 {#literature}

| 问题 | 代表性文献 |
| --- | --- |
| 原始延迟泛化现象 | [Grokking: Generalization Beyond Overfitting on Small Algorithmic Datasets](https://arxiv.org/abs/2201.02177)，2022 |
| 线性可解解释 | [Grokking in Linear Estimators](https://arxiv.org/abs/2310.16441)，2023 |
| 表示约束与相图 | [Towards Understanding Grokking](https://arxiv.org/abs/2205.10343)，2022 |
| 范数几何与更一般数据 | [Omnigrok: Grokking Beyond Algorithmic Data](https://arxiv.org/abs/2210.01117)，2022，ICLR 2023 |
| 简化网络的解析模运算解 | [Grokking modular arithmetic](https://arxiv.org/abs/2301.02679)，2023 |
| 非神经模型的特征迭代 | [Emergence in non-neural models](https://arxiv.org/abs/2407.20199)，2024 |
| 电路及其形成过程 | [Progress measures for grokking via mechanistic interpretability](https://arxiv.org/abs/2301.05217)，ICLR 2023 |
| 规则与记忆的效率竞争 | [Explaining grokking through circuit efficiency](https://arxiv.org/abs/2309.02390)，2023，ICLR 2024 |
| 晚期特征学习 | [Grokking as the Transition from Lazy to Rich Training Dynamics](https://arxiv.org/abs/2310.06110)，2023，ICLR 2024 |
| 早晚期偏置的理论分离 | [Dichotomy of Early and Late Phase Implicit Biases Can Provably Induce Grokking](https://arxiv.org/abs/2311.18817)，2023，ICLR 2024 |
| 自适应优化器的晚期行为 | [The Slingshot Mechanism](https://arxiv.org/abs/2206.04817)，2022 |
| 超出欧氏范数的结构偏好 | [Grokking Beyond the Euclidean Norm of Model Parameters](https://arxiv.org/abs/2506.05718)，ICML 2025 |
| 嵌入更新与分层动力学 | [Mechanistic Insights into Grokking from the Embedding Layer](https://arxiv.org/abs/2505.15624)，2025 |
| 语言预训练中的延迟类比 | [A Pre-Training Analogue of Grokking in Language Models](https://arxiv.org/abs/2606.00230)，2026 年预印本 |

## 与其他专题的联系 {#related-topics}

[神经网络理论与可解模型](/notes/neural_network_theory/)提供线性泛化、核方法、深度参数化和随机优化的背景。[Scaling Laws 机制与学习曲线](/notes/scaling_laws/)讨论慢谱、特征学习与资源效率，[层次结构数据与学习理论](/notes/hierarchical_data/)讨论可复用规则、组合结构与长程语言统计。Grokking 将这些问题放到训练拟合以后的时间轴上比较。

这里的数值内容包括可解模型计算、有限约束统计、已有模加法训练和末态函数分析，各图分别说明数据与假设。真实训练和显式算法构造承担不同结论；终点消融不替代尚未记录的中间表示历史。

[文献来源表](sources.json) · [数值记录与补充材料](reproducibility.zip)
