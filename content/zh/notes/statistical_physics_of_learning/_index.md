---
title: "神经网络的统计物理：好工作的近似方法"
date: 2026-10-30
summary: "以近似和可解极限为线索，阅读 scaling laws、扩散模型、层次数据与特征学习中的统计物理工作。"
layout: series
show_article_cards: false
article_class: neural-theory-note
math: true
tags:
  - 统计物理
  - 神经网络理论
  - 文献阅读
---

为什么高维学习问题可以被少数变量描述？一篇漂亮的理论工作常常先替换问题：固定特征、构造随机数据、取高维极限，或把学习过程投影到少数序参量上。这里收集这样的工作，重点不是罗列结论，而是辨认**哪一步让问题可解、保留了什么、又丢掉了什么**。文献不限期刊和会议；这是一份可继续扩展的阅读地图。

已有笔记可配合阅读：[神经网络理论与可解模型](/notes/neural_network_theory/)、[Scaling Laws](/notes/scaling_laws/)、[层次结构数据与 RHM](/notes/hierarchical_data/) 和 [扩散模型的动力学](/notes/diffusion_dynamics/)。

## 一张方法地图

| 想解释的问题 | 主要替换或极限 | 可计算的对象 | 首先要问的代价 |
| --- | --- | --- | --- |
| 数据与模型大小的标度 | 随机特征、核谱、平滑流形 | 模态误差、学习曲线 | 特征是否固定？谱和目标如何假设？ |
| 生成、泛化与记忆 | 最优经验 score、高维分布 | 谱不稳定性、熵与凝聚阈值 | 训练误差与正则化是否被忽略？ |
| 数据组合结构如何被学习 | 随机层次生成模型 | 局部相关、样本复杂度、层级表示 | 合成数据的生成规则是否决定了答案？ |
| 特征学习与相变 | teacher–student、贝叶斯推断、比例极限 | 重叠序参量、鞍点、泛化误差 | 是最优推断还是实际训练动力学？ |
| 训练中的集体动力学 | 动态平均场、响应函数与相关函数 | 自洽动力学、学习曲线 | 大宽度如何取极限？记忆与噪声是否保留？ |

## 1. Scaling laws：四个区间背后的机制

**Bahri, Dyer, Kaplan, Lee & Sharma, [Explaining neural scaling laws](https://doi.org/10.1073/pnas.2311878121)**，*PNAS* **121**, e2311878121 (2024)。

- **问题与模型：** 测试误差随数据量 $D$ 和模型大小 $P$ 的幂律从哪里来？作者先分辨数据量和模型大小各自的 **variance-limited** 与 **resolution-limited** 区间，共四种情形。
- **近似的关键：** 在随机特征 teacher–student 模型中，将非线性学习问题化为线性特征上的回归，并用核算子的特征值与目标在各模态上的分量描述误差。另一条论证把数据视为平滑低维流形，估计有限资源能解析的空间尺度。方差受限区间则考察趋近无限数据或无限宽度极限的有限规模修正。
- **得到什么：** 有限抽样或有限特征造成的涨落与尚未解析的细尺度结构给出不同标度机制；在随机特征设定中，模型大小与数据大小的 resolution-limited 指数存在联系。文中也用深网络和真实数据检验部分预测。
- **边界：** 随机特征不随任务学习；流形平滑性和谱衰减是解释指数的结构输入，不是所有网络自动满足的结论。深网络的参数数目与宽度的换算也依赖架构。**值得细读：** 将误差写成已解析和未解析谱模态的贡献，检查有限 $D$、$P$ 的极限顺序。

**Bordelon, Atanasov & Pehlevan, [A Dynamical Model of Neural Scaling Laws](https://proceedings.mlr.press/v235/bordelon24a.html)**，*ICML* (2024)。

- **近似：** 用梯度下降训练的随机特征模型代替完整网络，并在高维极限下跟踪训练时间、数据和特征数共同决定的学习动力学；数据重复使用带来的时间相关性不能简单当作每一步独立采样。
- **收获与边界：** 模型解释训练时间和模型大小指数为何可能不同，并给出计算最优分配的示例；它是固定特征模型，不能直接说明深网络的特征演化。可与 PNAS 文比较：一个强调收敛后四个区域，另一个把**时间**带回来。

## 2. 扩散生成：从经验 score 到玻璃凝聚

**Biroli, Bonnaire, De Bortoli & Mézard, [Dynamical regimes of diffusion models](https://www.nature.com/articles/s41467-024-54281-3)**，*Nature Communications* **15**, 9957 (2024)。

- **近似的第一步：** 假设模型学到训练集的**精确经验 score**，再同时考虑高维和大量样本。这样反向扩散可以由加噪经验分布来分析，而无需解神经网络训练本身。高斯混合提供可解析例子；更一般的数据再通过谱和熵分析并与数值实验比较。
- **物理图像：** 反向过程中依次出现噪声、类别分化和向具体训练点坍缩三个区域。前者的分化与协方差最大特征值有关；后者可用训练点权重的凝聚与随机能量模型类比来理解。论文时间 $t$ 从**正向**加噪起点计，因此生成过程沿时间反向走。
- **一个可检查的式子：** 在论文的约定下，分化时间满足 $\Lambda e^{-2t_S}=1$，其中 $\Lambda$ 是数据协方差的最大特征值。它展示了如何把高维动力学问题降到一个谱不稳定性判据。
- **边界：** 精确经验 score 最终回到经验分布，因此末期记忆是该设定内的结果；实际有限网络、提前停止和正则化可能偏离这一极限。玻璃凝聚类比不意味着所有实际扩散模型都会按同样方式坍缩。**值得细读：** 分离目标训练点与其他点的权重之和，为什么指数权重使普通中心极限定理不够用？

## 3. 层次数据：把复杂世界换成可控生成规则

**Cagnetta, Petrini, Tomasini, Favero & Wyart, [How Deep Neural Networks Learn Compositional Data: The Random Hierarchy Model](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.14.031001)**，*Physical Review X* **14**, 031001 (2024)。

- **模型替换：** 用随机层次模型（RHM）生成分类数据：高层特征通过规则逐层展开成低层特征，同一高层特征可有多个等价展开。它保留组合结构和局部不变性，却去掉真实语言、图像中的大量细节。
- **近似与观察：** 通过随机规则的统计性质估计标签与局部特征的相关及其有限样本涨落，再将可检测相关与学习所需样本数联系起来；训练网络并观察逐层表示及不变性。其力量在于把**数据结构**显式放入问题，而不仅仅取随机标签或各向同性高斯输入。
- **边界：** 样本复杂度估计依赖 RHM 的规则和所讨论的网络与学习机制；相关可检测不等于任意优化算法必然发现它。RHM 是结构化数据模型，并非自动适用于所有真实数据。站内已有 [RHM 与层级学习的系列推导](/notes/hierarchical_data/)。

## 4. 从固定特征到可学习表示

**Barbier 等, [Statistical Physics of Deep Learning: Optimal Learning of a Multilayer Perceptron near Interpolation](https://journals.aps.org/prx/abstract/10.1103/56sb-pdh6)**，*Physical Review X* **16**, 031014 (2026)。

- **设定与方法：** teacher–student 的贝叶斯最优学习；宽度与输入维度同比增长，样本数与参数数目处于可比较的插值尺度。作者结合 replica 与随机矩阵的 HCIZ 方法，保留这一区域中重要的矩阵结构及相关性。
- **物理内容：** 在所研究模型中，随着样本增加，先出现利用特征组合的阶段，之后出现单元对教师特征的特化；深层结构中也可讨论特化的逐层传播。
- **边界：** 这是指定生成模型和先验下的 Bayes-optimal 泛化，而非对实际 SGD 训练速度的直接计算；replica 鞍点及其稳定性、输入分布和激活函数都是阅读时需要检查的条件。

**Lauditi, Bordelon & Pehlevan, [Adaptive kernel predictors from feature-learning infinite limits of neural networks](https://proceedings.mlr.press/v267/lauditi25a.html)**，*ICML* (2025)。

- **方法：** 在特征学习的无限宽度极限研究带权重衰减的梯度流，以动态平均场理论（DMFT）的自洽方程描述表征和核；固定点定义任务适应的核预测器。
- **为什么值得比较：** 同样取无限宽度，不同参数化和极限可以得到固定的 lazy kernel，也可以得到适应任务的表示。极限宽度本身不等于没有特征学习。其固定点预测仍依赖指定动力学与极限，并不自动刻画有限网络所有训练路径。

**[Predictive Power of a Bayesian Effective Action for Fully Connected One Hidden Layer Neural Networks in the Proportional Limit](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.133.027301)**，*Physical Review Letters* **133**, 027301 (2024)。

- **方法：** 对单隐藏层贝叶斯网络采用样本数与宽度同比增长的有效作用量，并与 Langevin 训练的数值实验比较。该理论中无限宽核的整体重标度是重要机制。
- **边界：** 有效理论、先验、温度和采样动力学均属于问题定义；与实验吻合应理解为所测试范围内的支持，不能替代对其他架构和训练方式的验证。

**[Statistical Mechanics of Transfer Learning in Fully Connected Networks in the Proportional Limit](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.134.177301)**，*Physical Review Letters* **134**, 177301 (2025)。

- **方法：** 比例极限下用单样本 Franz–Parisi 形式描述源任务与目标任务之间的关联，得到重整化的跨任务核，判断迁移何时有益。
- **边界：** 在这一设定下有效的迁移依赖任务关系与贝叶斯模型；它给出可计算的机制，并不意味着所有迁移学习都可归结为同一个核。

## 5. 动力学闭合与现象论

**Kramp, Lindner & Helias, [Dynamics of neural scaling laws in random feature regression with powerlaw-distributed kernel eigenvalues](https://proceedings.mlr.press/v306/kramp26a.html)**，*ICML* (2026)。

- **近似与方法：** 对幂律核谱的随机特征回归建立 DMFT，用集体响应函数和有效噪声的自洽统计连接各谱模态，讨论梯度流、权重衰减和 Langevin 动力学。
- **边界：** 谱及目标结构是输入；谱模态在有效描述中仍通过集体量耦合。它展示怎样把许多自由度约化为少数时间相关函数，不是完整深网络特征学习的通解。

**Shi, Pan & Dokmanić, [Spring-Block Theory of Feature Learning in Deep Neural Networks](https://journals.aps.org/prl/abstract/10.1103/ys4n-2tj3)**，*Physical Review Letters* **134**, 257301 (2025)。

- **方法：** 从层间数据分离的宏观现象出发，构造具有非对称摩擦的随机弹簧—滑块链，组织噪声和非线性下的特征学习相图。
- **边界：** 这是自上而下的**现象模型**；它重现相关现象，但不是从每个神经元的微观训练方程严格推出的等价映射。与 DMFT 或 replica 对照，恰好能看到“拟合宏观机制”和“从微观模型推导”之间的不同。

**Biehl 等, [Dynamics and Generalization Ability of LVQ Algorithms](https://www.jmlr.org/beta/papers/v8/biehl07a.html)**，*Journal of Machine Learning Research* **8**, 323–360 (2007)。

- **方法：** 在线学习中用两个竞争原型处理高斯混合数据，在高维极限描述原型与数据结构的重叠，得到封闭的学习动力学和泛化曲线。
- **边界：** 它是早期且很简洁的原型分类模型，非现代深网络；适合作为入门练习，看清自平均、在线新样本和序参量闭合各自起什么作用。

## 读论文时留下的五个问题

1. **系统是什么？** 写出数据生成、模型、损失、训练或后验，明确随机性来自哪里。
2. **极限是什么？** $D$、宽度、输入维度、时间、温度各如何增长？极限的先后是否重要？
3. **近似发生在哪里？** 是固定特征、Gaussian 等价、replica 对称鞍点、自平均、动态闭合，还是现象论建模？不要把模型定义和后续近似混为一谈。
4. **哪一个量承载机制？** 可以是重叠、核谱、响应函数、自由能、熵或局部相关；追踪它如何连接到可观测的测试误差或生成轨迹。
5. **怎样证伪或拓展？** 对照有限尺寸、真实数据、不同训练算法和不同数据结构；检查理论失效的位置是否正好指向被忽略的机制。

这些论文分别研究不同模型和极限。跨论文比较时，尤其要区分**贝叶斯最优推断、梯度训练动力学、固定随机特征和现象论**：它们都能产生相图或标度律，但回答的不是同一个问题。
