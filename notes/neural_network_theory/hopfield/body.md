> **IV · 记忆与集体动力学 · Note 15**
>
> ← [Note 14](/notes/neural_network_theory/sgd/) · [系列目录](/notes/neural_network_theory/) · [Note 16](/notes/neural_network_theory/random_rnn/) →

**可解性层次：** 低负载平均场；广延负载需自洽修正。

## 1. 将记忆存成动力学的稳定状态

联想记忆的任务是：输入一个被破坏的模式，系统经过自身动力学恢复原模式。Hopfield 模型把记忆编码在相互作用里，将回忆过程实现为能量下降。这里没有用输出标签反向传播；学习规则与检索动力学是两个独立步骤。

经典模型非常适合展示序参量：一个微观状态有 $N$ 个自旋，而与某个存储模式的重叠可以直接描述回忆是否成功。

## 2. Hebbian 耦合与异步更新

存储模式 $\xi_i^\mu\in\{-1,1\}$，$\mu=1,\ldots,P$，设

$$
J_{ij}=\frac1N\sum_{\mu=1}^P\xi_i^\mu\xi_j^\mu,
\quad i\neq j,\qquad J_{ii}=0,
$$

$$
E(s)=-\frac12\sum_{i\neq j}J_{ij}s_is_j,
\qquad h_i=\sum_jJ_{ij}s_j.
$$

一次只更新一个自旋 $s_i\leftarrow\operatorname{sign}(h_i)$，$h_i=0$ 时保持原值。若发生翻转，能量变化为

$$
\Delta E=2s_i^{\rm old}h_i\leq0.
$$

有限状态空间与严格能量下降保证这种零温异步动力学最终进入不再翻转的固定点。这个证明使用了对称耦合、无自耦合与单点更新；并行同步更新可以产生二周期，非对称耦合也未必存在同样的能量函数。

## 3. 检索重叠与串扰

定义 $m_\mu=N^{-1}\sum_i\xi_i^\mu s_i$。局部场可重写为

$$
h_i=\sum_\mu\xi_i^\mu m_\mu-\frac PN s_i.
$$

最后一项正是去掉 $J_{ii}$ 的修正。低负载 $P/N\to0$ 时可忽略；广延负载 $P=\alpha N$ 时不能在中间计算中随意丢失。

若系统接近第一个模式，$m_1=O(1)$、其他重叠较小，$\xi_i^1m_1$ 提供信号，其余模式提供串扰。但这些串扰不是独立外加白噪声：自旋受到所有模式影响，又反过来进入所有 $m_\mu$，所以存在反馈相关。

## 4. 低负载单模式的自洽方程

先取 $P=1$、$N\to\infty$，采用符合详细平衡的随机单自旋更新，逆温度为 $\beta$。平衡条件下

$$
\langle s_i\rangle=\tanh(\beta\xi_i m),
\qquad
\boxed{m=\tanh(\beta m).}
$$

也可从每自旋自由能

$$
F(m)=\frac{m^2}2-\frac1\beta\log[2\cosh(\beta m)]
$$

求驻点得到。小 $m$ 展开为

$$
0=(\beta-1)m-\frac{\beta^3}{3}m^3+O(m^5),
$$

所以 $\beta_c=1$，$\beta\gt 1$ 时出现非零检索解，临界附近 $m^2\simeq3(\beta-1)/\beta^3$。这是低负载热噪声转变，不是随存储模式数量变化的容量转变。

## 5. 广延负载下，为什么简单 Gaussian 噪声不够

一种最朴素的近似把串扰看成方差为 $\alpha$ 的 Gaussian，写成 $m=\int Dz\tanh[\beta(m+\sqrt\alpha z)]$。它能说明噪声破坏检索，却忽略模式反馈，不能据此精确算出经典容量。

标准广延负载平衡分析在相应副本对称鞍点中引入 $q,r$：

$$
m=\int Dz\tanh[\beta(m+\sqrt{\alpha r}\,z)],
$$

$$
q=\int Dz\tanh^2[\beta(m+\sqrt{\alpha r}\,z)],
\qquad
r=\frac{q}{[1-\beta(1-q)]^2}.
$$

其中 $Dz=e^{-z^2/2}dz/\sqrt{2\pi}$，$r$ 反映串扰反馈放大。这些方程仍需比较自由能和检查稳定性；检索解存在、检索吸引域足够大、热平衡占优是不同判据。

经典常说的 $\alpha\simeq0.138$ 是独立无偏随机模式、Hebbian 耦合和零温极限下的检索容量尺度，带有具体稳定性判据。它不是所有联想记忆架构的信息容量，更不适用于稀疏模式、相关模式或修改过的学习规则。混合模式等伪吸引子也会占据状态空间。

有检索时某个 $m_\mu\neq0$；顺磁态没有宏观检索且平衡自旋玻璃序参量为零；自旋玻璃态可有冻结而无单个模式检索。有限系统的一次静态重叠，不能独自区分这些热力学相。

## 6. 数值图与实验

![Hopfield 异步检索能量与低负载自洽解。](figures/experiment.svg)

*图 1 · 两个不同极限的自制计算。* 左图是 $N=400,P=12$ 的零温异步检索，初始将模式 1 的 20% 自旋翻转；每个 sweep 随机遍历所有自旋，能量单调不增，本次重叠从 0.6 恢复至 1。右图独立求解低负载方程 $m=\tanh(\beta m)$。右图没有用于拟合左图的有限负载结果。

复现命令为 `python experiments.py --note 15`。有限尺寸下的检索成功率同时依赖负载和初始破坏比例，得到的是操作性容量。它与热力学相图中的解存在、稳定性及自由能判据需要分别讨论。

## 7. 相关文献

**相关文献。** [Hopfield（1982），Neural networks and physical systems with emergent collective computational abilities](https://www.pnas.org/doi/10.1073/pnas.79.8.2554)建立能量吸引子记忆；[Amit、Gutfreund 与 Sompolinsky（1985），Storing Infinite Numbers of Patterns in a Spin-Glass Model of Neural Networks](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.55.1530)分析广延存储。低负载方程与广延负载理论在本篇始终分别标明。

---

← [Note 14](/notes/neural_network_theory/sgd/) · [系列目录](/notes/neural_network_theory/) · [Note 16](/notes/neural_network_theory/random_rnn/) →
