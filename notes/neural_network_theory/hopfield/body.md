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

### 模式重叠为何成为集体变量

将 Hebbian 耦合代入能量，得到

$$
E(s)=-\frac N2\sum_{\mu=1}^P m_\mu^2+\frac P2.
$$

因此与任一模式形成宏观重叠都会降低能量。但同时与几个模式部分对齐也可能降低能量，能量式本身没有保证最终一定落到某个原始存储模式。混合态与伪吸引子是模型结构的一部分，而不是数值算法出了错。

检索的动力学意义是一个吸引域：许多不同的破坏输入落入同一个最终状态。能量低、固定点稳定和吸引域大是不同性质。一个存储模式即使对应局部稳定态，也可能只能纠正很少的位翻转；若实验任务要求从 30% 的破坏中恢复，其操作性容量会比只要求“存在检索解”更严格。

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

### 单模式的配分函数可以完整降到一维

对 $P=1$，作 gauge 变换 $\tau_i=\xi_i s_i$。独立随机模式被消去，能量成为 Curie–Weiss 形式

$$
E=-\frac1{2N}\left(\sum_i\tau_i\right)^2+\frac12.
$$

利用 Gaussian 积分恒等式，将平方耦合线性化：

$$
Z=e^{-\beta/2}\sqrt{\frac{N\beta}{2\pi}}
\int dm\,e^{-N\beta m^2/2}[2\cosh(\beta m)]^N.
$$

在 $N\to\infty$ 时用 Laplace 方法，配分函数由指数最大的 $m$ 控制，因此每自旋自由能就是前面的 $F(m)$ 的极小值。这个计算明确展示了熵与能量的竞争：$m^2/2$ 项来自集体场的代价，$\log\cosh$ 项来自给定场下独立自旋的状态求和。

在 $m=0$ 附近，

$$
F(m)=\mathrm{const}+\frac{1-\beta}{2}m^2
+\frac{\beta^3}{12}m^4+O(m^6).
$$

二次系数在 $\beta=1$ 变号，四次系数为正，给出连续转变。零场下正、负 $m$ 两个检索态对称；有限系统若在二者之间充分采样，总体平均 $\langle m\rangle$ 可以为零，但 $\langle|m|\rangle$ 不为零。研究自发对称破缺时，取热力学极限、外场趋零与长时间平均的顺序必须区分。

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

### 串扰的反馈如何改变方差

若把第一个模式当作信号，其他 $P-1$ 个小重叠各约为 $N^{-1/2}$，总串扰幅度可以为 $O(1)$。最初的 Gaussian 估计给出方差约 $\alpha$；但每个非检索模式又会扰动自旋，再通过自旋响应增强它自己的重叠。

在相应平均场响应近似中，记 $U=\beta(1-q)$，一次反馈把重叠乘上 $U$，反复响应形成 $1/(1-U)$ 的放大。方差因而乘上其平方，得到 $r=q/(1-U)^2$。这解释了自洽方程中的反馈分母从何而来，也说明为什么将串扰直接视为外加独立噪声会错估容量。严格推导还需要处理 quenched 平均与反应项，不能仅靠这个几何级数图像代替。

### 零温检索分支与 0.138 的来源

在零温检索分支上 $q\to1$，但 $U=\beta(1-q)$ 可以保持有限。有效场 $h=m+\sqrt{\alpha r}\,z$ 的符号平均给出

$$
m=\operatorname{erf}\!\left(\frac m{\sqrt{2\alpha r}}\right).
$$

另一方面，$\beta\operatorname{sech}^2(\beta h)$ 在分布意义下趋于 $2\delta(h)$，因此

$$
U=\sqrt{\frac2{\pi\alpha r}}
e^{-m^2/(2\alpha r)},\qquad r=(1-U)^{-2}.
$$

令 $y=m/\sqrt{2\alpha r}$，先由 $m=\operatorname{erf}y$ 得到
$U=2ye^{-y^2}/(\sqrt\pi\operatorname{erf}y)$，再消去 $r$，得到单参数检索分支

$$
\boxed{\alpha(y)=\left[
\frac{\operatorname{erf}y}{\sqrt2\,y}
-\sqrt{\frac2\pi}e^{-y^2}\right]^2.}
$$

这一分支的最大负载约为 $0.138$，对应经典零温检索解消失的容量尺度。这个数来自信号、串扰与反馈的自洽竞争，不是简单用每个模式占用 $N$ 个比特对总权重数作除法。若研究热平衡全局占优或更严格的稳定性区域，还需比较其他态与相应稳定性条件，不能把单参数分支的存在当作全部相图。

### 近似检索与精确存储每一位的不同容量

在一个原始模式上，忽略反馈的初步估计把局部场视为 $\xi_i^1h_i\simeq1+\sqrt\alpha z$，单个位不稳定概率约为 $\Phi(-1/\sqrt\alpha)$。对于固定 $\alpha\gt 0$，这个概率虽可很小，却不随 $N$ 消失，因此原模式并不必然逐位都是精确固定点。

如果要求一个给定模式的所有 $N$ 位都不出错，至少要让 $N\Phi(-1/\sqrt\alpha)$ 变小。Gaussian 尾部提示负载需按 $1/\log N$ 的量级缩小，远比允许少量错误但有宏观重叠的检索要求严格。这里是局部场近似的尺度判断，不把其前因子当作所有模式联合稳定的严格容量。

这一区分消除了两个表面矛盾：$P=O(N)$ 的联想检索能力，与更严格精确模式稳定性下较小的容量，可以同时成立。所谓“记住一个模式”，必须说明允许多少错误、从多大破坏恢复，以及观察多长时间。

![零温 Hopfield 检索分支的负载与重叠。](figures/retrieval_capacity.svg)

*补充图 · RS 检索分支的显式计算。* 左图直接计算 $\alpha(y)$，最大值为 0.137906；右图将同一分支参数化为 $(\alpha,m)$，转折点重叠约为 0.9674。曲线包括不同分支，不将整条曲线都视为稳定检索态。这里分析的是指定零温自洽方程的解结构，完整热平衡相图还需要自由能与稳定性判据。

## 6. 数值图与实验

![Hopfield 异步检索能量与低负载自洽解。](figures/experiment.svg)

*图 1 · 两个不同极限的自制计算。* 左图是 $N=400,P=12$ 的零温异步检索，初始将模式 1 的 20% 自旋翻转；每个 sweep 随机遍历所有自旋，能量单调不增，本次重叠从 0.6 恢复至 1。右图独立求解低负载方程 $m=\tanh(\beta m)$。右图没有用于拟合左图的有限负载结果。

复现命令为 `python experiments.py --note 15`。有限尺寸下的检索成功率同时依赖负载和初始破坏比例，得到的是操作性容量。它与热力学相图中的解存在、稳定性及自由能判据需要分别讨论。

## 7. 相关文献

**相关文献。** [Hopfield（1982），Neural networks and physical systems with emergent collective computational abilities](https://www.pnas.org/doi/10.1073/pnas.79.8.2554)建立能量吸引子记忆；[Amit、Gutfreund 与 Sompolinsky（1985），Storing Infinite Numbers of Patterns in a Spin-Glass Model of Neural Networks](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.55.1530)分析广延存储。低负载方程与广延负载理论在本篇始终分别标明。

---

← [Note 14](/notes/neural_network_theory/sgd/) · [系列目录](/notes/neural_network_theory/) · [Note 16](/notes/neural_network_theory/random_rnn/) →
