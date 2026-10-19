> **IV · 记忆与集体动力学 · Note 16**
>
> ← [Note 15](/notes/neural_network_theory/hopfield/) · [系列目录](/notes/neural_network_theory/) · [Note 17](/notes/neural_network_theory/time_scaling/) →

**可解性层次：** 非对称稠密随机耦合的无限规模极限。

## 1. 随机连接也会产生集体时间结构

循环网络的状态不断反馈，即使不接受外部输入，也可能衰减到静止态、进入其他固定点，或持续表现不规则活动。区分这些行为需要动力学稳定性、时间相关函数和扰动增长率；单看一条复杂曲线不能确认混沌。

本篇讨论经典的稠密、非对称 Gaussian 耦合网络。耦合若改成对称、稀疏或具有强结构，平均场与转变都可能改变。

## 2. 最小模型与零态失稳

令单元时间常数为 1：

$$
\dot x_i=-x_i+g\sum_{j=1}^N J_{ij}\phi(x_j),
\qquad J_{ij}\sim\mathcal N(0,1/N),
\quad \phi(x)=\tanh x.
$$

不同有向边独立，特别是 $J_{ij}$ 与 $J_{ji}$ 不相关。零态总是解，线性化为

$$
\delta\dot x=(-I+gJ)\delta x.
$$

大 $N$ 时 $J$ 的复特征值铺满单位圆盘，最大实部趋近 1。因此 $g\lt 1$ 的零态线性稳定，$g\gt 1$ 零态失稳。这个阈值仅是零态稳定性分析；单靠圆律没有证明失稳后一定混沌。

## 3. 从多体反馈到自洽 Gaussian 过程

对非对称稠密随机耦合，在合适的 $N\to\infty$ 极限中，典型单元可描述为

$$
(\partial_t+1)x(t)=\eta(t),\qquad
\mathbb E[\eta(t)\eta(t')]=g^2 C_\phi(t,t'),
$$

$$
C_\phi(t,t')=\mathbb E[\phi(x(t))\phi(x(t'))].
$$

$\eta$ 是由网络自洽产生的有色 Gaussian 场。它不是任意给定的外部白噪声；其协方差反过来取决于单元自身的统计。高度非对称性使反馈的某些反应项在该极限下消失，对称耦合通常需要额外记忆项。

若稳态均值为零且时间平移不变，记 $\Delta(\tau)=\mathbb E[x(t)x(t+\tau)]$。对相关函数作用两次线性滤波算子，得到

$$
\boxed{\Delta''(\tau)=\Delta(\tau)-g^2 C_\phi(\tau).}
$$

其中 $C_\phi(\tau)$ 由方差 $\Delta_0$、协方差 $\Delta(\tau)$ 的二元 Gaussian 积分给出：

$$
C_\phi(\tau)=\mathbb E_{u,v}\big[\tanh u\tanh v\big],
\quad \operatorname{Cov}(u,v)=
\begin{pmatrix}\Delta_0&\Delta(\tau)\\\Delta(\tau)&\Delta_0\end{pmatrix}.
$$

因此未知量是整个函数 $\Delta(\tau)$ 加上其初值 $\Delta_0$，不只是一个静态方差。

## 4. 势能类比与边界条件

对固定 $\Delta_0$ 定义

$$
V(\Delta;\Delta_0)=-\frac{\Delta^2}{2}
+g^2\int_0^\Delta C_\phi(u;\Delta_0)du.
$$

相关方程变成 $\Delta''=-\partial_\Delta V$。对可微平稳过程，$\Delta'(0)=0$。若考虑零均值混合解，还要求 $\Delta(\tau)\to0$、$\Delta'(\tau)\to0$，对应势能匹配 $V(\Delta_0;\Delta_0)=V(0;\Delta_0)$。

这个边界值问题用于确定非零活动方差及相关衰减。若有非零长期相关平台，边界条件应改变，不能强制设为零。

有时看到静态方程 $q=g^2\int Dz\tanh^2(\sqrt qz)$。它属于时间不变随机场的固定点 ansatz；在持续时变的混沌态中，滤波器会混合不同时间，不能用它替代完整 DMFT 方程求活动方差。

## 5. 混沌要用扰动增长检验

对同一条轨迹的微小扰动，

$$
\delta\dot x_i=-\delta x_i+
g\sum_jJ_{ij}\phi'(x_j(t))\delta x_j.
$$

最大 Lyapunov 指数定义为适当线性极限下的长期平均增长率

$$
\lambda_{\max}=\lim_{T\to\infty}\frac1T
\log\frac{\|\delta x(T)\|}{\|\delta x(0)\|}.
$$

数值上应反复归一化切向量并累加对数伸长，避免扰动饱和。正指数表示敏感依赖初值；相关衰减、宽频谱和复杂轨迹是辅助现象，不能替代这一检验。有限时间接近零的指数也可能来自周期轨道、慢暂态或有限尺寸，需改变时间、步长和规模。

## 6. 数值实验与记忆时间

![随机循环网络的状态轨迹和时间相关函数。](figures/experiment.svg)

*图 1 · 自制有限网络模拟。* 稠密非对称随机矩阵，比较稳定侧与高增益侧；脚本默认 $N=400$、$g=0.8,1.8$、Euler 步长 0.05。左图为单元轨迹，右图为暂态后跨单元平均的未去均值相关 $\Delta(\tau)/\Delta(0)$。稳定侧的幅度会趋近零，因此其归一化相关不能解释成一个有限幅度平稳活动的记忆函数。

复现命令为 `python experiments.py --note 16`，同时保存有限时间切向增长率。该估计依赖步长、运行时长和规模。靠近 $g=1$ 时线性慢时间尺度约为 $(1-g)^{-1}$，而信息记忆还依赖输入、读出和噪声的定义。

## 7. 相关文献与理论范围

**原始文献。** [Sompolinsky、Crisanti 与 Sommers（1988），Chaos in Random Neural Networks](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.61.259)给出经典动力学平均场理论及混沌转变。本文的有限网络轨迹是模型说明，理论极限与数值证据分别陈述。

---

← [Note 15](/notes/neural_network_theory/hopfield/) · [系列目录](/notes/neural_network_theory/) · [Note 17](/notes/neural_network_theory/time_scaling/) →
