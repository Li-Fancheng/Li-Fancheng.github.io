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

### 零态稳定性之外还有瞬态放大

线性扰动的精确传播子是 $e^{(-I+gJ)t}$。非对称随机矩阵通常非正规，特征向量并不正交，因此最大实部为负只保证长期衰减，不保证所有初始扰动的范数从第一刻就单调变小。短时输入响应和长期 Lyapunov 稳定性可能不同。

当讨论循环网络作为储备池时，这种差别很实际：系统可以在最终稳定的同时，对某些短时输入有较大增益；增益太大又会放大噪声。仅用谱半径或最长线性时间尺度判断“记忆能力”，没有包含读出方向和信噪比。

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

### 自洽场经过泄漏滤波后才成为状态

稳定的单元解可写为 $x(t)=\int_0^\infty e^{-u}\eta(t-u)du$，所以

$$
\Delta(\tau)=g^2\int_0^\infty\!du\int_0^\infty\!dv\,
e^{-u-v}C_\phi(\tau+u-v).
$$

对两个时间分别作用 $\partial_t+1$ 和 $\partial_{t'}+1$，在平稳条件下 $\partial_t\partial_{t'}\Delta=-\Delta''$，得到 $-\Delta''+\Delta=g^2C_\phi$。这一步的负号来自两个时间差导数方向相反。

特别在 $\tau=0$，卷积可以化为

$$
\Delta_0=g^2\int_0^\infty e^{-\tau}C_\phi(\tau)d\tau.
$$

它依赖整条相关函数，而不是只有 $C_\phi(0)$。静态方差 ansatz 相当于把时间依赖忽略，才会得到类似 $q=g^2\mathbb E\phi(\sqrt qz)^2$ 的方程。这个对照直接显示，时间相关不是在静态平均场上加一个附属指标，而是决定活动幅度本身的组成部分。

### 二元 Gaussian 积分的实际构造

给定 $\Delta_0$ 和 $c=\Delta/\Delta_0$，可取 $u=\sqrt{\Delta_0}z_1$、$v=\sqrt{\Delta_0}[cz_1+\sqrt{1-c^2}z_2]$，将 $C_\phi$ 化为两个独立 Gaussian 的积分。端点 $|c|=1$ 需按连续极限处理。任何候选 $\Delta(\tau)$ 都应满足 $|\Delta(\tau)|\leq\Delta_0$，更强地还应作为一个正定相关函数；仅数值满足一个二阶 ODE，不自动保证它是可实现的随机过程。

## 4. 势能类比与边界条件

对固定 $\Delta_0$ 定义

$$
V(\Delta;\Delta_0)=-\frac{\Delta^2}{2}
+g^2\int_0^\Delta C_\phi(u;\Delta_0)du.
$$

相关方程变成 $\Delta''=-\partial_\Delta V$。对可微平稳过程，$\Delta'(0)=0$。若考虑零均值混合解，还要求 $\Delta(\tau)\to0$、$\Delta'(\tau)\to0$，对应势能匹配 $V(\Delta_0;\Delta_0)=V(0;\Delta_0)$。

这个边界值问题用于确定非零活动方差及相关衰减。若有非零长期相关平台，边界条件应改变，不能强制设为零。

有时看到静态方程 $q=g^2\int Dz\tanh^2(\sqrt qz)$。它属于时间不变随机场的固定点 ansatz；在持续时变的混沌态中，滤波器会混合不同时间，不能用它替代完整 DMFT 方程求活动方差。

### 机械能匹配并不是任意初值积分

将 $\Delta''=-\partial_\Delta V$ 乘以 $\Delta'$，得到守恒量 $\frac12(\Delta')^2+V(\Delta;\Delta_0)$。混合解从 $\Delta(0)=\Delta_0$、$\Delta'(0)=0$ 出发，并在无穷远趋近零，所以必须同时满足两端机械能相同。

$\Delta_0$ 不只是可随意设置的初始振幅，因为它还出现在 $C_\phi(\Delta;\Delta_0)$ 的边缘方差中，因而改变整个势函数。求解流程实际上是：对一个候选方差构造势，再检查能量匹配和相关函数边界；只有自洽方差才能给出物理解。任意从某个方差出发积分出的振荡或发散曲线，不一定代表网络的一种稳态。

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

### 两副本稳定性把 Lyapunov 问题化为谱问题

考虑两条极近轨迹，在相同自洽活动背景上研究它们的交叉相关扰动 $k(t,t')$。在非对称大规模 DMFT 的相应线性化中，Gaussian 协方差求导给出

$$
(\partial_t+1)(\partial_{t'}+1)k(t,t')
=g^2 C_{\phi'}(t-t')k(t,t'),
$$

其中 $C_{\phi'}(\tau)=\mathbb E[\phi'(x(t))\phi'(x(t+\tau))]$。它仍依赖两个时刻，不能只用 $\mathbb E\phi'^2$ 替代。

令中心时间 $T=(t+t')/2$、相对时间 $\tau=t-t'$，并取 $k=e^{\omega T}\psi(\tau)$。利用 $\partial_t\partial_{t'}=\frac14\partial_T^2-\partial_\tau^2$，得到

$$
\left[-\frac{d^2}{d\tau^2}+1-g^2C_{\phi'}(\tau)\right]\psi
=\left(-\omega-\frac{\omega^2}{4}\right)\psi.
$$

因此一个一维 Schrödinger 型算子的最低本征值 $E_0$ 控制最强相关扰动增长。相关扰动对应扰动幅度的平方，故最大 Lyapunov 指数为

$$
\boxed{\lambda_{\max}=-1+\sqrt{1-E_0}.}
$$

这是在所述平稳大规模自洽理论中的关系；有限网络与有限时间估计需要另行比较。$E_0\lt 0$ 时指数为正，零态失稳和持续混沌于是由不同层次的谱问题联系起来。

### 非平凡衰减相关为何对应不稳定轨迹

对相关方程求一次 $\tau$ 导数，并使用 Price 恒等式 $dC_\phi/d\Delta=C_{\phi'}$，得到

$$
\left[-\frac{d^2}{d\tau^2}+1-g^2C_{\phi'}(\tau)\right]\Delta'(\tau)=0.
$$

因此平稳相关的平移导数是零能解。对一条平滑、偶对称且向零衰减的非平凡相关函数，$\Delta'$ 是有节点的奇函数。在适用的全实线束缚态问题中，最低能态无节点，能量低于这个零能奇态，因此 $E_0\lt 0$。这个论证解释了相关函数形状与正 Lyapunov 指数之间的物理联系，而不只是观察一条不规则轨迹后给它贴上混沌标签。

若相关有不衰减平台、外部噪声使原点边界出现尖点，或耦合对称性引入记忆反应项，上述边界和谱论证都要重新检查。随机循环网络的可解性来自完整的模型条件，不是任意复杂反馈系统都能照搬同一个 Schrödinger 算子。

## 6. 数值实验与记忆时间

![随机循环网络的状态轨迹和时间相关函数。](figures/experiment.svg)

*图 1 · 自制有限网络模拟。* 稠密非对称随机矩阵，比较稳定侧与高增益侧；脚本默认 $N=400$、$g=0.8,1.8$、Euler 步长 0.05。左图为单元轨迹，右图为暂态后跨单元平均的未去均值相关 $\Delta(\tau)/\Delta(0)$。稳定侧的幅度会趋近零，因此其归一化相关不能解释成一个有限幅度平稳活动的记忆函数。

复现命令为 `python experiments.py --note 16`，同时保存有限时间切向增长率。该估计依赖步长、运行时长和规模。靠近 $g=1$ 时线性慢时间尺度约为 $(1-g)^{-1}$，而信息记忆还依赖输入、读出和噪声的定义。

## 7. 相关文献与理论范围

**原始文献。** [Sompolinsky、Crisanti 与 Sommers（1988），Chaos in Random Neural Networks](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.61.259)给出经典动力学平均场理论及混沌转变。本文的有限网络轨迹是模型说明，理论极限与数值证据分别陈述。

---

← [Note 15](/notes/neural_network_theory/hopfield/) · [系列目录](/notes/neural_network_theory/) · [Note 17](/notes/neural_network_theory/time_scaling/) →
