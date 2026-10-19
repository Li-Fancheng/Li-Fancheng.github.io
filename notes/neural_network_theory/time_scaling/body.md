> **V · Scaling Laws · Note 17**
>
> ← [Note 16](/notes/neural_network_theory/random_rnn/) · [系列目录](/notes/neural_network_theory/) · [Note 18](/notes/neural_network_theory/data_scaling/) →

**可解性层次：** 固定谱模型精确；连续谱为渐近。

## 1. 指数模态怎样叠加成幂律

一个模态按指数衰减，不代表总损失也只能是单个指数。若模型含有跨越许多量级的松弛速率，快模态先消失、慢模态逐渐接管，总损失可以在一段很宽的时间窗口里近似幂律。

本篇从精确固定谱动力学出发，推导损失加权速率谱与时间指数的关系。它解释优化曲线的形状，不直接给出样本规模或模型规模的泛化指数。

## 2. 有限模态的精确表达

在线性回归、固定随机特征或冻结 NTK 下，残差可写成

$$
r_k(t)=a_ke^{-\lambda_kt},\qquad
\mathcal L(t)=\frac12\sum_k a_k^2e^{-2\lambda_kt}.
$$

常数 $1/D$ 可吸收进 $a_k^2$；必须在整篇计算中保持一致。零模态给出永久平台 $\mathcal L_0=\frac12\sum_{\lambda_k=0}a_k^2$。

定义损失加权的谱测度

$$
d\mu(\lambda)=\sum_{\lambda_k\gt 0}a_k^2\delta(\lambda-\lambda_k)d\lambda.
$$

于是 $\mathcal L(t)-\mathcal L_0=\frac12\int e^{-2\lambda t}d\mu(\lambda)$，是谱测度的 Laplace 变换。仅知道特征值计数密度还不够：初始化或目标在各模态上的能量同样进入权重。

## 3. 连续谱小速率端决定长时间指数

假设相关窗口内 $d\mu/d\lambda\simeq A\lambda^{\nu-1}$，$\nu\gt 0$。若下边缘为零且其他部分可控，则

$$
\mathcal L(t)-\mathcal L_0\simeq
\frac A2\int_0^\infty\lambda^{\nu-1}e^{-2\lambda t}d\lambda.
$$

令 $u=2\lambda t$，得到

$$
\boxed{\mathcal L(t)-\mathcal L_0
\sim\frac{A\Gamma(\nu)}{2^{\nu+1}}t^{-\nu}.}
$$

指数来自小速率端的**加权**谱。若计数密度 $\rho(\lambda)\sim\lambda^a$、每模态平均初始能量 $\mathbb E[a_k^2\mid\lambda]\sim\lambda^b$，则 $\nu=a+b+1$。忽略目标权重会把时间指数算错。

当 $\nu\leq0$ 时，该连续密度在零端不可积；有限初始损失需要截断或不同权重，不能沿用这个公式再声称得到负指数增长。

## 4. 有限窗口、谱隙与极限顺序

若谱仅在 $[\lambda_{\min},\lambda_{\max}]$ 上近似该幂律，那么幂律通常需要

$$
\lambda_{\max}^{-1}\ll t\ll\lambda_{\min}^{-1}.
$$

太早时，所有模态都未明显变化，损失接近初值；太晚时，谱下边缘或最慢单模态支配。连续密度若在正下边缘附近为 $B(\lambda-\lambda_{\min})^a$，长时间会出现 $e^{-2\lambda_{\min}t}t^{-(a+1)}$，而非无限延续的纯幂律。

任何有限正定矩阵最终都有最慢指数。先取系统规模无穷、再取时间无穷，才可能得到无谱隙幂律；反过来则一般得到最慢模态指数。两个极限不应默默交换。

## 5. 用模态编号也能推导指数

若 $\lambda_k\sim k^{-p}$、$a_k^2\sim k^{-q}$，$p\gt 0,q\gt 1$，则

$$
\mathcal L(t)\simeq\frac12\int_1^\infty
k^{-q}e^{-2t k^{-p}}dk.
$$

变量替换 $u=2t k^{-p}$ 给出

$$
\boxed{\mathcal L(t)\sim
\frac{\Gamma((q-1)/p)}{2p}(2t)^{-(q-1)/p}.}
$$

直观上，时间 $t$ 大致学会 $k\lesssim t^{1/p}$ 的快模态，剩余尾部 $\sum_{k\gt t^{1/p}}k^{-q}\sim t^{-(q-1)/p}$。精确积分与这个移动截止尺度给出同一指数。

同一核谱的时间指数还依赖目标系数；后续数据指数又包含估计噪声与正则化，二者不能只因共享一个谱就直接相等。

## 6. 数值图与诊断

![连续松弛速率谱产生的损失幂律及截断效应。](figures/experiment.svg)

*图 1 · 自制积分。* 左图取 $d\mu/d\lambda=\lambda^{\nu-1}$、$0\lt \lambda\lt 1$，$\nu=0.5,1,2$，用不完全 Gamma 函数计算完整曲线；右图在 $\nu=1$ 时加入 $\lambda_{\min}=10^{-3}$，用 1500 点谱积分比较无下截断结果。幂律窗口和最终偏离都来自同一个谱模型。

复现命令为 `python experiments.py --note 17`。$\nu=1$ 时积分简化为 $\mathcal L(t)=(1-e^{-2t})/(4t)$，短时极限为 $1/2$，长时为 $1/(4t)$，与完整数值曲线一致。

实验拟合应同时报告局部斜率 $\beta_{\rm eff}(t)=-d\log(\mathcal L-\mathcal L_0)/d\log t$。若未知平台未减去、指数窗过短或时间横轴采用更新数而学习率在变化，直线拟合可能给出没有机制意义的指数。

## 7. 特征学习怎样超出固定谱解释

若核随时间改变，模态会旋转，能量在方向间转移，不能一般写成同一组 $a_k^2e^{-2\lambda_kt}$。即使瞬时谱变化缓慢，还需控制基底旋转和非对易性。观察到幂律既不能证明固定谱机制，也不能证明特征学习；两者需要独立测量区分。

**相关文献。** 固定谱优化可与 [Advani 与 Saxe](https://arxiv.org/abs/1710.03667)的高维学习动力学联系；更广的幂律机制比较见 [Bahri 等，Explaining Neural Scaling Laws](https://arxiv.org/abs/2102.06701)。本篇的 Laplace 推导独立给出，不把某个时间指数宣称为所有网络的普适定律。

---

← [Note 16](/notes/neural_network_theory/random_rnn/) · [系列目录](/notes/neural_network_theory/) · [Note 18](/notes/neural_network_theory/data_scaling/) →
