> **III · 深度与优化动力学 · Note 10**
>
> ← [Note 09](/notes/neural_network_theory/deep_linear/) · [系列目录](/notes/neural_network_theory/) · [Note 11](/notes/neural_network_theory/nonlinear_teacher_student/) →

**可解性层次：** 线性模型定理；因子化结论依赖条件。

## 1. 训练误差为零以后，优化还在选择什么

过参数化平方回归有无穷多个插值解；线性可分分类的 logistic 损失则没有有限范数的最小点。两种情形都需要问：优化算法沿什么方向前进，它究竟在解空间中偏好什么？

隐式偏置指目标函数未显式写出正则项时，由算法、参数化、初始化和损失共同造成的解选择。这个词没有指定某个普遍范数，更不自动意味着偏置总有利于真实任务。

## 2. 平方回归：为何零初始化选择最小欧氏范数

对 $\mathcal L=\|Xw-y\|^2/(2D)$，梯度始终在 $X$ 的行空间，故零空间分量不变。若从零出发并收敛，

$$
\hat w=X^+y.
$$

当数据可插值时，这是满足 $Xw=y$ 的欧氏最小范数解。证明很直接：任意其他解为 $X^+y+z$，$z\in\ker X$，两项正交，所以范数平方相加。

若用预条件更新 $\dot w=-P\nabla\mathcal L$，$P\succ0$ 固定，从零出发选择的通常是 $w^\top P^{-1}w$ 最小的插值解。算法的度量已经改变，不能继续沿用相同的欧氏结论。非零初始化则保留对应不可观测分量。

### 初始化与优化度量同时进入解选择

用 $w=P^{1/2}z$ 重写预条件流，可将其变为 $z$ 上的欧氏梯度流。最小 $\|z\|^2$ 等价于最小 $w^\top P^{-1}w$。如果学习率矩阵依时间或参数而变，这个固定度量变换一般不存在，不能将 Adam 或自适应预条件方法直接归结成同一个范数最小化问题。

非零初始化下，平方回归的长期解是 $X^+y+(I-X^+X)w_0$。因此“梯度下降有最小范数偏置”完整的表述还应包括：针对哪一组参数坐标，使用什么度量，从什么初值出发，是否达到插值，以及损失是否在有限点取到最小值。省掉其中任一项，都可能把一个正确的特例变成错误的普遍陈述。

## 3. 可分 logistic 分类：最优点在无穷远

标签 $y_\mu\in\{-1,1\}$，模型无单独偏置项，损失

$$
\mathcal L(w)=\sum_\mu\log(1+e^{-y_\mu w^\top x_\mu}).
$$

若存在 $v$ 使所有 $y_\mu v^\top x_\mu\gt 0$，则沿 $w=cv$、$c\to\infty$，损失趋于零。不存在有限 $w$ 达到这个下确界。梯度流为

$$
\dot w=\sum_\mu\frac{y_\mu x_\mu}{1+e^{y_\mu w^\top x_\mu}}.
$$

大间隔点的梯度指数变小，靠近当前分类边界的样本主导长期更新。因此训练分类准确率已经为 100%，参数方向仍可缓慢转动。

## 4. 最大间隔方向的渐近结构

定义硬间隔解

$$
\hat w=\arg\min_w\frac12\|w\|^2
\quad\text{s.t.}\quad y_\mu w^\top x_\mu\geq1.
$$

其 KKT 条件为

$$
\hat w=\sum_\mu\alpha_\mu y_\mu x_\mu,\quad
\alpha_\mu\geq0,\quad
\alpha_\mu(y_\mu\hat w^\top x_\mu-1)=0.
$$

用渐近 ansatz $w(t)=\hat w\log t+\rho(t)$ 理解其来源。支持向量的梯度约为 $t^{-1}$，非支持点的梯度约为 $t^{-m_\mu}$，其中 $m_\mu=y_\mu\hat w^\top x_\mu\gt 1$，衰减更快。主导向量只能由支持向量的非负组合构成，与 KKT 条件一致。

这段配平给出机制解释，而完整定理还要控制余项。对可分数据、相应指数尾损失和合适步长的梯度下降，可证明方向趋向 $\hat w/\|\hat w\|$；普通情形余项较温和，退化数据可能出现 $\log\log t$ 等修正。因此不要把“余项永远有界”写成一般定理。

范数主项为 $O(\log t)$，方向误差通常只按对数速度缩小。分类准确率的饱和远早于方向收敛，解释了继续训练仍可能改变测试预测，但不保证所有测试分布都因此改善。

### 指数尾怎样选择支持向量

令 $z_\mu=y_\mu x_\mu$。假设渐近展开 $w(t)=\hat w\log t+\rho(t)$，且先考虑 $\rho$ 有界的非退化情形。第 $\mu$ 个梯度项为

$$
\frac{z_\mu}{1+e^{z_\mu^\top w}}
\simeq z_\mu\,t^{-z_\mu^\top\hat w}
e^{-z_\mu^\top\rho}.
$$

左边的参数速度主项是 $\hat w/t$。要使右边与它配平，最小的 $z_\mu^\top\hat w$ 必须归一化到 1。大于 1 的项积分后仅产生次领先修正，等于 1 的项则构成主导组合 $\hat w=\sum_{\mu\in\mathcal S}c_\mu z_\mu$，其中 $c_\mu\geq0$。这恰好是最大间隔 KKT 条件的活跃集结构。

退化数据中，某些支持向量的系数和几何约束不能由一个有限 $\rho$ 同时满足，需要更慢发散的余项。这不改变方向收敛的核心结论，却会显著影响有限时间拟合。将有限步数下的方向差异解释成不同隐式偏置之前，需要先排除这种对数级收敛造成的长暂态。

### 最大间隔方向不等于任意分布下最优分类器

最大间隔是相对训练样本和指定参数度量的几何性质。测试分布若在不同区域集中，或者教师边界与最大间隔方向不同，继续增大训练间隔未必单调降低测试误差。该定理解释算法选择哪条分离方向，并没有独立解决这条方向为何适合真实数据的统计问题。

加入可训练偏置、不同样本权重、非可分噪声、显式正则化或非指数尾损失，都可能改变长期结构。特别是有噪声且不可分时，logistic 最优点可能是有限的，参数范数不再按 $\log t$ 发散；把可分渐近套到所有分类训练尾部，会失去问题本身的条件。

## 5. 参数化改变了有效偏置

考虑正坐标因子化 $\beta_i=u_i^2$。若损失只通过 $\beta$ 依赖 $u$，则

$$
\dot u_i=-2u_i\partial_{\beta_i}\mathcal L,
\qquad
\dot\beta_i=-4\beta_i\partial_{\beta_i}\mathcal L.
$$

有效坐标的迁移率与自身大小成正比，小坐标增长慢；这与直接优化 $\beta$ 的欧氏梯度流不同。若允许正负分支、改变初始化或增加层数，偏置还会改变，所以一个因子化例子不足以证明通用稀疏性定理。

矩阵分解还有一个常用的变分恒等式：

$$
\|M\|_*=\inf_{UV^\top=M}
\frac12(\|U\|_F^2+\|V\|_F^2).
$$

用 $M=A\Sigma B^\top$ 的 SVD，取 $U=A\Sigma^{1/2}$、$V=B\Sigma^{1/2}$ 得到等号；反向不等式来自核范数与因子 Frobenius 范数的界。这个式子说明**显式**惩罚因子范数与核范数的关系，却不能直接证明无正则化因子梯度下降总选择最小核范数解。初始化、观测算子、可交换性及极限顺序都可能重要。

### 一个能精确看见偏置连续变化的有符号对角模型

取 $\beta_i=u_i^2-v_i^2$，损失仅通过线性预测 $X\beta$ 依赖参数，记 $g_i=\partial_{\beta_i}\mathcal L$。欧氏梯度流为 $\dot u_i=-2u_ig_i$、$\dot v_i=2v_ig_i$，因此 $c_i=u_iv_i$ 守恒。

若初始化 $u_i=v_i=\epsilon_i\gt 0$，则 $c_i=\epsilon_i^2$，初始预测为零，但参数并非零。利用 $(u_i^2+v_i^2)^2=\beta_i^2+4c_i^2$，有效流为

$$
\boxed{\dot\beta_i=-4\sqrt{\beta_i^2+4c_i^2}\,g_i.}
$$

定义严格凸势函数

$$
\Phi_c(\beta)=\frac14\sum_i\left[
\beta_i\operatorname{arsinh}\!\left(\frac{\beta_i}{2c_i}\right)
-\sqrt{\beta_i^2+4c_i^2}\right].
$$

它的二阶导数为 $\partial_i^2\Phi_c=1/[4\sqrt{\beta_i^2+4c_i^2}]$，所以流可以改写为

$$
\frac d{dt}\nabla\Phi_c(\beta)=-\nabla_\beta\mathcal L.
$$

这是一种 mirror 型动力学：在原坐标中的状态依赖迁移率，变成对偶坐标 $\nabla\Phi_c$ 的线性受迫运动。对于平方回归，$\nabla\mathcal L=X^\top(X\beta-y)/D$，因此对偶坐标变化始终落在 $X$ 的行空间。

若轨迹收敛到一个有限插值解，令 $\beta_0=0$，则其极限满足

$$
X\beta_\infty=y,\qquad
\nabla\Phi_c(\beta_\infty)-\nabla\Phi_c(\beta_0)
\in\operatorname{im}X^\top.
$$

这正是约束最小化 Bregman 距离
$D_{\Phi_c}(\beta,\beta_0)=\Phi_c(\beta)-\Phi_c(\beta_0)-\nabla\Phi_c(\beta_0)^\top(\beta-\beta_0)$
的 KKT 条件。这里解选择可以完整推出来，而不是只给出一个“倾向稀疏”的定性判断；但仍须保留插值收敛这个前提。

### 小初始化与大初始化的两端

当全部 $c_i=c$ 很大、相关解满足 $|\beta_i|\ll c$ 时，$\Phi_c$ 在原点附近近似 $\sum_i\beta_i^2/(16c)$，因而接近 L2 型解选择。小初始化 $c=\epsilon^2\to0$ 时，对固定非零 $\beta_i$，

$$
\Phi_c(\beta)-\Phi_c(0)
\simeq\frac14\sum_i|\beta_i|
\left[\log\frac{|\beta_i|}{c}-1\right].
$$

若所有坐标以同一速度缩小，领先项是 $\frac14\log(1/c)\|\beta\|_1$，表现出 L1 型选择；余项仍会在多个 L1 最优解之间作进一步区分。坐标初始化若不均匀，则对应不同的有效权重。

函数类始终是线性预测，但初始化尺度已经能把偏置从近似 L2 连续推向稀疏选择。这一例子与 [Woodworth 等，Kernel and Rich Regimes in Overparametrized Models](https://arxiv.org/abs/2002.09277)讨论的尺度依赖机制相关。它同时解释为什么“同架构、同损失、同训练误差”不足以决定最终泛化。

![有符号对角模型的初始化尺度与最终插值解。](figures/initialization_bias.svg)

*补充图 · 初始化尺度控制解选择。* 约束为 $\beta_1+2\beta_2=1$，两个因子均从 $u_i=v_i=\epsilon$ 出发。曲线由 mirror KKT 方程求得，另用原始有效梯度流独立核对，最大差约 $7.2\times10^{-10}$。大初始化接近 L2 解 $(0.2,0.4)$，极小初始化向 L1 解 $(0,0.5)$ 移动。有限小初始化仍有次领先修正，并不在某个有限阈值突然等于 L1 解。

## 6. 数值图与实验

![可分 logistic 梯度流的损失、范数与最大间隔方向。](figures/experiment.svg)

*图 1 · 自制二维可分数据实验。* 用三种已吸收标签的向量 $z_\mu=y_\mu x_\mu=(1,0),(1,1),(2,-0.2)$。硬间隔解为 $(1,0)$，因为第一约束要求 $w_1\geq1$，其余约束在该点满足。脚本数值积分 logistic 梯度流；左图比较损失和参数范数，右图测量与 $(1,0)$ 的方向夹角。横轴为优化时间，未加入权重衰减。

复现命令为 `python experiments.py --note 10`。训练准确率饱和后，参数方向仍继续变化。加入非零 L2 正则化后，最优点变为有限范数，长期行为不再是无正则化的 $\log t$ 发散。

## 7. 相关文献与理论边界

**相关文献。** [Soudry 等，The Implicit Bias of Gradient Descent on Separable Data](https://arxiv.org/abs/1710.10345)给出线性分类方向收敛理论；[Gunasekar 等，Implicit Regularization in Matrix Factorization](https://arxiv.org/abs/1705.09280)提出并分析矩阵因子化的偏置。后者原文包含猜想与特定证据，不应被简化成任意矩阵观测下的普遍最小核范数定理。

---

← [Note 09](/notes/neural_network_theory/deep_linear/) · [系列目录](/notes/neural_network_theory/) · [Note 11](/notes/neural_network_theory/nonlinear_teacher_student/) →
