> **III · 深度与优化动力学 · Note 09**
>
> ← [Note 08](/notes/neural_network_theory/lazy_and_mean_field/) · [系列目录](/notes/neural_network_theory/) · [Note 10](/notes/neural_network_theory/implicit_bias/) →

**可解性层次：** 对齐和平衡条件下精确可解。

## 1. 同一个线性函数类，可以有完全不同的学习时间

单层线性回归沿每个模态指数接近目标；把同一个线性映射写成多层矩阵乘积后，优化可以先停在平台，再快速下降。激活函数没有非线性，动力学却已经非线性，因为梯度分别作用在每个因子上。

本篇的精确解有清楚条件：输入白化，学习模态对齐，标量因子处在平衡分支。一般矩阵初始化会发生模态旋转与耦合，不能直接把所有奇异值独立代入 Logistic 方程。

## 2. 两层标量模型与平衡守恒量

令有效权重 $p=w_2w_1$，目标 $\sigma\gt 0$，损失

$$
\mathcal L=\frac12(w_2w_1-\sigma)^2.
$$

单位迁移率的梯度流为

$$
\dot w_1=(\sigma-p)w_2,\qquad
\dot w_2=(\sigma-p)w_1.
$$

因此

$$
\frac d{dt}(w_1^2-w_2^2)=0,\qquad
\dot p=(w_1^2+w_2^2)(\sigma-p).
$$

只有 $w_1^2=w_2^2$ 且选取同号、$p\gt 0$ 分支时，才有 $w_1^2+w_2^2=2p$。于是

$$
\boxed{\dot p=2p(\sigma-p).}
$$

若初始不平衡常数为 $c=w_1^2-w_2^2$，精确的有效方程是 $\dot p=\sqrt{c^2+4p^2}(\sigma-p)$。平台与初始化尺度会受 $c$ 影响；将任意初始化都设成平衡分支会丢掉这个物理量。

## 3. Logistic 解、平台与快速学习

分离变量：

$$
\int\frac{dp}{p(\sigma-p)}=2t+\mathrm{const},
$$

得到 $0\lt p_0\lt \sigma$ 下的精确解

$$
\boxed{p(t)=\frac\sigma{1+(\sigma/p_0-1)e^{-2\sigma t}}.}
$$

早期 $p\ll\sigma$ 时，$p\simeq p_0e^{2\sigma t}$。它的**相对**增长率并不小，但绝对函数变化仍很小，所以损失看起来近乎不动。到 $p=O(\sigma)$ 时，输出快速改变；接近目标后，误差又指数衰减。

达到一半目标的时间为

$$
t_{1/2}=\frac1{2\sigma}\log\left(\frac\sigma{p_0}-1\right),
$$

适用于 $p_0\lt \sigma/2$。强信号 $\sigma$ 较早学会，小初始化使平台按对数延长。若 $w_1=w_2=0$，梯度严格为零，不会自行离开；“极小但非零”与“恰好为零”不同。

例如 $\sigma=1,p_0=10^{-4}$ 时，$t_{1/2}\simeq4.605$；若把 $p_0$ 再缩小 100 倍，平台仅再延长约 $2.303$。这些数值使用的是当前损失与迁移率归一化。

## 4. 多层标量网络与更强的启动瓶颈

深度 $L$、$p=\prod_{\ell=1}^Lw_\ell$，同样使用 $\mathcal L=(p-\sigma)^2/2$。所有 $w_\ell^2-w_k^2$ 在梯度流下守恒。在平衡、正权重分支 $w_\ell=p^{1/L}$ 上，

$$
\boxed{\dot p=L p^{2-2/L}(\sigma-p).}
$$

$L=1$ 是线性松弛，$L=2$ 是 Logistic；$L\gt 2$ 时早期 $\dot p\simeq L\sigma p^{2-2/L}$，从 $p_0$ 长到一个固定的小阈值 $p_c$ 的时间为

$$
t_{\rm esc}\simeq
\frac{p_0^{-(1-2/L)}-p_c^{-(1-2/L)}}
{L\sigma(1-2/L)}.
$$

对 $p_0$ 的依赖从对数变成幂次。比较深度时还应说明固定的是单层初始化 $w_0$ 还是乘积初始化 $p_0=w_0^L$；两种比较给出不同的量级。

## 5. 矩阵模态何时可以独立求解

设输入已白化 $\mathbb E[xx^\top]=I$，输入输出交叉矩阵 $S=\mathbb E[yx^\top]$，两层映射 $M=W_2W_1$。忽略与 $M$ 无关的常数，总体损失为 $\|M-S\|_F^2/2$，梯度流

$$
\dot W_1=W_2^\top(S-W_2W_1),\qquad
\dot W_2=(S-W_2W_1)W_1^\top.
$$

直接相减可验证

$$
\frac d{dt}(W_2^\top W_2-W_1W_1^\top)=0.
$$

若 $S=U\operatorname{diag}(\sigma_i)V^\top$，初始化恰能写成 $W_2=UA$、$W_1=BV^\top$，且 $A,B$ 为匹配的对角矩阵，流保持这一子空间。每个 $p_i=A_{ii}B_{ii}$ 因而独立满足前面的标量方程；平衡时各自具有 Logistic 解。

这给出分阶段学习的可解模型：不同 $\sigma_i$ 产生不同启动时间。一般未对齐初始化、非白输入或瓶颈秩不足，会改变这些方程。精确模态分离是一个有用的理论控制条件，而非所有深度线性训练的默认事实。

## 6. 数值结果

![深度线性网络的 Logistic 学习和深度启动时间。](figures/experiment.svg)

*图 1 · 自制精确解与 ODE。* 左图用 $p_0=10^{-4}$、$\sigma=0.5,1,2$ 展示两层精确解，圆点为独立数值积分。右图固定 $\sigma=1$ 与有效乘积初始化 $p_0=10^{-3}$，比较 $L=1,2,3,4$ 的平衡 ODE；横轴按需使用对数尺度。两幅图没有加入标签噪声或特征旋转。

复现命令为 `python experiments.py --note 9`。完整因子方程中的平衡守恒量提供了独立的数值检查；离散梯度下降只近似守恒，较大学习率还会引入连续时间模型之外的效应。

## 7. 物理解释与相关文献

深度在这里改变参数空间的有效迁移率：$p$ 的移动速度乘上 $L p^{2-2/L}$。靠近小乘积处，函数坐标上的运动很慢。平台可以来自参数化造成的瓶颈，不必意味着模型先在记忆训练样本。

**原始文献。** [Saxe、McClelland 与 Ganguli，Exact solutions to the nonlinear dynamics of learning in deep linear neural networks](https://arxiv.org/abs/1312.6120)是本主题的核心阅读。本文将其可解思想落实到明确的白化、对齐和平衡设定，不宣称任意深度线性矩阵流都具有相同闭式解。

---

← [Note 08](/notes/neural_network_theory/lazy_and_mean_field/) · [系列目录](/notes/neural_network_theory/) · [Note 10](/notes/neural_network_theory/implicit_bias/) →
