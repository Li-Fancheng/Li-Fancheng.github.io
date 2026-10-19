> **II · 无限宽与特征学习 · Note 06**
>
> ← [Note 05](/notes/neural_network_theory/gaussian_process/) · [系列目录](/notes/neural_network_theory/) · [Note 07](/notes/neural_network_theory/kernel_learning_curves/) →

**可解性层次：** 固定核动力学精确；冻结核为条件性极限。

## 1. 非线性网络与线性动力学并不矛盾

一个网络可以是输入的非线性函数，却在训练期间只沿初始化处的切线特征变化。此时预测的演化服从线性方程。真正需要解释的是：哪一步永远成立，哪一步依赖无限宽及参数化，而不是把所有网络直接称为核机器。

本篇先写精确的时变核方程，再明确加入冻结核假设。这样能看见理论从哪里进入，也能为有限网络设计诊断。

## 2. 精确恒等式：训练集上的时变 NTK

对标量输出 $f_\theta(x)$，平方损失为

$$
\mathcal L(\theta)=\frac1{2D}\sum_{\mu=1}^D
[f_\theta(x_\mu)-y_\mu]^2.
$$

记 $J_{\mu a}=\partial_{\theta_a}f_\theta(x_\mu)$，残差 $r=f_X-y$。梯度流满足 $\dot\theta=-J^\top r/D$，链式法则给出

$$
\boxed{\dot r=-K_t r,\qquad K_t=\frac1D J_tJ_t^\top.}
$$

任意向量 $v$ 满足 $v^\top K_tv=\|J_t^\top v\|^2/D\geq0$，所以训练损失单调不增：$\dot{\mathcal L}=-r^\top K_tr/D$。但只要 $K_t$ 随参数改变，这仍然不是闭合的常系数线性系统。

文献有时把 $JJ^\top$ 本身叫 NTK，此时方程右侧多一个 $1/D$。本文把它吸收进 $K_t$，所有时间尺度都以这个约定为准。

## 3. 切线模型与冻结核条件

在初始化处展开

$$
f_\theta(x)\simeq f_{\theta_0}(x)+
\nabla_\theta f_{\theta_0}(x)^\top(\theta-\theta_0).
$$

如果训练过程中 Jacobian 的变化可以忽略，便得到 $K_t\simeq K_0$。适当 NTK 参数化下，固定数据集、固定深度、足够宽网络可以在相应条件下趋向这个极限。其成立需要控制参数移动、Jacobian 稳定性、训练时间和核的退化；“宽度很大”不是一个脱离这些条件的独立证明。

例如两层网络 $f=N^{-1/2}\sum_j a_j\phi(w_j^\top x/\sqrt d)$，训练 $a,w$ 时，未除以 $D$ 的极限核为

$$
\Theta(x,x')=\mathbb E[\phi(u)\phi(v)]
+\frac{x^\top x'}d\mathbb E[a^2\phi'(u)\phi'(v)].
$$

第一项来自读出权重，第二项来自隐藏权重。因此它通常不同于仅有前向协方差的 NNGP。激活非光滑时，导数与收敛条件也需相应处理。

## 4. 残差的精确谱解

在冻结核模型内，分解 $K_0=U\operatorname{diag}(\kappa_i)U^\top$，便有

$$
r(t)=e^{-K_0t}r(0),\qquad
r_i(t)=e^{-\kappa_it}r_i(0),
$$

$$
\boxed{\mathcal L(t)=\frac1{2D}\sum_i r_i(0)^2e^{-2\kappa_it}.}
$$

大 $\kappa_i$ 模态快，小模态慢；$\kappa_i=0$ 的标签分量完全不学。核正定时训练残差最终归零，但训练集外误差还需要核的总体谱和目标对齐，不能只从 $\kappa_{\min}\gt 0$ 推出泛化好。

若 $K_t$ 不同时间之间不对易，真实解需要时间有序指数，不能将它写成 $\exp[-\int K_tdt]r(0)$。这也是用“平均一个核”替代完整特征学习动力学的潜在错误。

## 5. 训练集外预测也可以解出

定义训练点与测试点的归一化交叉核

$$
k_0(x,X)_\mu=\frac1D
\nabla f_{\theta_0}(x)^\top\nabla f_{\theta_0}(x_\mu).
$$

由于 $\dot f(x)=-k_0(x,X)^\top r(t)$，积分得到

$$
\boxed{f_t(x)=f_0(x)+k_0(x,X)^\top K_0^+
(I-e^{-K_0t})(y-f_0(X)).}
$$

交叉核对训练 Gram 矩阵零空间正交，所以伪逆表达没有遗漏一项线性增长。若 $f_0=0$、$K_0$ 可逆，长期预测就是无正则化核插值。更一般地，学到的是相对随机初始函数的最小切线范数修正，而非自动得到零均值 GP 后验。

这一公式把两种问题分开：训练集上的模态拟合速度由经验核决定；测试预测是否正确，取决于测试点如何与这些模态耦合。插值只是第一步。

## 6. 怎样检验有限网络是否接近该理论

至少比较三个可观测量：真实网络与初始化切线模型的预测差 $\|f_t-f_t^{\rm lin}\|$、核漂移 $\|K_t-K_0\|_F/\|K_0\|_F$、隐藏特征的任务对齐。总参数范数的相对变化很小，仍可能在某些低维关键方向造成显著函数变化。

核整体变化小，也不自动保证逆核预测变化小；若初始核接近奇异，小特征值方向尤其敏感。应同时报告谱和正则化条件。

![冻结 NTK 模态松弛和训练点外预测。](figures/experiment.svg)

*图 1 · 自制冻结核计算。* 训练点为 $[-\pi,\pi]$ 上 24 点，目标 $\sin x+0.3\cos3x$，RBF 核长度尺度 0.8；采用 $K=\Theta/D$、零初始函数。左图为经验核上几个非零模态的衰减，右图直接用测试公式计算 $t=1,10,1000$ 的预测。这里计算的是精确固定核模型，不把 RBF 核冒称为某个特定架构的 NTK。

复现命令为 `python experiments.py --note 6`。有限网络与切线模型的比较需要共享初始化、训练数据与优化时间；改变学习率而未换算时间，会混入与宽度极限无关的差异。

## 7. 适用边界与相关文献

NTK 的价值是给出一个训练可解析的基准。它解释固定切线特征能够学什么、学多快；不能仅凭该极限描述隐藏单元向教师特征专门化。下一篇研究固定核泛化，Note 08 则改变极限以保留表示演化。

**原始文献。** [Jacot、Gabriel 与 Hongler（2018），Neural Tangent Kernel](https://arxiv.org/abs/1806.07572)提出 NTK 极限；这里的有限训练集恒等式与冻结核解均从损失直接推导，极限成立的架构和初始化条件需与原论文逐项匹配。

---

← [Note 05](/notes/neural_network_theory/gaussian_process/) · [系列目录](/notes/neural_network_theory/) · [Note 07](/notes/neural_network_theory/kernel_learning_curves/) →
