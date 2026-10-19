> **II · 无限宽与特征学习 · Note 08**
>
> ← [Note 07](/notes/neural_network_theory/kernel_learning_curves/) · [系列目录](/notes/neural_network_theory/) · [Note 09](/notes/neural_network_theory/deep_linear/) →

**可解性层次：** 相应参数化下的平均场极限。

## 1. “无限宽”为什么不唯一

增加隐藏单元数时，输出归一化、初始化尺度、学习率和观察时间必须一起指定。同样是两层网络，一种极限中每个特征几乎不动，另一种极限中单个特征会移动有限距离，经验参数分布发生演化。

所以宽度不是唯一控制量。固定核与特征学习可以分别成为不同缩放下的自洽理论；不能用其中一个极限否认另一个极限的存在。

## 2. 两种参数化的尺度比较

NTK 型参数化采用 $f_N=N^{-1/2}\sum_jh(x;\theta_j)$，在常规学习率和相应初始化下，单粒子梯度通常是 $O(N^{-1/2})$。大量微小贡献可以改变函数，Jacobian 却在合适条件下保持接近初始化。

平均场型参数化采用

$$
f_N(x)=\frac1N\sum_{j=1}^N h(x;\theta_j),
\qquad h(x;\theta)=a\phi(w^\top x),\quad\theta=(a,w).
$$

对总体风险 $\mathcal R_N=\frac12\mathbb E[(f_N(x)-y)^2]$，有

$$
\nabla_{\theta_j}\mathcal R_N=
\frac1N\mathbb E[(f_N-y)\nabla h(x;\theta_j)].
$$

若仍用固定学习率，有限时间粒子会冻结。选择 $N$ 倍粒子学习率，等价于重标时间，得到

$$
\dot\theta_j=-\mathbb E[(f_N-y)\nabla h(x;\theta_j)].
$$

现在速度为 $O(1)$，粒子可以移动有限距离。还可以给 $a,w$ 选择不同迁移率；它们会改变极限度量，不能默默吸收成同一个常数。

## 3. 从粒子系统推导分布方程

定义经验分布 $\rho_N=N^{-1}\sum_j\delta_{\theta_j}$，极限预测写为

$$
f_\rho(x)=\int h(x;\theta)\rho(d\theta).
$$

风险泛函的一阶变分为

$$
\Psi(\theta;\rho)=\frac{\delta\mathcal R}{\delta\rho}(\theta)
=\mathbb E[(f_\rho(x)-y)h(x;\theta)],
$$

差一个与 $\theta$ 无关的常数不会影响速度。对光滑测试函数 $g(\theta)$，

$$
\frac d{dt}\int g\,d\rho_N
=-\int\nabla g\cdot\nabla\Psi(\theta;\rho_N)\,d\rho_N.
$$

若粒子经验分布在相应正则性和矩条件下收敛，弱形式对应

$$
\boxed{\partial_t\rho
=\nabla_\theta\cdot[\rho\nabla_\theta\Psi(\theta;\rho)].}
$$

这是参数空间的连续性方程：概率质量由速度 $v=-\nabla\Psi$ 输运。它常被称作 Wasserstein 梯度流。写成偏微分方程后仍需要边界衰减和可积性，才能合法地分部积分。

## 4. 风险耗散与时变函数核

在无通量边界或充分衰减条件下，

$$
\frac{d\mathcal R}{dt}
=-\int\|\nabla\Psi(\theta;\rho)\|^2\rho(d\theta)\leq0.
$$

这证明该流降低风险，却不单独证明收敛到全局最优。全局结论还需要激活、参数支持、齐次性或初始分布等额外条件。

对预测求导还得到

$$
\partial_t f_\rho(x)
=-\mathbb E_{x',y'}[K_\rho(x,x')(f_\rho(x')-y')],
$$

$$
K_\rho(x,x')=\int
\nabla h(x;\theta)\cdot\nabla h(x';\theta)\,\rho(d\theta).
$$

形式上仍像核动力学，但核由正在移动的分布决定。固定核解不能直接使用，因为 $\rho_t$ 与残差互相耦合。闭合在一个分布上，并不等于闭合在有限个均值或重叠上。

## 5. 一个能看见粒子移动的特征学习实验

取 $x\sim\mathrm{Unif}[-1,1]$，教师为 $\tanh(3x)$，学生为 $N^{-1}\sum_ja_j\tanh(w_jx)$。初始化 $a_j$ 为小随机数、$w_j$ 从窄 Gaussian 取样。全批量训练时

$$
\dot a_j=-\mathbb E[r(x)\tanh(w_jx)],
\quad
\dot w_j=-\mathbb E[r(x)a_jx\operatorname{sech}^2(w_jx)].
$$

固定 $w_j$ 只训练 $a_j$ 是本实验的随机特征对照。它与同时训练 $a,w$ 共享输入、初始化、损失和时间约定；两者的区别有直接的动力学来源。

![平均场粒子演化与固定特征控制组。](figures/experiment.svg)

*图 1 · 自制有限粒子计算。* $N=128$、201 点固定积分网格、步长 0.05，更新使用上面的 $N$ 倍学习率约定。左图比较全参数与固定隐藏层风险；右图用相同粒子连线展示 $(w_j,a_j)$ 的初末位置。它是平均场方程的粒子近似，不是一个无限宽定理的数值证明。

复现命令为 `python experiments.py --note 8`。不同宽度下自然比较的是经验分布，粒子坐标之间没有天然的一一对应。严格对称初始化也可能使流停留在不变子流形，逃逸时间依赖初始扰动及有限尺寸涨落。

## 6. 隐式噪声、熵与理论边界

若给粒子加独立 Brownian 噪声 $\sqrt{2T}\,dB_j$，分布方程增加 $T\Delta\rho$，可以解释为风险加熵泛函的梯度流。真实小批量 SGD 的噪声通常状态依赖且在粒子间相关，不能无条件等同于这里的独立各向同性扩散；Note 14 会详细讨论。

**相关文献。** [Mei、Montanari 与 Nguyen，A Mean Field View of the Landscape of Two-Layers Neural Networks](https://arxiv.org/abs/1804.06561)和 [Chizat 与 Bach，On the Global Convergence of Gradient Descent for Over-parameterized Models using Optimal Transport](https://arxiv.org/abs/1805.09545)给出相关平均场理论。本文只在明确缩放下推导输运方程，不将特定全局收敛条件推广到任意深层网络。

---

← [Note 07](/notes/neural_network_theory/kernel_learning_curves/) · [系列目录](/notes/neural_network_theory/) · [Note 09](/notes/neural_network_theory/deep_linear/) →
