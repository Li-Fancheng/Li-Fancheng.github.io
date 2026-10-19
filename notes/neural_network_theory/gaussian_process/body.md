> **II · 无限宽与特征学习 · Note 05**
>
> ← [Note 04](/notes/neural_network_theory/teacher_student/) · [系列目录](/notes/neural_network_theory/) · [Note 06](/notes/neural_network_theory/ntk/) →

**可解性层次：** 固定有限输入集上的无限宽极限。

## 1. 无限多个随机神经元，为什么留下一个随机函数

随机初始化网络对不同输入给出相关的随机输出。网络很宽时，这些输出往往接近联合 Gaussian；输入之间的相关性并没有消失，而是被协方差核保存下来。无限宽极限因此给出的是一套函数先验，不是一个确定函数，也还没有描述训练。

Gaussian process（GP）的定义是：对任意固定有限输入集，函数值组成多元 Gaussian。每个输入上的单点直方图近似正态，尚不足以建立 GP 结论；还必须处理不同输入的联合分布。

## 2. 从两层网络开始的中心极限定理

设

$$
f_N(x)=\frac{\sigma_a}{\sqrt N}\sum_{j=1}^N
a_j\phi(w_j^\top x/\sqrt d)+b,
$$

其中 $a_j\sim\mathcal N(0,1)$、$w_j\sim\mathcal N(0,\sigma_w^2I_d)$、$b\sim\mathcal N(0,\sigma_b^2)$，各组独立。对固定 $x_1,\ldots,x_D$，每个神经元贡献一个 $D$ 维随机向量。若其二阶矩有限，多元中心极限定理给出

$$
(f_N(x_1),\ldots,f_N(x_D))
\Rightarrow\mathcal N(0,K),
$$

$$
K(x,x')=\sigma_b^2+\sigma_a^2
\mathbb E_w[\phi(w^\top x/\sqrt d)\phi(w^\top x'/\sqrt d)].
$$

输出均值为零来自读出权重的零均值，并不要求 ReLU 特征自身零均值。独立性、有限矩与 $N^{-1/2}$ 归一化都是结论的一部分。

## 3. 多层网络：条件 Gaussian 与协方差递推

对于全连接网络，令第 $\ell+1$ 层预激活为

$$
h_i^{\ell+1}(x)=\frac{\sigma_w}{\sqrt{n_\ell}}
\sum_jW_{ij}^{\ell+1}\phi(h_j^\ell(x))+b_i^{\ell+1},
$$

其中 $W$ 为单位方差独立 Gaussian，偏置方差为 $\sigma_b^2$。给定上一层激活，下一层预激活严格是条件 Gaussian，其条件协方差为一个经验平均。宽度趋于无穷时，该平均集中到确定值，从而得到递推

$$
\boxed{q^{\ell+1}(x,x')=\sigma_b^2+
\sigma_w^2\mathbb E_{(u,v)\sim\mathcal N(0,Q^\ell)}
[\phi(u)\phi(v)],}
$$

$$
Q^\ell=\begin{pmatrix}
q^\ell(x,x)&q^\ell(x,x')\\
q^\ell(x',x)&q^\ell(x',x')
\end{pmatrix}.
$$

这一步的要点是先条件化再取集中极限。直接说“每层都是许多独立量相加”会忽略不同神经元共享前层随机性的事实。固定深度、固定有限输入集是最干净的极限设置；深度、样本数也增长时需要额外控制。

## 4. ReLU 核的显式计算

若 $u,v$ 方差分别为 $q,q'$，相关系数 $c=\cos\theta$，对二维 Gaussian 的极坐标积分给出

$$
\mathbb E[\operatorname{ReLU}(u)\operatorname{ReLU}(v)]
=\frac{\sqrt{qq'}}{2\pi}
[\sin\theta+(\pi-\theta)\cos\theta].
$$

端点结果为：$c=1$ 时等于 $\sqrt{qq'}/2$；$c=-1$ 时为零；$c=0$ 时为 $\sqrt{qq'}/(2\pi)$。最后一个数不是零，因为 ReLU 有正均值。

取零偏置和 He 方差 $\sigma_w^2=2$，若两输入方差相等，方差逐层保持，相关系数满足

$$
c_{\ell+1}=F(c_\ell)
=\frac{\sqrt{1-c_\ell^2}+(\pi-\arccos c_\ell)c_\ell}{\pi}.
$$

$F(1)=1$，不同输入的相关性逐渐趋近 1，深层先验可能越来越难区分输入。这里 $F'(1)=1$，不能套用严格收缩映射的指数收敛。展开得到

$$
F(c)-c\simeq\frac{2\sqrt2}{3\pi}(1-c)^{3/2},
$$

所以连续深度近似给出 $1-c_\ell\sim\ell^{-2}$。这是本例的边缘情况，并不是所有激活的普适深度指数。

## 5. 信号传播与训练之间隔着什么

对一般激活，在方差固定点 $q_*$ 附近，小相关差异的线性放大率常写为

$$
\chi=\sigma_w^2\mathbb E_{z\sim\mathcal N(0,1)}
[\phi'(\sqrt{q_*}z)^2].
$$

在适用的递推条件下，$\chi\lt 1$ 表示相关差异收缩，$\chi\gt 1$ 表示其在该点附近增长。这个前向相关分析也与梯度尺度有关，但它没有单独证明优化成功或泛化优良；还需看任务、读出核谱以及训练造成的分布变化。

若对 NNGP 先验做 Gaussian 噪声回归，后验均值为

$$
m(x)=K(x,X)[K(X,X)+\sigma_n^2I]^{-1}y,
$$

后验方差为 $K(x,x)-K(x,X)[K+\sigma_n^2I]^{-1}K(X,x)$。这里的矩阵正则化量是观测噪声方差；若损失采用 $1/D$ 归一化，与 ridge 参数比较时必须带上 $D$。GP 后验也不等于训练所有网络参数的梯度下降结果；后者的固定核极限通常由下一篇的 NTK 描述。

## 6. 数值图与实验

![ReLU Gaussian process 核与逐层相关性递推。](figures/experiment.svg)

*图 1 · 自制核计算。* 左图比较二维 Gaussian Monte Carlo 积分与上述 ReLU 解析式，方差均为 1，使用 200000 个基础随机样本。右图取零偏置、$\sigma_w^2=2$，对 $c_0=-0.5,0,0.5,0.9$ 迭代 60 层。曲线表示随机初始化的先验相关，不是训练损失。

复现命令为 `python experiments.py --note 5`。$c=0$ 时的积分值为 $1/(2\pi)$。有限宽网络的经验核还随初始化波动，平均核的收敛与单个网络的函数分布收敛需要分别考察。

## 7. 适用边界与相关文献

GP 描述的是随机函数的有限维分布。非 Gaussian、相关或重尾初始化，残差结构、卷积权重共享以及深宽同时增长，都可能需要不同论证。先验相似也不能保证训练后的分布相同。

**相关文献。** [Lee 等，Deep Neural Networks as Gaussian Processes](https://arxiv.org/abs/1711.00165)建立深层 NNGP 递推；[Schoenholz 等，Deep Information Propagation](https://arxiv.org/abs/1611.01232)研究方差和相关性传播。本文将核积分、极限条件和训练解释分别处理。

---

← [Note 04](/notes/neural_network_theory/teacher_student/) · [系列目录](/notes/neural_network_theory/) · [Note 06](/notes/neural_network_theory/ntk/) →
