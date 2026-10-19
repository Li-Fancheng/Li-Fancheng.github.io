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

### 条件 Gaussian 与混合分布的区别

由于 $a_j$ 本身是 Gaussian，固定全部隐藏权重以后，有限宽输出向量已经严格为 Gaussian，条件协方差为一个随机经验核。但对隐藏权重再平均，得到的是协方差随机变化的 Gaussian 混合，一般不再严格 Gaussian。无限宽的作用，是让这个随机协方差集中到确定核，而不只是让单点输出“看上去像钟形”。

单个输入的高阶累积量可以显示有限宽修正。写 $b_j=\phi(w_j^\top x/\sqrt d)$，暂去偏置并取 $\sigma_a=1$，则

$$
\mathbb E f_N^2=\mathbb E b^2,\qquad
\operatorname{cum}_4(f_N)
=\frac3N\left[\mathbb E b^4-(\mathbb E b^2)^2\right].
$$

这里使用了 $\mathbb E a^4=3$ 和独立求和时累积量可加。有限四阶矩下，非 Gaussian 性通常以 $1/N$ 衰减；如果四阶矩发散，这条有限宽估计就不成立。对一组输入，还需要对应的联合高阶累积量一起消失，不能仅检查一个输入的峰度。

“固定有限输入集”也有实际含义。宽度无穷时，对任意事先给定的有限样本集都有联合分布极限；若样本数、网络深度或所观察的极端输入随宽度增长，需要对误差作一致控制。有限维分布收敛本身没有给出所有这些同时增长极限。

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

### 递推闭合究竟省略了哪些相关性

在有限宽网络中，不同层的经验核存在涨落，后层神经元也通过前层经验分布发生统计耦合。对固定层数，逐层条件化后，经验协方差的集中使这些涨落在极限中消失，从而只需传播 $q^\ell$。闭合在二阶统计量上，依赖的是预激活 Gaussian 极限及权重独立性，而不是任意非线性随机过程只由协方差决定。

这个区别在训练后尤其重要。全参数训练会使权重依赖样本标签，层间独立结构不再保持初始化时的形式；把训练后的权重方差代入同一前向递推，通常不足以描述其函数分布。NNGP 的直接对象是初始化先验，后验推断或训练动力学需要分别建立。

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

### 二维极坐标积分的展开

将标准二维 Gaussian 写成半径 $r$ 和角度 $\alpha$，$\alpha$ 均匀、$\mathbb E r^2=2$。令 $u=\sqrt q\,r\cos\alpha$、$v=\sqrt{q'}\,r\cos(\alpha-\theta)$，其中 $\theta\in[0,\pi]$ 是两个输入方向的夹角。ReLU 乘积只在两个余弦同时为正的交叠角区非零，区间长度为 $\pi-\theta$。

在该区间直接积分，

$$
\int_{\rm overlap}\cos\alpha\cos(\alpha-\theta)d\alpha
=\frac12[\sin\theta+(\pi-\theta)\cos\theta].
$$

再乘上 $\sqrt{qq'}\mathbb E r^2/(2\pi)$，就得到前面的核公式。它的非线性来自截断后的角区重叠：两个方向越接近，同一批神经元同时激活的比例越大。核并不是对输入距离随意拟合的函数，而是激活区域几何的平均。

### 深度极限的代数退化

对 He 归一化 ReLU，令 $\epsilon_\ell=1-c_\ell$。从 $c=\cos\theta$ 的小角展开，$\sin\theta-\theta\cos\theta\simeq\theta^3/3$，以及 $\theta\simeq\sqrt{2\epsilon}$，得到

$$
\epsilon_{\ell+1}-\epsilon_\ell
\simeq-a\epsilon_\ell^{3/2},\qquad
a=\frac{2\sqrt2}{3\pi}.
$$

把深度视为缓变变量，积分 $d\epsilon/d\ell=-a\epsilon^{3/2}$，给出

$$
\epsilon_\ell^{-1/2}\simeq\epsilon_0^{-1/2}+\frac a2\ell,
\qquad \epsilon_\ell\sim\frac{9\pi^2}{2\ell^2}.
$$

这是足够深处的渐近式，常数不能拿到很浅的层数直接使用。相关性趋于 1 意味着归一化核接近所有输入共享的常量方向，小的区分模态逐渐受压制；它不意味着有限深度时函数严格为常数，而是核的条件数和可分辨方向可能越来越不利。

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

### 方差稳定性与相关稳定性不是一个条件

若等范数输入的对角方差趋向 $q_*$，它满足

$$
q_* =\sigma_b^2+\sigma_w^2\int Dz\,\phi(\sqrt{q_*}z)^2.
$$

对角方差映射的导数，在光滑及可积条件下为 $\sigma_w^2\mathbb E[\phi'^2+\phi\phi'']$。相关系数接近 1 时的导数则是前面的 $\chi=\sigma_w^2\mathbb E\phi'^2$。二者来自不同的扰动：前者同时改变两个变量的方差，后者固定边缘方差而改变共同协方差。

这里可用 Gaussian 分部积分，也称 Price 恒等式：对固定方差的二元 Gaussian，$\partial_C\mathbb E[\phi(u)\phi(v)]=\mathbb E[\phi'(u)\phi'(v)]$。将 $C=q_*c$ 代入，归一化时的 $q_*$ 抵消，才得到 $\chi$。因此“方差不爆炸”并不足以保证输入之间的差异能跨很多层保留。

文献中前向传播的有序/混沌区域，指输入差异随层数的放大或收缩；这里层数扮演迭代方向。Note 16 的循环网络混沌则是同一个网络状态随真实动力学时间的敏感依赖。两种分析共享平均场思想，但不能将临界线、Lyapunov 指数或记忆时间直接互换。

### GP 后验还保留不确定度

Gaussian 观测模型下，测试函数与训练标签组成联合 Gaussian。对块协方差矩阵作条件化，就得到前面的后验均值和 Schur 补方差。均值提供一个核回归预测，方差则衡量在该先验下数据尚未约束的函数方向。

若测试点与所有训练点的交叉核都很小，后验均值接近先验均值，但方差接近先验方差。这与“确定地预测零”是不同的判断。相反，先验若错误地把两个不同任务区域看得非常相似，后验也可能低估真正的不确定性；Bayesian 计算的严格性始终以所选择的先验和观测模型为条件。

## 6. 数值图与实验

![ReLU Gaussian process 核与逐层相关性递推。](figures/experiment.svg)

*图 1 · 自制核计算。* 左图比较二维 Gaussian Monte Carlo 积分与上述 ReLU 解析式，方差均为 1，使用 200000 个基础随机样本。右图取零偏置、$\sigma_w^2=2$，对 $c_0=-0.5,0,0.5,0.9$ 迭代 60 层。曲线表示随机初始化的先验相关，不是训练损失。

复现命令为 `python experiments.py --note 5`。$c=0$ 时的积分值为 $1/(2\pi)$。有限宽网络的经验核还随初始化波动，平均核的收敛与单个网络的函数分布收敛需要分别考察。

## 7. 适用边界与相关文献

GP 描述的是随机函数的有限维分布。非 Gaussian、相关或重尾初始化，残差结构、卷积权重共享以及深宽同时增长，都可能需要不同论证。先验相似也不能保证训练后的分布相同。

**相关文献。** [Lee 等，Deep Neural Networks as Gaussian Processes](https://arxiv.org/abs/1711.00165)建立深层 NNGP 递推；[Schoenholz 等，Deep Information Propagation](https://arxiv.org/abs/1611.01232)研究方差和相关性传播。本文将核积分、极限条件和训练解释分别处理。

---

← [Note 04](/notes/neural_network_theory/teacher_student/) · [系列目录](/notes/neural_network_theory/) · [Note 06](/notes/neural_network_theory/ntk/) →
