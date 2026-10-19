> **I · 泛化与高维模型 · Note 02**
>
> ← [Note 01](/notes/neural_network_theory/linear_regression/) · [系列目录](/notes/neural_network_theory/) · [Note 03](/notes/neural_network_theory/random_features/) →

**可解性层次：** Gaussian 模型精确期望与高维渐近。

## 1. 先固定任务，再讨论参数数量

当模型刚好能插值训练数据时，测试风险有时出现尖峰；继续增大模型，风险又下降。这个现象不能只用“参数多更复杂”描述，因为插值附近最重要的是数据矩阵的小奇异值。

本篇先研究一个干净的模型族：每个维数 $d$ 下，输入都是标准 Gaussian，教师信号范数固定。这样能精确计算两侧风险。它不是“在同一个无限维任务上依次增加前 $d$ 个坐标”的完整嵌套模型；后者还要计入遗漏信号，才会自然出现最左侧的经典下降段。

## 2. Gaussian 设计与风险的条件分解

设 $X_{\mu i}\sim\mathcal N(0,1)$ 独立，$y=Xw_*+\varepsilon$，$\|w_*\|^2=r^2$，噪声方差 $\sigma^2$。测试输入 $x\sim\mathcal N(0,I_d)$，所以超额风险等于参数误差平方。

最小范数估计为 $\hat w=X^+y$。令 $P=X^+X$，则

$$
\hat w-w_*=-(I-P)w_*+X^+\varepsilon,
$$

$$
\mathbb E_\varepsilon[R_{\rm ex}\mid X]
=\|(I-P)w_*\|^2+\sigma^2\operatorname{Tr}[(X^+)(X^+)^\top].
$$

第一项是教师在不可观测空间中的投影；第二项是反演小奇异值时放大的噪声。两者都由几何决定。

## 3. 欠参数侧：逆 Wishart 给出发散

当 $D\gt d$ 且矩阵满列秩时，$P=I_d$，普通最小二乘没有条件偏差。

$$
\hat w=w_*+(X^\top X)^{-1}X^\top\varepsilon.
$$

若 $D\gt d+1$，逆 Wishart 的一阶矩为

$$
\mathbb E[(X^\top X)^{-1}]=\frac{I_d}{D-d-1},
$$

因而得到有限维的精确期望

$$
\boxed{\mathbb E R_{\rm ex}=\sigma^2\frac{d}{D-d-1}.}
$$

注意这里没有遗漏变量偏差，所以这一侧风险随 $d$ 增加而上升。将它直接画成一条先下降再上升的曲线，会掩盖任务设定发生了变化。

## 4. 过参数侧：未恢复信号与噪声竞争

当 $d\gt D$，

$$
\hat w=X^\top(XX^\top)^{-1}y.
$$

Gaussian 行空间在 $d$ 维中均匀取向，故 $\mathbb E P=(D/d)I_d$。对固定教师，平均不可恢复信号为 $r^2(1-D/d)$。若 $d\gt D+1$，噪声部分由另一侧的逆 Wishart 得到

$$
\boxed{\mathbb E R_{\rm ex}
=r^2\left(1-\frac Dd\right)+\sigma^2\frac{D}{d-D-1}.}
$$

增大 $d$ 降低噪声放大，却增加未恢复教师分量。过参数侧并不保证始终下降，也不保证最终优于所有小模型。取 $\gamma=d/D$ 固定并令规模趋于无穷，

$$
R_{\rm ex}\longrightarrow
\begin{cases}
\sigma^2\gamma/(1-\gamma),&\gamma\lt 1,\\
r^2(1-\gamma^{-1})+\sigma^2/(\gamma-1),&\gamma\gt 1.
\end{cases}
$$

总测试 MSE 还要加不可约噪声 $\sigma^2$。在 $\gamma=1$ 附近，某次有限矩阵的风险通常有限，但它的重尾可能使期望发散。有限次重复的中位数、均值与理论期望不是同一个量。

## 5. Marchenko–Pastur 谱怎样讲同一个故事

令 $H=X^\top X/D$，按 $d$ 个特征值归一化的极限谱为

$$
\rho_\gamma(s)=\left(1-\gamma^{-1}\right)_+\delta_0(s)
+\frac{\sqrt{(s_+-s)(s-s_-)}}{2\pi\gamma s}
\mathbf1_{[s_-,s_+]}(s),
\quad s_\pm=(1\pm\sqrt\gamma)^2.
$$

当 $\gamma\to1$，连续谱下边缘碰到零，逆谱矩变大。过参数侧的零特征值不参与噪声反演，却留下偏差；真正被反演的是非零谱，不能把 $\delta_0/s$ 写成一个普通积分。

ridge 通过滤波 $s/(s+\lambda)$ 避免直接反演小 $s$。若教师方向相对 $X$ 各向同性，风险可写成

$$
\mathbb E R_{\rm ex}\simeq
r^2\int\frac{\lambda^2}{(s+\lambda)^2}\rho_\gamma(s)\,ds
+\sigma^2\gamma\int\frac{s}{(s+\lambda)^2}\rho_\gamma(s)\,ds.
$$

这里 $\lambda\gt 0$，积分连零原子一起定义。记 resolvent $m(\lambda)=\int(s+\lambda)^{-1}\rho_\gamma(s)ds$，两项分别涉及 $-\lambda^2m'(\lambda)$ 与 $m(\lambda)+\lambda m'(\lambda)$。因此“风险取决于谱”可以转化为求 resolvent 及其导数的技术问题。

## 6. 数值图：把发散与有限样本分开看

![Gaussian 最小范数回归的风险与小谱端。](figures/experiment.svg)

*图 1 · 自制 Gaussian 回归实验。* 左图为 $r^2=1,\sigma^2=0.25$ 的高维公式，并叠加脚本中 $D=80$、每个维数 24 次独立设计的风险中位数。误差直接用 $\|\hat w-w_*\|^2$ 测量，不引入有限测试集误差。右图为三个 $\gamma$ 的 MP 连续密度；零原子单独在图例标注，未画成无限高竖线。

复现命令为 `python experiments.py --note 2`。无噪声时，插值尖峰消失，而过参数偏差仍然存在。ridge 改变的是估计器及其整个谱滤波，不能等同于将原风险曲线的尖峰截平。

一个清楚的算例是 $D=100,d=200$：偏差为 $0.5r^2$，噪声为 $100\sigma^2/99$。若 $d=1000$，噪声降至 $100\sigma^2/899$，偏差却增到 $0.9r^2$。这说明“更多参数有益”需要说明信噪比和比较区间。

## 7. 什么能推广，什么还需要证明

小奇异值噪声放大是可推广的机制；Gaussian 逆 Wishart 公式不是任意相关、重尾输入的通式。特征学习也会改变谱与教师对齐，不能把训练后的非线性网络直接替换成独立 Gaussian 矩阵。

**相关文献。** [Hastie、Montanari、Rosset 与 Tibshirani，Surprises in High-Dimensional Ridgeless Least Squares Interpolation](https://arxiv.org/abs/1903.08560)给出更一般的高维风险结果；[Mei 与 Montanari，随机特征回归的精确渐近](https://arxiv.org/abs/1908.05355)研究非线性固定特征。这里的线性推导不冒充后者的全部理论。

---

← [Note 01](/notes/neural_network_theory/linear_regression/) · [系列目录](/notes/neural_network_theory/) · [Note 03](/notes/neural_network_theory/random_features/) →
