> **II · 无限宽与特征学习 · Note 07**
>
> ← [Note 06](/notes/neural_network_theory/ntk/) · [系列目录](/notes/neural_network_theory/) · [Note 08](/notes/neural_network_theory/lazy_and_mean_field/) →

**可解性层次：** 总体代理精确；随机设计学习曲线需附加假设。

## 1. 一样的核，为什么有的目标容易、有的很难

核的特征值描述它偏好哪些函数方向；目标函数在这些方向上的系数描述任务需要哪些方向。仅给出核谱，而不说明目标对齐，通常无法决定泛化曲线。

本篇先建立总体核谱与有限样本 Gram 矩阵的区别，再用一个可精确计算的高斯序列代理推导学习曲线。随机设计回归中额外出现的耦合单独说明，避免把“总体协方差代替经验协方差”误写成普遍精确解。

## 2. Mercer 展开与目标能量

输入分布为 $\rho$，总体核积分算子为

$$
(Tg)(x)=\int k(x,x')g(x')\,d\rho(x').
$$

对适当平方可积正定核，取 $L^2(\rho)$ 正交本征函数 $\psi_i$，则

$$
k(x,x')=\sum_i\eta_i\psi_i(x)\psi_i(x'),
\qquad f_*(x)=\sum_i a_i\psi_i(x),
$$

其中 $\int\psi_i\psi_jd\rho=\delta_{ij}$、$\eta_i\geq0$。$\eta_i$ 是总体算子的谱；训练矩阵 $k(X,X)/D$ 的特征值是随机经验谱，有限 $D$ 时二者不等同。

若目标在核零空间有能量，该部分无法由该核表示，形成不可约近似误差。对于正谱方向，RKHS 范数为 $\sum_i b_i^2/\eta_i$，所以小特征值系数受到更强的正则化。

### RKHS 与总体平方误差使用不同的度量

在正谱方向上写 $f=\sum_i b_i\psi_i$。由于 $\psi_i$ 在输入分布下正交，函数误差是 $\|f-f_*\|_{L^2(\rho)}^2=\sum_i(b_i-a_i)^2$；但 RKHS 范数是 $\sum_i b_i^2/\eta_i$。同一个坐标在两种度量下具有不同权重，正则化的谱选择正是从这个差异产生。

若目标有一部分 $f_\perp$ 正交于所有正谱本征函数，它会给总体误差增加 $\|f_\perp\|_{L^2}^2$。这与标签观测噪声不同：它可能是确定而有结构的函数，只是当前核无法表达。将不可表示目标都并入一个常数噪声项，只有在另有独立性或统计等效论证时才合理。

核回归的 representer 结构也可从几何上得到。将函数分成训练核截面的线性组合与其 RKHS 正交补，$f=f_\parallel+f_\perp^{\mathcal H}$。由再生性质，正交补在每个训练点的值都是零；它不会降低训练误差，却增加范数。因此带正则化的最优解位于 $\operatorname{span}\{k(x_\mu,\cdot)\}$。这只是训练解所在的函数子空间，没有保证该子空间已经覆盖所有总体重要模态。

## 3. 可解代理：直接观察带噪谱系数

考虑 Gaussian 序列实验

$$
z_i=a_i+\frac\sigma{\sqrt D}\xi_i,
\qquad \xi_i\sim\mathcal N(0,1).
$$

最小化每模的 $(b_i-z_i)^2/2+\lambda b_i^2/(2\eta_i)$，得到

$$
\hat a_i=\frac{\eta_i}{\eta_i+\lambda}z_i.
$$

这是一个有明确定义的代理观测问题。在有限个正交经验特征且经验 Gram 恰好为单位阵时，它也对应普通回归；无限多个独立白噪系数本身不在 $\ell_2$ 中，需要滤波后的方差和收敛。

平均超额风险为

$$
\boxed{E(D,\lambda)=\sum_i
\underbrace{\left(\frac\lambda{\eta_i+\lambda}\right)^2a_i^2}_{\text{偏差}}
+\underbrace{\frac{\sigma^2}{D}
\left(\frac{\eta_i}{\eta_i+\lambda}\right)^2}_{\text{方差}}.}
$$

若对总体损失做梯度流而不加正则化，同样有 $\hat a_i(t)=(1-e^{-\eta_it})a_i$。所以核谱既控制优化速度，也控制正则化下可恢复的模态；但改变 $t$ 与改变样本数 $D$ 是不同操作。

### 从真实随机采样到代理，中间缺少的矩阵

暂时截断到 $M$ 个总体本征函数，令 $Z_{\mu i}=\psi_i(x_\mu)$、$E=\operatorname{diag}(\eta_i)$。真实样本给出的是 $y=Za+\varepsilon$，不是直接给每个 $a_i$ 加一个独立噪声。核回归等价于

$$
\min_b\frac1{2D}\|Zb-y\|^2+\frac\lambda2 b^\top E^{-1}b.
$$

记 $S=Z^\top Z/D$、$A=S+\lambda E^{-1}$，则

$$
\hat a=A^{-1}Z^\top y/D,
\qquad
\hat a-a=-\lambda A^{-1}E^{-1}a+A^{-1}Z^\top\varepsilon/D.
$$

给定设计矩阵、对标签噪声平均后，有限维模型的总体误差恰为

$$
\boxed{\mathbb E_\varepsilon[E_g\mid Z]
=\|\lambda A^{-1}E^{-1}a\|^2
+\frac{\sigma^2}{D}\operatorname{Tr}(A^{-1}SA^{-1}).}
$$

右边的条件平均记号表示先固定 $Z$。关键障碍是 $S$ 随采样波动，而且一般与 $E$ 不对易。总体正交性只给出 $\mathbb ES=I$，并没有给出 $S=I$，也没有给出 $\mathbb E A^{-1}=(\mathbb EA)^{-1}$。

如果恰好 $S=I$，各模态才解耦成 $\hat a_i=\eta_i z_i/(\eta_i+\lambda)$，上一节的序列风险随之恢复。这明确指出代理省略了什么：经验模态之间的混合、有限秩约束，以及矩阵反演对这些涨落的放大。即使每个 $S_{ij}$ 的误差只有 $D^{-1/2}$，当模态数与 $D$ 同阶时，整个算子的谱误差也未必小。

对于被截断的无限尾部，$y=Za+f_{\rm tail}(X)+\varepsilon$。尾部与已保留本征函数在总体上正交，但有限样本的相关不为零；它会像混叠项进入估计。这也是从有限 $M$ 公式推广到无限核时，需要控制截断而不能直接丢掉尾部的原因。

## 4. 有效维数与任务对齐

定义两个相关但不同的量

$$
d_{\rm eff}(\lambda)=\sum_i\frac{\eta_i}{\eta_i+\lambda},
\qquad d_2(\lambda)=\sum_i\left(\frac{\eta_i}{\eta_i+\lambda}\right)^2.
$$

本代理的方差恰为 $\sigma^2d_2/D$，不是无条件等于 $\sigma^2d_{\rm eff}/D$。对于常见规则谱，两者常有相同标度，但常数与精确公式不同。

一个两模态例子很说明问题。令 $\eta_1=1,\eta_2=0.01$，取 $\lambda=0.1$。目标若全在第一模态，信号恢复比例约 $0.909$；若全在第二模态，只恢复 $0.091$。两个任务使用相同核、相同数据预算，却具有完全不同的偏差。

因此“低频先学”只有在核本征函数恰好按某个频率指标排序时才成立。在不平移不变的分布或核中，谱模态可能是局域、几何或语义结构，不能靠本征值大小直接命名空间频率。

### 两种有效维数的物理含义

$d_{\rm eff}$ 是平滑算子 $T(T+\lambda I)^{-1}$ 的迹，可以理解成模型保留的连续自由度：每个模态按恢复比例贡献 0 到 1。$d_2$ 则是该平滑算子平方的迹，表示独立噪声经过滤波后的总能量。前者计数响应，后者计数涨落，二者结构相近而不相同。

若 $\eta_i\sim c i^{-p}$、$p\gt 1$，将求和变成积分并以 $i=(c/\lambda)^{1/p}u$ 重标，可得

$$
d_{\rm eff}\sim(c/\lambda)^{1/p}
\int_0^\infty\frac{du}{1+u^p},
\quad
d_2\sim(c/\lambda)^{1/p}
\int_0^\infty\frac{du}{(1+u^p)^2}.
$$

第一个积分等于 $\pi/[p\sin(\pi/p)]$，第二个是前者的 $(1-1/p)$ 倍。所以在这个规则谱里 $d_2/d_{\rm eff}\to1-1/p$。它们有相同幂次，却连领先常数都不同；用“有效参数数”一个词代替所有谱迹，会把这种区别抹平。

## 5. 随机设计学习曲线为什么多一个自洽尺度

真实训练样本只观察 $f_*(x_\mu)+\varepsilon_\mu$，不同谱系数通过随机采样混合。对某类随机特征设计或统计物理近似，可使用等效尺度 $\kappa$ 与涨落因子 $\gamma$：

$$
\kappa=\Lambda+\sum_i\frac{\kappa\eta_i}{\kappa+D\eta_i},
\qquad
\gamma=\sum_i\frac{D\eta_i^2}{(\kappa+D\eta_i)^2},
$$

$$
E\simeq\frac1{1-\gamma}
\sum_i\frac{\kappa^2a_i^2+\sigma^2D\eta_i^2}
{(\kappa+D\eta_i)^2}.
$$

这里使用**未除以 $D$ 的平方损失**加 $\Lambda\|f\|_{\mathcal H}^2/2$；与前一节参数比较时，$\Lambda=D\lambda$。上述方程不是对所有 Mercer 本征函数分布都无条件精确的定理。它展示经验采样涨落如何改变有效正则化，$1/(1-\gamma)$ 如何放大误差；使用时需检查设计假设、极限和解的稳定分支。

本系列的数据幂律计算将明确采用上一节可解代理，不把这套自洽方程的适用性扩展到未经论证的数据集。

### 自洽方程的尺度结构

将 $\kappa$ 方程除以 $\kappa$，可写成

$$
1=\frac\Lambda\kappa+\sum_i\frac{\eta_i}{\kappa+D\eta_i}.
$$

对于 $\Lambda\gt 0$，右边随 $\kappa$ 单调降低，因此正解唯一。$\kappa/D$ 扮演的并不是手动设置的 $\lambda$，而是将有限采样反馈吸收后的等效谱截止。在样本数有限时，即使显式正则化趋零，等效截止仍可能非零。

涨落因子 $\gamma$ 包含模态响应平方。可以把 $1/(1-\gamma)$ 形式地看成反复反馈的几何级数 $1+\gamma+\gamma^2+\cdots$：一次样本扰动改变估计，估计改变又通过同一份设计反馈到其他方向。这个图像有助于理解误差放大，但不能替代严格的 cavity、replica 或随机矩阵推导；对具体设计，正是这些方法需要证明反馈项如何集中。

### 等谱特例把核曲线接回双下降

取 $M$ 个模态、全部 $\eta_i=1$、目标能量 $\sum_i a_i^2=r^2$，考察 $\Lambda\downarrow0$ 的分支。当 $M\gt D$ 时，自洽方程有 $\kappa=M-D$，且 $\gamma=D/M$。代入得到

$$
E_g=\left(1-\frac DM\right)r^2
+\sigma^2\frac{D}{M-D}.
$$

当 $D\gt M$ 时 $\kappa\to0$、$\gamma=M/D$，于是 $E_g=\sigma^2 M/(D-M)$。这正好恢复 Note 02 的高维各向同性风险，有限样本逆 Wishart 的 $-1$ 修正则不在该渐近方程里。

这一检查很重要：独立序列代理若在 $\lambda=0$ 时只给出 $\sigma^2M/D$，它没有插值分母；自洽采样理论则恢复 $D-M$ 或 $M-D$ 的病态。这两套曲线回答不同观测模型的问题，不能在推导数据指数时随意交换。

### 谱与本征函数统计都属于问题定义

在某些高维随机设计极限中，谱和目标能量分布足以决定风险的领先项。但对一般输入，本征函数可能高度局域化，少数样本才落入关键区域；即使本征值完全相同，采样覆盖也可能不同。两个核积分算子同谱，不等于其随机采样 Gram 矩阵具有相同统计。

因此将谱理论应用到实际网络时，除了测 $\eta_i$ 和 $a_i$，还需要检查本征函数的集中性、训练/测试分布变化以及经验模态的稳定性。谱偏置给出了清晰语言，但只有把采样几何纳入，才可能成为完整的泛化解释。

## 6. 数值图与复现

![核谱滤波、目标对齐和学习曲线。](figures/experiment.svg)

*图 1 · 自制高斯序列代理。* 使用 10000 个模态、$\eta_i=i^{-2}$，对比 $a_i^2\propto i^{-2}$ 与 $i^{-4}$，均归一化为单位目标能量。左图给出两种目标的谱能量与核谱；右图对每个 $D$ 在对数网格上选择最小解析风险，噪声 $\sigma=0.3$。这是在已知教师下计算的最优正则化基线，实际实验必须用验证集选择超参数。

复现命令为 `python experiments.py --note 7`。核谱和目标系数是影响学习曲线的两个独立因素。有限模态截断会在极大数据处引入新的区间，因此单一指数只适用于远离截断的尺度窗口。

## 7. 相关文献

**原始文献。** [Bordelon、Canatar 与 Pehlevan（2020），Spectrum Dependent Learning Curves in Kernel Regression and Wide Neural Networks](https://proceedings.mlr.press/v119/bordelon20a.html)系统研究谱与目标对齐。本文的独立高斯序列实验用于给出透明的逐模推导；随机采样的自洽描述则保留其额外假设。

---

← [Note 06](/notes/neural_network_theory/ntk/) · [系列目录](/notes/neural_network_theory/) · [Note 08](/notes/neural_network_theory/lazy_and_mean_field/) →
