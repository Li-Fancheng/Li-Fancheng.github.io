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

## 4. 有效维数与任务对齐

定义两个相关但不同的量

$$
d_{\rm eff}(\lambda)=\sum_i\frac{\eta_i}{\eta_i+\lambda},
\qquad d_2(\lambda)=\sum_i\left(\frac{\eta_i}{\eta_i+\lambda}\right)^2.
$$

本代理的方差恰为 $\sigma^2d_2/D$，不是无条件等于 $\sigma^2d_{\rm eff}/D$。对于常见规则谱，两者常有相同标度，但常数与精确公式不同。

一个两模态例子很说明问题。令 $\eta_1=1,\eta_2=0.01$，取 $\lambda=0.1$。目标若全在第一模态，信号恢复比例约 $0.909$；若全在第二模态，只恢复 $0.091$。两个任务使用相同核、相同数据预算，却具有完全不同的偏差。

因此“低频先学”只有在核本征函数恰好按某个频率指标排序时才成立。在不平移不变的分布或核中，谱模态可能是局域、几何或语义结构，不能靠本征值大小直接命名空间频率。

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

## 6. 数值图与复现

![核谱滤波、目标对齐和学习曲线。](figures/experiment.svg)

*图 1 · 自制高斯序列代理。* 使用 10000 个模态、$\eta_i=i^{-2}$，对比 $a_i^2\propto i^{-2}$ 与 $i^{-4}$，均归一化为单位目标能量。左图给出两种目标的谱能量与核谱；右图对每个 $D$ 在对数网格上选择最小解析风险，噪声 $\sigma=0.3$。这是在已知教师下计算的最优正则化基线，实际实验必须用验证集选择超参数。

复现命令为 `python experiments.py --note 7`。核谱和目标系数是影响学习曲线的两个独立因素。有限模态截断会在极大数据处引入新的区间，因此单一指数只适用于远离截断的尺度窗口。

## 7. 相关文献

**原始文献。** [Bordelon、Canatar 与 Pehlevan（2020），Spectrum Dependent Learning Curves in Kernel Regression and Wide Neural Networks](https://proceedings.mlr.press/v119/bordelon20a.html)系统研究谱与目标对齐。本文的独立高斯序列实验用于给出透明的逐模推导；随机采样的自洽描述则保留其额外假设。

---

← [Note 06](/notes/neural_network_theory/ntk/) · [系列目录](/notes/neural_network_theory/) · [Note 08](/notes/neural_network_theory/lazy_and_mean_field/) →
