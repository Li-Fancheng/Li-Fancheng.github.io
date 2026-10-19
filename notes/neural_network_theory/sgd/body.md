> **IV · 记忆与集体动力学 · Note 14**
>
> ← [Note 13](/notes/neural_network_theory/grokking/) · [系列目录](/notes/neural_network_theory/) · [Note 15](/notes/neural_network_theory/hopfield/) →

**可解性层次：** 局部线性噪声模型精确；SGD 扩散为近似。

## 1. 小批量训练里的随机性来自哪里

同一个参数点，抽取不同批次会得到不同梯度。靠近稳定极小点时，可以把平均梯度线性化、把批次噪声暂时冻结，得到 Ornstein–Uhlenbeck（OU）过程。这个近似能精确解释局部协方差，却不能把整段神经网络训练自动变成热平衡。

尤其要分清三个量：学习率、批次梯度协方差，以及真正满足涨落耗散关系的温度。它们并非同义词。

## 2. 从离散更新到局部扩散

设总体经验损失 $\mathcal L=D^{-1}\sum_\mu\ell_\mu$，单样本梯度 $g_\mu$ 的均值为 $g$，定义

$$
C(w)=\frac1D\sum_\mu(g_\mu-g)(g_\mu-g)^\top.
$$

大小为 $B$、有放回独立采样的批次均值满足 $\operatorname{Cov}(g_B)=C/B$。无放回采样则为 $(D-B)C/[B(D-1)]$。SGD 更新为

$$
w_{k+1}=w_k-\eta g(w_k)+\eta\xi_k,
\qquad \operatorname{Cov}(\xi_k\mid w_k)=C(w_k)/B.
$$

取时间 $t=k\eta$，匹配每单位时间的一、二阶矩，得到小步长扩散近似

$$
dw=-\nabla\mathcal L(w)dt+
\sqrt{\frac\eta B C(w)}\,dW_t.
$$

它省略高阶批次累积量、离散误差与批次间相关。Brownian 近似要求相应小步长与局部统计条件，不能只因算法名字是 SGD 就直接使用。

## 3. 二次极小点附近的精确 OU 模型

令 $u=w-w_*$，$\mathcal L\simeq\mathcal L_*+u^\top Hu/2$，$H\succ0$。冻结扩散矩阵 $A=\eta C(w_*)/(2B)$，得到

$$
du=-Hu\,dt+\sqrt{2A}\,dW_t.
$$

均值满足 $\dot m=-Hm$。协方差 $S=\mathbb E[(u-m)(u-m)^\top]$ 用 Itô 公式得到

$$
\dot S=-HS-SH+2A.
$$

长期驻态满足连续 Lyapunov 方程

$$
\boxed{HS+SH=2A,\qquad
S=2\int_0^\infty e^{-Ht}Ae^{-Ht}dt.}
$$

在 $H$ 的本征基里，$S_{ij}=2A_{ij}/(h_i+h_j)$。若两矩阵同时对角化，则 $S_{ii}=A_{ii}/h_i$：软方向涨落大，还要看该方向的噪声注入有多强。仅凭 Hessian 不能决定驻态占据。

## 4. 什么时候可以说“有效温度”

如果 $A=TI$，驻态为 $p_*(u)\propto\exp[-u^\top Hu/(2T)]$，$S=TH^{-1}$，满足详细平衡。一般各向异性 $A$ 产生的 Gaussian 驻态仍然存在，但未必是关于原损失的单温度 Gibbs 分布。

Fokker–Planck 概率流为

$$
j(u)=(-H+A S^{-1})u\,p_*(u).
$$

详细平衡要求 $HS=A$。对称 $H,A$ 不对易时，这通常不成立，驻态可以持续有旋转概率流；“分布不再随时间变”不等于“没有非平衡流”。即使 $H,A$ 对易但 $A\neq TI$，也往往只有逐方向温度，不能用一个标量温度乘原损失。

若 $H$ 有零模、而 $A$ 沿零模仍注入噪声，该方向是扩散，没有可归一化的 OU 驻态。若有负曲率，则局部线性模型不稳定。将这些方向硬塞进 $S=A/H$ 没有意义。

## 5. 离散时间的修正

对于冻结噪声的真实离散线性更新

$$
u_{k+1}=(I-\eta H)u_k+\eta\xi_k,
$$

驻态满足离散 Lyapunov 方程

$$
S=(I-\eta H)S(I-\eta H)+\eta^2 C/B.
$$

若同时对角化，

$$
\boxed{S_{ii}=\frac{\eta C_{ii}}{B h_i(2-\eta h_i)}.}
$$

在 $\eta h_i\ll1$ 时恢复连续近似；接近 $\eta h_i=2$ 的稳定边界时涨落显著增强。把 OU 近似用于大步长训练，可能恰好错过这种离散效应。

## 6. 数值图与实验

![二维 OU 驻态的协方差与各向异性涨落。](figures/experiment.svg)

*图 1 · 自制 OU 模型。* 使用 $H=\operatorname{diag}(1,4)$、$A=\begin{pmatrix}1&0.8\\0.8&2\end{pmatrix}$，理论 $S=\begin{pmatrix}1&0.32\\0.32&0.5\end{pmatrix}$。左图为 5000 条独立 Euler–Maruyama 轨迹在 $t=7$ 的终点及理论两倍标准差椭圆；右图比较协方差。步长 0.002，故存在小离散误差和有限采样误差。

复现命令为 `python experiments.py --note 14`。本次相对协方差误差约为 $0.28\%$。当 $A$ 为标量单位阵时恢复单温度 Gibbs 关系；当前非对易的 $H,A$ 则需要同时考察协方差和驻态概率流。

## 7. 理论边界与相关文献

真实 SGD 的 $C(w)$ 会变，噪声可能低秩、非 Gaussian，与曲率耦合；在插值且每个样本梯度都为零的点，批次噪声甚至消失。学习率与批次大小之比只能在明确近似下充当噪声尺度，不能无条件当作热力学温度。

**相关文献。** [Mandt、Hoffman 与 Blei，Stochastic Gradient Descent as Approximate Bayesian Inference](https://arxiv.org/abs/1704.04289)使用局部近似研究 SGD 的分布。本文的 OU 模型在所给条件下精确，SGD 到 OU 的替换则保留近似身份。

---

← [Note 13](/notes/neural_network_theory/grokking/) · [系列目录](/notes/neural_network_theory/) · [Note 15](/notes/neural_network_theory/hopfield/) →
