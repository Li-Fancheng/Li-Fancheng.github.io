> **I · 泛化与高维模型 · Note 04**
>
> ← [Note 03](/notes/neural_network_theory/random_features/) · [系列目录](/notes/neural_network_theory/) · [Note 05](/notes/neural_network_theory/gaussian_process/) →

**可解性层次：** Hebbian 基线渐近可解；副本为方法框架。

## 1. 让生成规则已知，才能分清“恢复了什么”

真实数据的潜在规则通常不可见。Teacher–student 模型先给定一个教师，再让学生用有限样本学习：研究者知道教师，训练算法只见到样本。这使参数恢复、预测泛化和算法动力学可以分开测量。

它不是一种特定网络架构，而是一类问题设定。改变教师、学生约束、噪声或采样方式，就会改变样本复杂度乃至相变性质。先从一个没有突变阈值的可解例子开始，有助于避免把“相变”当成所有学习曲线的必然结论。

## 2. 感知机教师与方向重叠

取 $x\sim\mathcal N(0,I_d)$，单位教师向量 $u_*$ 给出 $y=\operatorname{sign}(u_*^\top x)$。学生方向 $u$ 也归一化，定义

$$
R=u^\top u_*\in[-1,1].
$$

测试输入的局部场 $h=u^\top x$、$h_*=u_*^\top x$ 是联合 Gaussian，协方差矩阵为 $\begin{pmatrix}1&R\\R&1\end{pmatrix}$。输入在两方向张成的平面里旋转对称；两条分类直线夹角 $\theta=\arccos R$，符号不同的区域占角度比例 $\theta/\pi$。因此

$$
\boxed{\epsilon_g(R)=\mathbb P(hh_*\lt 0)=\frac1\pi\arccos R.}
$$

$R=0$ 对应随机猜测，$R=1$ 对应完全恢复方向。分类误差只依赖方向，与正数倍权重范数无关。这已经说明参数欧氏距离与泛化不必等价。

### 角度公式的概率积分解释

也可把学生局部场写成 $h=Rh_*+\sqrt{1-R^2}z$，其中 $z,h_*$ 为独立标准 Gaussian。给定 $h_*\gt 0$，错误概率为 $\Phi[-Rh_*/\sqrt{1-R^2}]$。利用正负对称性，

$$
\epsilon_g=2\int_0^\infty\frac{e^{-v^2/2}}{\sqrt{2\pi}}
\Phi\!\left(-\frac{Rv}{\sqrt{1-R^2}}\right)dv.
$$

这个积分等于角度表达式，但保留条件概率形式很有意义：分类错误主要来自教师场较小、靠近决策边界的输入。参数方向已经非常接近教师时，远离边界的样本几乎不会改变标签，剩余误差集中在一个狭窄的角楔区域。

由 $R=\cos\theta$，小角度时 $1-R\simeq\theta^2/2$，所以 $\epsilon_g\simeq\sqrt{2(1-R)}/\pi$。重叠误差和分类误差的指数相差一半。讨论学习曲线时，必须先说明测的是哪一个序参量的偏离。

## 3. 一个完整可解的样本学习器

考虑简单 Hebbian 估计，不做损失优化：

$$
v=\frac1D\sum_{\mu=1}^D y_\mu x_\mu,
\qquad u=v/\|v\|.
$$

将输入分解为 $x=zu_*+x_\perp$，$z\sim\mathcal N(0,1)$，便有

$$
\mathbb E[yx]=\mathbb E|z|\,u_*=c\,u_*,
\qquad c=\sqrt{2/\pi}.
$$

每个正交方向的样本均值方差为 $1/D$，共约 $d$ 个这样的方向。令 $D/d\to\alpha\gt 0$，大数集中给出

$$
v^\top u_*\to c,\qquad
\|v_\perp\|^2\to\alpha^{-1},
$$

$$
\boxed{R(\alpha)=\frac{c}{\sqrt{c^2+\alpha^{-1}}},
\qquad\epsilon_g(\alpha)=\frac1\pi\arccos R(\alpha).}
$$

这是特定学习规则的高维渐近曲线，不是 Bayes 最优误差，也不是感知机梯度下降的通式。小 $\alpha$ 时 $R\simeq c\sqrt\alpha$；大 $\alpha$ 时，$\epsilon_g\simeq1/(\pi c\sqrt\alpha)$。幂律在这里来自平均掉 $d$ 个正交噪声分量。

### 有限维涨落的分布可以写得更精确

在教师方向与其正交补的分解中，

$$
v_\parallel=\frac1D\sum_\mu|z_\mu|,
\qquad
v_\perp=\frac1D\sum_\mu\operatorname{sign}(z_\mu)x_{\perp,\mu}.
$$

条件在全部 $z_\mu$ 上，$v_\perp$ 仍为协方差 $I_{d-1}/D$ 的 Gaussian，而且其条件分布不依赖这些 $z_\mu$ 的具体值。因此 $v_\parallel$ 与 $v_\perp$ 独立，且

$$
\|v_\perp\|^2\overset{d}=\frac1D\chi^2_{d-1},
\quad \mathbb E v_\parallel=c,
\quad\operatorname{Var}(v_\parallel)=\frac{1-c^2}{D}.
$$

有限维重叠是这两个随机量的比值，而不是“均值之比”：$R=v_\parallel/\sqrt{v_\parallel^2+\chi^2_{d-1}/D}$。高维公式来自分子和正交范数的同时集中。正交噪声虽然每个坐标只有 $D^{-1/2}$，却在 $d$ 个坐标上累积成 $O(\sqrt{d/D})$ 的总长度，这正是高维统计中不能逐坐标忽略小误差的原因。

Hebbian 学习保留的是所有样本的带标签一阶矩。它没有利用每个样本还施加了一个半空间约束这一更强信息，因此它的曲线是一条可解算法基线，而不是对任务信息极限的求解。比较不同算法需要区分：样本是否含有教师信息、某种估计量是否使用了这些信息，以及优化能否找到相应解。

## 4. 标签噪声怎样改变有效信号

如果标签以概率 $p\lt 1/2$ 独立翻转，令 $y=\zeta\operatorname{sign}(u_*^\top x)$，其中 $\mathbb E\zeta=1-2p$，则 $c$ 被替换为 $c_p=(1-2p)c$，正交噪声方差仍为 $1/D$。

相对干净教师的分类误差仍是 $\arccos R/\pi$；相对新的有噪声测试标签则是

$$
\epsilon_{\rm noisy}=p+(1-2p)\frac{\arccos R}{\pi}.
$$

两种风险的不可约下限不同。若把有噪声测试误差误当方向恢复误差，就会错误判断学生永远无法接近教师。

## 5. 从可解估计走向版本空间与学习相变

对球面学生 $\|w\|^2=d$，可定义零训练错误的版本空间体积

$$
Z(X,y)=\int d\mu(w)\prod_{\mu=1}^D
\Theta\!\left(y_\mu\frac{w^\top x_\mu}{\sqrt d}\right).
$$

其中 $d\mu$ 是球面归一化测度。若希望研究典型样本而非罕见样本主导的平均，应研究 $\mathbb E\log Z$。副本方法形式上使用

$$
\mathbb E\log Z=\lim_{n\to0}\frac{\mathbb E Z^n-1}{n}.
$$

对正整数 $n$，引入副本间重叠 $q_{ab}=w^a\cdot w^b/d$ 和教师重叠 $R_a=w^a\cdot w_*/d$，Gaussian 输入平均只依赖这些协方差，于是高维积分可转为序参量鞍点问题。

这一步解释了自洽方程为什么出现，但**从整数副本继续到 $n=0$ 并不是由该恒等式自动证明的**。还要选择对称性假设、检查稳定性，必要时处理副本对称破缺。本篇不给未推导的鞍点方程贴上“严格精确”的标签。

无噪声、匹配的教师数据保证教师自己始终可行；任意随机标签的可满足性容量是另一个问题。对连续球面、零间隔感知机，经典随机标签阈值 $D/d=2$ 不能直接解释为“有结构教师在两倍样本前完全不能泛化”。离散学生、稀疏先验、非单调激活等设定才可能产生不同的恢复阈值或突变。

### 一个确实产生阈值的对照：随机标签存储

为了看到学习相变的数学来源，可以暂时换成没有教师结构的随机标签任务。连续球面权重满足 $\|w\|^2=d$，要求所有 $D=\alpha d$ 个样本都有间隔 $y_\mu w^\top x_\mu/\sqrt d\geq\kappa$。独立随机标签可吸收到 Gaussian 输入符号里，所以约束等价于全部局部场大于 $\kappa$。

版本空间体积为 $Z=\int d\mu(w)\prod_\mu\Theta(h_\mu-\kappa)$。对整数副本数 $n$，固定副本重叠矩阵 $Q$ 后，各样本的局部场协方差就是 $Q$；权重构型的几何体积在指数阶贡献 $\exp[(d/2)\log\det Q]$。这一步解释了熵项与样本约束项为什么有相同的 $O(d)$ 规模，二者可以竞争产生阈值。

在副本对称 ansatz 中，$Q_{aa}=1$、$Q_{ab}=q$。共同场分解为 $h_a=\sqrt q\,z+\sqrt{1-q}\,t_a$，其中 $z,t_a$ 独立 Gaussian。记 $H(u)=\int_u^\infty Dz$，单个样本的副本约束平均为

$$
\int Dz\left[H\!\left(\frac{\kappa-\sqrt q\,z}{\sqrt{1-q}}\right)\right]^n.
$$

另一方面，$\det Q=(1-q)^{n-1}[1+(n-1)q]$。形式上展开到 $n$ 的一阶，得到每维熵的 RS 鞍点函数（省略与 $q$ 无关的常数）

$$
s_{\rm RS}(q)=\frac12\log(1-q)+\frac{q}{2(1-q)}
+\alpha\int Dz\log H\!\left(
\frac{\kappa-\sqrt q\,z}{\sqrt{1-q}}\right).
$$

几何项偏好较大的可行空间，最后一项是加入样本造成的体积损失。容量边界上，所有可行解挤向彼此，$q\to1$。令 $\epsilon=1-q$，对 $z\lt \kappa$ 使用 Gaussian 尾部 $\log H(u)\simeq-u^2/2$，得到最强发散项

$$
s_{\rm RS}\simeq\frac1{2\epsilon}
\left[1-\alpha\int_{-\infty}^{\kappa}Dz(\kappa-z)^2\right]
+O(\log\epsilon).
$$

临界配平因而给出

$$
\boxed{\alpha_c(\kappa)=
\left[\int_{-\kappa}^\infty Dz(z+\kappa)^2\right]^{-1}.}
$$

在零间隔处，该积分为 $1/2$，得到 $\alpha_c(0)=2$。正间隔使约束更严格，容量降低。这里列出的是经典球面感知机非负间隔问题的 RS 推导；整数副本到零副本的解析延拓和鞍点稳定性仍不是上述代数自动证明的，对负间隔等不同区域更不能直接沿用同一稳定性结论。

### 存储阈值、恢复阈值和算法阈值

这个 $\alpha_c$ 回答随机标签是否仍能被所有约束同时满足。无噪声教师标签始终由教师本身满足，不会在 $\alpha=2$ 突然变得不可满足。若问的是教师方向是否能被唯一恢复，需要研究后验或版本空间沿教师方向的集中；若问某个梯度算法是否能找到恢复态，还要研究动力学盆地与时间尺度。

因此至少有三种不同的阈值概念：解存在的可满足性阈值，给定数据允许达到的统计恢复阈值，以及指定算法实际达到的动力学阈值。简单模型里它们可能重合，也可能分离。Teacher–student 设定的价值，正是让这些问题有独立的定义，而不是把所有陡峭曲线都称作同一种学习相变。

## 6. 数值图：平滑恢复也是重要基线

![Hebbian 学生的方向恢复和分类泛化。](figures/experiment.svg)

*图 1 · 自制 teacher–student 实验。* 输入维数 $d=300$，独立 Gaussian 样本，20 次重复。左图比较高维 $R(\alpha)$ 与有限维测量；右图比较干净标签与 $p=0.2$ 标签翻转时的理论泛化曲线。图中的噪声曲线按有噪声测试标签评分，极限为 $p$。

复现命令为 `python experiments.py --note 4`。当 $\alpha=1$，无噪声理论给出 $R\simeq0.624$、$\epsilon_g\simeq0.286$。有限尺寸比较需要固定 $D/d$；固定 $D$ 而增大 $d$，同时也改变了样本比例。

训练准确率与独立测试准确率也对应不同的局部场分布。相同样本既参与构造 $v$ 又参与打分时，局部场不再服从独立测试的二元 Gaussian 分布，角度公式不能不加修改地用于训练误差。

## 7. 开放问题

最小模型将泛化压缩成一个重叠，但真实表征可能有多个可交换隐藏单元；届时需要重叠矩阵，而非单个 $R$。这在 Note 11 会变成专门化动力学。

**相关文献。** [Gardner（1988），The space of interactions in neural network models](https://doi.org/10.1088/0305-4470/21/1/030)是版本空间方法的经典来源；[Biehl、Riegler 与 Stechert（1995），Learning from noisy data: An exactly solvable model](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.52.R4624)研究有噪声感知机在线学习。这里的批量 Hebbian 推导明确独立于那些算法的具体动力学。

---

← [Note 03](/notes/neural_network_theory/random_features/) · [系列目录](/notes/neural_network_theory/) · [Note 05](/notes/neural_network_theory/gaussian_process/) →
