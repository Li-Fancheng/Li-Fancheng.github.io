> **I · 泛化与高维模型 · Note 01**
>
> [系列目录](/notes/neural_network_theory/) · [Note 02](/notes/neural_network_theory/kernel_methods/) →

**可解性层次：** 有限设计精确求解；Gaussian 随机设计的精确期望与高维渐近。

## 前言：先解一份数据，再研究高维数据的随机性

线性回归与 double descent 实际上在研究同一个反演问题的两个层次。固定一份设计矩阵，可以精确追踪每个方向怎样恢复信号、怎样吸收噪声；再让设计矩阵随机，便能研究这些方向的统计分布，以及插值附近为什么出现风险峰。将两层计算连起来，小特征值就同时获得了动力学与统计含义：它既使恢复缓慢，也使长期噪声反演危险。

这一篇因此先保留有限维答案，再进入 Gaussian 高维模型。前半部分把训练、参数与测试误差分开，后半部分计算随机投影和逆 Wishart 如何进入风险。早停、ridge、零空间与最小范数不是旁支，它们共同决定给定数据允许观察什么，以及算法怎样使用这些信息。

我关心的不是用一条曲线概括“参数多好不好”，而是明确曲线中哪些变化来自信号覆盖，哪些来自噪声放大，哪些又只是训练时间与归一化不同。只有把任务、估计规则与平均方式固定，插值两侧的比较才有清楚含义。

## 1. 固定设计下的模型、风险与可辨识方向 {#regression}

### 三个误差与条件偏差–方差分解

训练集有 $D$ 个样本、$d$ 个输入坐标，设计矩阵 $X\in\mathbb R^{D\times d}$ 的第 $\mu$ 行是 $x_\mu^\top$。设

$$
y=Xw_*+\varepsilon,\qquad
\mathbb E\varepsilon=0,\qquad
\mathbb E\varepsilon\varepsilon^\top=\sigma_\varepsilon^2 I_D,
$$

$$
\mathcal L_\lambda(w)=\frac{1}{2D}\|Xw-y\|^2+
\frac\lambda2\|w\|^2,\qquad \lambda\geq0.
$$

噪声与输入独立，教师参数 $w_*$ 固定。测试输入独立于训练集，二阶矩 $\Sigma=\mathbb E[xx^\top]$。于是未加 $1/2$ 的测试均方误差为

$$
R_{\rm test}=\mathbb E_{x,\varepsilon_{\rm new}}
[(x^\top w-y_{\rm new})^2]
=\underbrace{(w-w_*)^\top\Sigma(w-w_*)}_{R_{\rm ex}}
+\sigma_\varepsilon^2.
$$

$R_{\rm ex}$ 是超额预测风险，$\|w-w_*\|^2$ 是参数误差；只有 $\Sigma=I_d$ 时两者相等。训练损失则用同一份参与拟合的数据计算，不能把上式的独立测试平均直接搬过去。

例如，测试输入几乎不沿某个方向变化，那么该方向参数即使估计得很差，对预测的影响仍可很小。反过来，测试分布偏重训练数据很少覆盖的方向时，训练误差可能严重低估实际风险。

#### 风险分解中究竟平均了什么

上面的测试风险是给定已经训练出的 $w$ 后，对一个新样本取平均。若还要讨论算法在不同训练集上的典型表现，就有三层随机性：测试输入、训练标签噪声、训练设计矩阵。顺序可以写为 $\mathbb E_X\mathbb E_{\varepsilon\mid X}\mathbb E_{x,\varepsilon_{\rm new}}$。条件化在 $X$ 上很有用，因为此时数据谱是确定的，可以先完成动力学计算，再研究随机谱的统计。

以 $e=w-w_*$ 为例，测试误差展开为

$$
\mathbb E[(x^\top e-\varepsilon_{\rm new})^2]
=e^\top\mathbb E[xx^\top]e
-2e^\top\mathbb E[x\varepsilon_{\rm new}]
+\mathbb E\varepsilon_{\rm new}^2.
$$

交叉项为零用到了新标签噪声与输入独立且均值为零。若噪声与输入相关、教师模型遗漏了系统性结构，第二项就不能被自动删除。这不是无关紧要的技术条件：我们往往把所有不能解释的残差都叫噪声，但只有真正条件均值为零的部分才符合这里的设定。

进一步记 $\bar w=\mathbb E_{\varepsilon\mid X}w$、$C_w=\operatorname{Cov}_{\varepsilon\mid X}(w)$，则

$$
\mathbb E_{\varepsilon\mid X}R_{\rm ex}
=(\bar w-w_*)^\top\Sigma(\bar w-w_*)
+\operatorname{Tr}(\Sigma C_w).
$$

这是给定设计矩阵后的偏差–方差分解。偏差不是一个与算法无关的“模型不够复杂”的量：有限时间、ridge、不可观测零空间和初始化都会进入 $\bar w-w_*$。同样，方差也不是参数数量本身，而是噪声通过学习算子后，在测试度量 $\Sigma$ 下的放大。

#### 损失为何是一个二次型

把经验损失逐项展开，得到

$$
\mathcal L_\lambda(w)
=\frac12w^\top(H+\lambda I)w-b^\top w
+\frac{y^\top y}{2D}.
$$

因此 $H+\lambda I$ 同时是 Hessian 和梯度流的线性恢复算子。$H$ 的大本征值表示沿该方向移动一点就会显著改变训练预测，损失曲率大，梯度也大；小本征值表示训练数据对该方向不敏感。将“小曲率”理解成“容易沿它快速走远”并不适用于未经预条件的梯度流，因为速度还由梯度本身决定。这里的小曲率恰恰意味着慢恢复。

## 2. 从梯度流到谱滤波：哪些方向先被恢复 {#spectral-dynamics}

先固定 $X$，每个经验模态就是一条确定的恢复通道。它的速率和信噪比可以分别求出，之后才需要随机矩阵理论。

### 矩阵方程、积分因子与 SVD

记 $H=X^\top X/D$、$b=X^\top y/D$，梯度流为

$$
\dot w=-(H+\lambda I)w+b.
$$

分解 $H=V\operatorname{diag}(s_i)V^\top$，其中 $s_i\geq0$。令 $z_i=v_i^\top w$、$b_i=v_i^\top b$，便有

$$
\dot z_i=-(s_i+\lambda)z_i+b_i,
$$

$$
\boxed{z_i(t)=e^{-(s_i+\lambda)t}z_i(0)+
\frac{1-e^{-(s_i+\lambda)t}}{s_i+\lambda}b_i.}
$$

若 $s_i=\lambda=0$，由 $Xv_i=0$ 知 $b_i=0$，该方向保持初始化。不要把上式写成 $0/0$ 后随意赋值。用矩阵指数与伪逆可以统一表达，但逐模写法更清楚地显示了可学习子空间。

无正则化时，非零模态的松弛时间 $\tau_i=s_i^{-1}$。大特征值方向先学会，小特征值方向很慢；这首先是**经验协方差方向**上的顺序，不等于任意输入坐标或空间频率上的顺序。

#### 积分因子与 SVD 的同一解释

单模态方程乘上 $e^{(s_i+\lambda)t}$ 后，左边成为 $d[e^{(s_i+\lambda)t}z_i]/dt$。从 0 积分到 $t$，再除回指数，才得到上面的解。这个过程也说明：初值和标签驱动分别产生齐次解与受迫解，不能只保留长期固定点而丢掉初始化。

令 $X=U\operatorname{diag}(\sqrt{Ds_i})V^\top$ 为薄 SVD。对非零模态，

$$
b_i=\sqrt{\frac{s_i}{D}}\,u_i^\top y,
\qquad
z_i(\infty)=\frac{\sqrt{s_i/D}}{s_i+\lambda}u_i^\top y.
$$

从参数坐标看，是用 $1/(s_i+\lambda)$ 反演曲率；从标签坐标看，是先把 $y$ 投影到左奇异向量，再映射回右奇异向量。二者是同一线性算子在两个空间中的表示。参数零空间和标签中不可拟合的子空间，也分别由右、左奇异向量补空间给出。

### 早停、ridge 与单模态停止时间

取 $w(0)=0$、$\lambda=0$，并令

$$
\eta_i=\frac{v_i^\top X^\top\varepsilon}{D},
\qquad b_i=s_i z_i^*+\eta_i,
\qquad
\mathbb E[\eta_i\eta_j\mid X]
=\frac{\sigma_\varepsilon^2}{D}s_i\delta_{ij}.
$$

于是

$$
z_i(t)-z_i^*=-e^{-s_it}z_i^*
+\frac{1-e^{-s_it}}{s_i}\eta_i,
$$

$$
\boxed{\mathbb E[(z_i(t)-z_i^*)^2\mid X]
=e^{-2s_it}(z_i^*)^2+
\frac{\sigma_\varepsilon^2}{Ds_i}(1-e^{-s_it})^2.}
$$

第一项是尚未恢复的信号，第二项是已经拟合的噪声。小 $s_i$ 使最终噪声误差很大，但也使噪声进入得很慢。这就是早停能够有益的最小机制。

若 $\Sigma$ 与 $H$ 共用本征向量，测试风险是这些式子乘以对应测试权重后求和。一般情况下应保留 $V^\top\Sigma V$：给定 $X$ 的偏差向量可在测试度量下产生交叉项。简单地把每个参数误差相加，已经额外假设了各向同性测试分布。

早停滤波器为 $g_t(s)=1-e^{-st}$，最终 ridge 的信号滤波器为 $g_\lambda(s)=s/(s+\lambda)$。在 $s\ll t^{-1}$ 时前者约为 $st$，在 $s\ll\lambda$ 时后者约为 $s/\lambda$，所以 $\lambda\sim t^{-1}$ 能匹配低谱端的量级；两种完整滤波曲线并不相同，也不存在适用于所有谱与目标的精确换算。

#### 噪声协方差怎样穿过学习算子

定义 $F_{t,\lambda}(H)$ 的本征值为
$[1-e^{-(s_i+\lambda)t}]/(s_i+\lambda)$，退化点按连续极限取值。零初始化解为 $w=F_{t,\lambda}(H)b$。因为 $b=Hw_*+X^\top\varepsilon/D$，所以

$$
\bar e=[F_{t,\lambda}(H)H-I]w_*,
\qquad
C_w=\frac{\sigma_\varepsilon^2}{D}
F_{t,\lambda}(H)H F_{t,\lambda}(H).
$$

令 $A=V^\top\Sigma V$，偏差在此基底的分量为 $\bar e_i$，便得到一般测试度量下的完整表达：

$$
\mathbb E[R_{\rm ex}\mid X]
=\sum_{ij}A_{ij}\bar e_i\bar e_j
+\frac{\sigma_\varepsilon^2}{D}\sum_i A_{ii}s_i
\left[\frac{1-e^{-(s_i+\lambda)t}}{s_i+\lambda}\right]^2.
$$

噪声在经验谱基底里不相关，但偏差是一个确定向量；因此一般测试协方差的非对角元仍会耦合偏差。条件化以后只把风险写成单个谱积分，需要额外使用各向同性、共同本征基或教师方向平均，不能仅凭 $H$ 可对角化就完成这一步。

长期 ridge 的逐模误差为

$$
\mathbb E[(z_i-z_i^*)^2\mid X]
=\frac{\lambda^2(z_i^*)^2}{(s_i+\lambda)^2}
+\frac{\sigma_\varepsilon^2s_i}{D(s_i+\lambda)^2}.
$$

增加 $\lambda$ 同时压低信号和噪声。它并不识别哪个模态“有用”，只按照曲率统一过滤；实际效果取决于目标能量是否集中在被保留的谱方向。

#### 单模态最优停止时间

对于一个正特征值，记 $A_i=(z_i^*)^2$、$B_i=\sigma_\varepsilon^2/(Ds_i)$，无 ridge 风险是 $E_i=A_iu^2+B_i(1-u)^2$，其中 $u=e^{-s_it}$。直接对 $u$ 求导，得到

$$
u_* =\frac{B_i}{A_i+B_i},\qquad
\boxed{t_i^*=\frac1{s_i}\log\left(1+
\frac{Ds_i(z_i^*)^2}{\sigma_\varepsilon^2}\right).}
$$

这个表达同时包含学习速度和该模态的信噪比。没有教师信号时 $t_i^*=0$；噪声趋于零时最优时间趋于无穷。多模态模型只有一个全局停止时间，一般无法同时取到所有 $t_i^*$，因此早停是一种受限的谱正则化，而不是逐模最优去噪。

若 $s_i$ 很小，同时 $Ds_iA_i\ll\sigma_\varepsilon^2$，对数展开给出 $t_i^*\simeq DA_i/\sigma_\varepsilon^2$。这一极限中，最优时刻本身未必随 $1/s_i$ 发散，但该时刻的信号恢复比例仍很小。慢模态可能根本不值得完整学习，这与“所有方向最终都要学完”的优化直觉不同。

#### Gaussian 先验、ridge 与 Bayes 最优滤波

若教师先验是 $w_*\sim\mathcal N(0,\tau^2I)$，标签噪声是方差 $\sigma_\varepsilon^2$ 的 Gaussian，负对数后验为

$$
\frac{\|Xw-y\|^2}{2\sigma_\varepsilon^2}
+\frac{\|w\|^2}{2\tau^2}.
$$

与本文平均平方损失比较，后验均值对应
$\lambda=\sigma_\varepsilon^2/(D\tau^2)$，后验协方差是
$\sigma_\varepsilon^2(X^\top X+\sigma_\varepsilon^2I/\tau^2)^{-1}$。
同一个 $\lambda$ 对所有模态都是合理的，来自教师各向同性先验这一额外信息；若教师方差依方向变化，Bayes 滤波也会变成方向依赖。频率或谱方向上的偏好，可以被理解为对教师结构的先验，而不只是一个技术性的稳定化参数。

### 零空间、训练残差与离散步长

令 $r=Xw-y$，在无正则化梯度流下

$$
\dot r=-\frac{XX^\top}{D}r,
\qquad r(t)=e^{-Kt}r(0),\qquad K=XX^\top/D.
$$

无法落入 $\operatorname{im}X$ 的标签分量不会消失。若 $X$ 满行秩、$d\geq D$，则训练集可以插值，但参数零空间里的教师分量依然没有被数据确定。

离散更新 $w_{n+1}=w_n-\eta(Hw_n-b)$ 中，每个齐次模态乘以 $1-\eta s_i$。收敛条件是 $0\lt \eta\lt 2/s_{\max}$；要让所有模态单调松弛还需 $\eta\leq1/s_{\max}$。连续时间解不能用于解释越过稳定边界后的振荡或发散。

### 有限谱的数值结果

![线性回归各模态的恢复与早停风险。](figures/linear_regression.svg)

*图 1 · 精确谱模型，自制计算。* 参数取 $s=(1,0.1,0.01)$、$z_i^*=1$、$D=100$、$\sigma_\varepsilon=0.5$。左图为信号滤波 $1-e^{-s_it}$；右图将条件期望中的偏差和方差分别求和，测试二阶矩取单位阵。右图没有重新抽训练集，也不是深度网络的实测曲线。

最慢方向是否含有教师信号，决定了晚期拟合是在恢复有用结构，还是主要吸收噪声。

作为数值核对，取 $s=0.01$、$t=100$：信号恢复比例约为 $0.632$，噪声项为 $0.25(1-e^{-1})^2\simeq0.100$。损失是否除以样本数会使噪声项相差 $D$ 倍。

## 3. 滤波器之外：参数度量与预测敏感度 {#regularization}

同一份数据上的曲率谱，还不足以独立决定统计结论。测试协方差规定错误的代价，参数坐标规定梯度使用的度量，而同一标签被用于拟合和评分，又会使训练误差产生系统性乐观。

### 小奇异值为什么既慢又危险

令 $A=X/\sqrt D$、$\tilde y=y/\sqrt D$。最小二乘估计实际要解 $Aw\simeq\tilde y$。写成奇异值展开，$A=\sum_i\sqrt{s_i}u_iv_i^\top$，无正则化的非零模态解为

$$
\widehat w_i=\frac{u_i^\top\tilde y}{\sqrt{s_i}}
=w_{*,i}+\frac{u_i^\top\varepsilon}{\sqrt{Ds_i}}.
$$

因此小 $s_i$ 有双重含义。它使梯度沿这个方向小，恢复速度慢；同时使反演因子 $s_i^{-1/2}$ 大，标签误差被放大。这两件事不是独立巧合，而是同一个观测算子在正向和逆向使用时的结果。训练得更久所获得的最后一部分拟合能力，恰好可能属于最不可靠的方向。

若只保留 $s_i\geq\tau$ 的模态，得到截断 SVD 估计，其滤波器是硬阈值；梯度流的滤波器是 $g_t(s)=1-e^{-st}$；ridge 的滤波器是 $g_\lambda(s)=s/(s+\lambda)$。三者的预测均可写成

$$
\widehat w=\sum_{s_i\gt 0}g(s_i)
\left(w_{*,i}+\frac{u_i^\top\varepsilon}{\sqrt{Ds_i}}\right)v_i.
$$

在测试协方差与 $H$ 共用本征向量、其本征值为 $\nu_i$ 的情况下，条件风险统一为

$$
R_{\rm ex}(g)=\sum_{s_i\gt 0}\nu_i
\left[(1-g(s_i))^2w_{*,i}^2+
\frac{\sigma_\varepsilon^2}{Ds_i}g(s_i)^2\right]
+\sum_{s_i=0}\nu_iw_{*,i}^2.
$$

这使“选择训练时间”变成“选择一个滤波函数”。三个方法都抑制弱观测方向，但并不逐模等价。$\lambda\sim1/t$ 只是转折尺度相似：在 $s=1/t$，梯度流的保留比例为 $1-e^{-1}$，取 $\lambda=1/t$ 的 ridge 则为 $1/2$。把早停叫作隐式正则化是合理的，把它当作完全相同的 ridge 则会丢掉滤波形状的信息。

### 改变坐标怎样改变路径与隐式几何

设先对白化坐标 $z=\Sigma^{-1/2}x$ 建模。预测函数保持相同时，系数变为 $\beta=\Sigma^{1/2}w$。在 $\beta$ 上做欧氏梯度下降，再变回 $w$，得到

$$
\dot w=-\Sigma^{-1}\nabla_w\mathcal L.
$$

这与原来 $\dot w=-\nabla_w\mathcal L$ 有相同的无约束最优预测，却有不同的路径与零空间选择。一般取对称正定迁移率 $M$，令 $u=M^{-1/2}w$，则

$$
\dot u=-M^{1/2}HM^{1/2}u+M^{1/2}b,
\qquad w=M^{1/2}u.
$$

真正控制时间尺度的于是变成 $M^{1/2}HM^{1/2}$ 的谱。零初始化的插值解还最小化 $w^\top M^{-1}w$。这说明数据预处理、参数归一化、预条件器与隐式偏置之间有直接的联系；它们并非只影响数值效率，有时也改变最终选择的函数。

当 $M=H^{-1}$ 且 $H$ 可逆时，各个经验模态的确定性收敛速度相同。然而这并没有消除观测噪声，长期解仍是同一个最小二乘解。加速弱方向也可能更快地把它们的噪声拟合进去。因此优化条件数变好与统计风险变好，仍需要分别计算。

### 帽子矩阵、自由度与训练误差的乐观偏差

对零初始化的线性平滑器，训练预测可写成 $\widehat y=S_ty$，其中

$$
S_t=U\operatorname{diag}(1-e^{-s_it})U^\top,
\qquad \operatorname{df}(t)=\operatorname{Tr}S_t
=\sum_{s_i\gt 0}(1-e^{-s_it}).
$$

这里的有效自由度从 0 连续增长到训练设计矩阵的秩，而不是一开始就等于参数数目。对于 ridge，$S_\lambda=X(X^\top X+D\lambda I)^{-1}X^\top$，相应自由度是 $\sum_i s_i/(s_i+\lambda)$。它度量预测对每个训练标签的总敏感度，能够刻画有限时间或有限正则化时实际使用了多少数据方向。

固定设计，令真实无噪声训练均值 $\mu=Xw_*$，再生成独立的同设计标签 $y'=\mu+\varepsilon'$。对任意固定的对称平滑矩阵 $S$，分别计算训练残差与独立重复标签误差：

$$
\frac1D\mathbb E\|Sy-y\|^2
=\frac1D\|(S-I)\mu\|^2
+\frac{\sigma_\varepsilon^2}{D}\operatorname{Tr}(S-I)^2,
$$

$$
\frac1D\mathbb E\|Sy-y'\|^2
=\frac1D\|(S-I)\mu\|^2
+\frac{\sigma_\varepsilon^2}{D}\operatorname{Tr}S^2
+\sigma_\varepsilon^2.
$$

两式相减，得到精确的乐观偏差

$$
\frac1D\mathbb E\big[\|Sy-y'\|^2-\|Sy-y\|^2\big]
=\frac{2\sigma_\varepsilon^2}{D}\operatorname{Tr}S.
$$

这给“训练误差为什么偏乐观”一个与动力学相连的答案：随着时间增加，平滑器对标签的响应增强，同一份噪声被同时用于拟合与评价，训练误差便少付出一项与自由度成正比的代价。它并不要求参数数目小于样本数，也不要求所有模态已经收敛。

不过这个结论只对应在原输入点重新采样标签的 fixed-design 预测，尚未处理新输入带来的设计变化。若测试协方差与训练覆盖不同，还要回到前面的 $\Sigma$ 加权风险。用自由度修正训练误差很有价值，但不能因此省去分布外方向的分析。

留一法还揭示了单个样本的杠杆作用。对固定 $\lambda\gt 0$ 的 ridge，删除第 $\mu$ 个样本时应保持求和形式的惩罚系数 $D\lambda$ 不变，才有

$$
y_\mu-\widehat y_\mu^{(-\mu)}
=\frac{y_\mu-\widehat y_\mu}{1-(S_\lambda)_{\mu\mu}}.
$$

分母接近零说明该点对自身拟合贡献很大，极小训练残差仍可能对应较大留一误差。在非均匀设计中，平均自由度会隐藏这些高杠杆样本。由此看，谱描述给出全局方向，杠杆值则补充这些方向在具体样本上的分布；二者对于理解数据覆盖都不可缺少。

## 4. 随机设计与插值：条件风险如何变成双下降 {#interpolation}

现在让输入独立 Gaussian，研究设计矩阵随机性。上一部分中的逆特征值会进入风险的设计平均；阈值附近的病态不再只是某个方向的偶然困难，而可能形成集体谱效应。以下仍固定任务和信号尺度，不将改变任务难度混入参数增长。

### Gaussian 设计与条件风险

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

### 欠参数侧的逆 Wishart 发散

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

#### 逆 Wishart 一阶矩可以直接推导

令 $W=X^\top X$。对 $D\geq d$ 的标准 Gaussian 设计，正定矩阵上的 Wishart 密度为

$$
p(W)\propto |W|^{(D-d-1)/2}e^{-\operatorname{Tr}W/2}.
$$

记 $a=(D-d-1)/2$。对一个对角矩阵元作积分分部，边界项在相应可积条件下消失，因而

$$
0=\int\frac{\partial p(W)}{\partial W_{ii}}dW
=a\,\mathbb E[(W^{-1})_{ii}]-\frac12.
$$

旋转不变性使 $\mathbb E W^{-1}$ 必为单位阵的倍数，于是得到 $\mathbb E W^{-1}=I/(D-d-1)$。分母中额外的 $-1$ 是有限样本效应；只从“自由度约为 $D-d$”的计数无法恢复它。

同时，这个积分并不是在 $D\gt d$ 时总收敛。矩阵可以几乎奇异，小概率事件对应极大的逆矩阵元；$D\leq d+1$ 时，逆的一阶矩不再有限。一次计算通常得到有限参数，不意味着跨设计矩阵的平均误差有限。插值附近需要明确报告均值、中位数或某个分位数，而不是将它们当作同一个极限。

更直接的几何图像来自最小奇异值。对标签噪声的一维投影 $\varepsilon_i$，反演产生的参数噪声为 $\varepsilon_i/\sqrt{Ds_i}$。当一个训练方向几乎由其他方向线性组合得到时，数据几乎无法区分沿该方向的大幅参数变化，却仍要求精确拟合标签；巨大的参数系数就是这个不适定反演的代价。

### 过参数侧的投影偏差与噪声方差

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

#### 为什么“最小范数”仍然留下偏差

把教师写成 $w_*=Pw_*+(I-P)w_*$。前一部分可由训练行空间观察，后一部分对所有训练输入都没有输出。最小范数规则将后一部分设为零，不是因为知道它在真实教师里不存在，而是因为训练数据没有提供关于它的信息。

对 Gaussian 设计，行空间是均匀随机的 $D$ 维子空间。由投影矩阵 $P^2=P$、$\operatorname{Tr}P=D$ 和旋转不变性，$\mathbb EP=(D/d)I$。因此

$$
\mathbb E\|(I-P)w_*\|^2
=w_*^\top(I-\mathbb EP)w_*
=r^2(1-D/d).
$$

这一步没有把教师随机化，只平均了随机训练子空间。若教师与设计矩阵存在共同生成的结构，或者某些输入方向方差特别大，均匀投影的结论就会改变。良好的过参数泛化常常与“信号集中在可恢复方向、剩余方向主要吸收噪声”有关，而不是仅由零空间维数决定。

高维过参数风险对 $\gamma$ 的导数为

$$
\frac{dR_{\rm ex}}{d\gamma}
=\frac{r^2}{\gamma^2}-\frac{\sigma^2}{(\gamma-1)^2}.
$$

若 $r\leq\sigma$，这一侧单调下降；若 $r\gt \sigma$，内部极小点为 $\gamma_*=r/(r-\sigma)$，之后风险又向 $r^2$ 上升。极端增大模型会恢复到接近“什么都没有学到”的风险基线。这说明即使存在第二次下降，最优点也未必在无限参数处。

## 5. 随机矩阵谱：峰值、零空间与正则化交叉区 {#random-matrices}

有限样本逆矩给出精确平均，连续谱则说明大量方向如何共同产生这个平均。二者的极限和修正要分别保留。

### Marchenko–Pastur 谱与 resolvent

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

#### Resolvent 的分支与插值奇点

这里将 resolvent 定义在正实正则化参数上：$m(\lambda)=d^{-1}\operatorname{Tr}(H+\lambda I)^{-1}$。MP 极限满足二次方程

$$
\gamma\lambda m^2+(\lambda+1-\gamma)m-1=0.
$$

取 $m\gt 0$ 且 $m\sim1/\lambda$ 的分支，得到

$$
m(\lambda)=\frac{-(\lambda+1-\gamma)
+\sqrt{(\lambda+1-\gamma)^2+4\gamma\lambda}}
{2\gamma\lambda}.
$$

平方根的选取是物理条件的一部分，另一根不能表示正定矩阵的逆迹。对 $\gamma\lt 1$，$m(0)=1/(1-\gamma)$；对 $\gamma\gt 1$，$m(\lambda)\sim(1-1/\gamma)/\lambda$，发散系数恰好是零谱原子的质量。$\gamma=1$ 则有 $m(\lambda)\sim\lambda^{-1/2}$，对应连续谱直接接触原点。

把 $m$ 代回 ridge 风险，可以在一个表达式中同时保留偏差、噪声与零空间。先令 $\lambda\to0$ 再靠近 $\gamma=1$，和固定 $\lambda\gt 0$ 穿越插值阈值，得到不同的曲线。正则化不是事后修饰尖峰，而是改变了对整个小谱端的探测分辨率。

#### 同一个嵌套任务怎样产生完整双下降

前面固定信号范数的模型族，在欠参数侧只有方差上升。若要研究同一个任务随可用特征增多的变化，可以取一组无限正交 Gaussian 坐标，教师为 $y=\sum_{j\geq1}\theta_jx_j+\varepsilon$，学生只用前 $d$ 个坐标。记尾部能量 $T_d=\sum_{j\gt d}\theta_j^2$。

未使用坐标与已使用坐标独立，所以对该学生而言，训练有效噪声方差是 $v_d=\sigma^2+T_d$。当 $D\gt d+1$ 时，相对完整教师的超额预测风险为

$$
\boxed{R_d=T_d+(\sigma^2+T_d)\frac{d}{D-d-1}.}
$$

第一个项随模型增大下降，第二项在插值附近上升。这里经典的左侧下降段来自遗漏信号减少，插值尖峰来自同时拟合不可预测的有效噪声，两者有不同来源。

当 $d\gt D+1$ 时，同一任务的风险变成

$$
R_d=T_d+\left(1-\frac Dd\right)\sum_{j\leq d}\theta_j^2
+(\sigma^2+T_d)\frac{D}{d-D-1}.
$$

现在信号尾部、已表示信号的不可恢复投影和噪声放大同时竞争。这一式子解释为什么增加模型、增加独立数据、减少噪声不能简单地视为同一个控制方向。它也说明双下降曲线的具体形状需要指定教师系数序列，而不是由“过参数化”四个字自动给出。

#### 泛化良好与风险下降的区别

一个风险从很大的插值尖峰降到一个仍然很大的平台，也属于第二次下降，但不能据此称为良好泛化。接近噪声底、低于零预测器、接近 Bayes 风险，是更强且不同的要求。各向同性模型中，未恢复信号比例通常不会因为无限加宽自行消失；具有结构的协方差谱才可能改变这个结论。随机矩阵理论真正提供的是可计算的比较标准，而不是“越大越好”的普遍口号。

### 均值风险为什么不等于典型一次运行

给定设计矩阵，噪声风险含有 $\operatorname{Tr}(H^{-1})$ 或非零谱上的伪逆迹。把它再对随机设计取期望，得到的是总体平均曲线。这条平均曲线不必等于一组中等规模实验中最常看到的曲线：在阈值附近，少数异常小的最小奇异值就能给均值带来很大贡献。

以 $d=1$、$D$ 个标准 Gaussian 输入为例，$X^\top X\sim\chi_D^2$，于是纯噪声估计的条件方差为 $\sigma_\varepsilon^2/(X^\top X)$。其期望存在当且仅当 $D\gt 2$，并有

$$
\mathbb E\frac{1}{\chi_D^2}=\frac{1}{D-2}.
$$

但对任意一份非零设计，条件方差都是有限的。期望无穷说的是随机设计分布的尾部足够重，并不是每次实验都会输出无穷大。高维 Wishart 的 $D-d-1$ 分母正是同一种现象在多个方向上的推广。

实验中因此应同时保留均值、中位数和分位区间，并说明是否截去了极端值。只报告平均值会隐藏典型性，只报告中位数又无法检验均值公式。二者回答不同问题：典型一次运行有多大风险，和重复随机设计后平均会付出多少风险。

### 正则化怎样控制插值附近的小谱

令 $H$ 的本征值为 $s_i$，各向同性测试和 ridge 下的噪声项为

$$
R_{\rm var}(\lambda\mid X)
=\frac{\sigma_\varepsilon^2}{D}
\sum_i\frac{s_i}{(s_i+\lambda)^2}.
$$

对任意 $s\geq0$，函数 $s/(s+\lambda)^2$ 在 $s=\lambda$ 达到最大值 $1/(4\lambda)$，因此有限维有界估计为

$$
R_{\rm var}(\lambda\mid X)
\leq\frac{\sigma_\varepsilon^2\operatorname{rank}(X)}{4D\lambda}.
$$

这不是最紧的泛化界，但清楚说明 $\lambda\gt 0$ 切断了无正则化的小谱发散。相应偏差则为 $\lambda^2w_*^\top(H+\lambda I)^{-2}w_*$，包含零空间中无法恢复的信号。降低方差需要付出怎样的偏差，仍取决于教师在哪些方向上有能量。

在高维等比例极限，若 $\rho_\gamma(s)$ 是经验协方差的极限谱，定义

$$
m_\gamma(\lambda)=\int\frac{\rho_\gamma(s)}{s+\lambda}\,ds,
\qquad -m_\gamma'(\lambda)=
\int\frac{\rho_\gamma(s)}{(s+\lambda)^2}\,ds.
$$

于是噪声风险的谱积分为

$$
R_{\rm var}(\lambda)
\simeq\sigma_\varepsilon^2\gamma
\left[m_\gamma(\lambda)+\lambda m_\gamma'(\lambda)\right].
$$

这个表达比只记插值两侧的分式更容易推广。给定一个可求的 resolvent，就可以得到正则化曲线；同时还可以检查 $\lambda\to0$ 与 $\gamma\to1$ 的极限是否交换。对插值临界谱 $\rho(s)\sim s^{-1/2}$，代入 $s=\lambda u$ 可见方差积分按 $\lambda^{-1/2}$ 增长。有限正则化平滑了峰，却没有抹去临界小谱的痕迹。

### 高维公式与有限样本的数值比较

![Gaussian 最小范数回归的风险与小谱端。](figures/double_descent.svg)

*图 2 · 自制 Gaussian 回归实验。* 左图为 $r^2=1,\sigma^2=0.25$ 的高维公式，并叠加 $D=80$、每个维数 24 次独立设计的风险中位数。误差直接用 $\|\hat w-w_*\|^2$ 测量，不引入有限测试集误差。右图为三个 $\gamma$ 的 MP 连续密度；零原子单独在图例标注，未画成无限高竖线。

无噪声时，插值尖峰消失，而过参数偏差仍然存在。ridge 改变的是估计器及其整个谱滤波，不能等同于将原风险曲线的尖峰截平。

一个清楚的算例是 $D=100,d=200$：偏差为 $0.5r^2$，噪声为 $100\sigma^2/99$。若 $d=1000$，噪声降至 $100\sigma^2/899$，偏差却增到 $0.9r^2$。这说明“更多参数有益”需要说明信噪比和比较区间。

## 6. 比较曲线之前，先比较任务与观测量 {#risk-interpretation}

### 参数数量的横轴究竟改变了什么

固定 $\|w_*\|$ 并改变 $d$，与固定一个无限维教师、逐步纳入更多坐标，不是同一个实验。前者把相同总信号放进不同维数的空间；后者使遗漏目标能量随模型规模下降。在第二种设定下，欠参数侧可能出现熟悉的偏差下降，而在前一种设定下未必存在那一段。

还有一种常见比较是固定每个坐标的信号方差。这时 $\|w_*\|^2$ 随 $d$ 增长，若输入归一化也不改变，任务本身的输出方差就增加。风险的纵轴会同时包含信号尺度变化。将输入方差设为 $1/d$ 可以控制输出尺度，却又改变协方差谱和优化时间的单位。因此跨模型的曲线应先核对信号、输入、噪声和损失的归一化。

同样，横轴若是神经网络参数数目，未必对应线性模型中的有效特征秩。重复神经元、参数对称性、隐藏层冗余和数据限制都会使二者偏离。插值阈值是关于实际约束与有效自由方向的关系，不能仅凭参数数目等于样本数来定位。

### 同谱数据仍可能有不同目标难度

Gaussian 各向同性模型之所以能将风险写成几个标量，是因为教师方向相对随机子空间没有特殊地位。真实任务中，大特征值方向可能携带更多信号，也可能主要携带干扰。把目标能量分布写成谱测度

$$
\nu_*(ds)=\sum_i(v_i^\top w_*)^2\delta_{s_i}(ds)
$$

以后，偏差与总谱密度分离：噪声往往由谱迹控制，信号偏差则由 $\nu_*$ 控制。两个协方差矩阵可以有相同特征值，却因教师方向不同而产生不同风险。由此看，double descent 是数据谱与估计规则能产生的一类现象，而不是排除了数据语义以后的普适单参数定律。

有限时间、无限时间与高维极限还存在次序问题。在给定数据上，有限时间会抑制慢模态；长期最小范数解则完全释放所有可观测方向。若临界小谱越来越接近零，先让维数无穷与先训练到收敛，未必得到同一风险。后面的 Scaling Laws 将以此讨论训练与测试可能具有不同时间指数。

## 相关文献

这里精确解释的是固定设计矩阵、平方损失下的优化与噪声拟合。它没有产生新特征，没有解释非线性网络的表征重组。随机训练集的风险平均，使用本篇后半部分的随机矩阵工具。

**相关文献。** [Advani 与 Saxe（2017），High-dimensional dynamics of generalization error in neural networks](https://arxiv.org/abs/1710.03667)把线性可解动力学与高维谱联系起来。本文的有限维计算直接从上述损失推导；该预印本的作者为两人，后续扩展版本的书目信息不与之混用。

小奇异值噪声放大是可推广的机制；Gaussian 逆 Wishart 公式不是任意相关、重尾输入的通式。特征学习也会改变谱与教师对齐，不能把训练后的非线性网络直接替换成独立 Gaussian 矩阵。

**相关文献。** [Hastie、Montanari、Rosset 与 Tibshirani，Surprises in High-Dimensional Ridgeless Least Squares Interpolation](https://arxiv.org/abs/1903.08560)给出更一般的高维风险结果；[Mei 与 Montanari，随机特征回归的精确渐近](https://arxiv.org/abs/1908.05355)研究非线性固定特征。这里的线性推导不冒充后者的全部理论。

## 结语：从一条恢复通道到随机高维几何

有限维梯度流提供逐模答案，随机设计平均解释这些模态怎样集体组织。它们共享同一个观测算子，却回答不同层次的问题。早停控制尚未释放的噪声，ridge 改变每个方向的收缩，过参数零空间留下未恢复信号，测试分布再决定这些错误是否重要。

因此 double descent 可以放回完整学习过程理解：插值并没有消除统计问题，而是改变了约束与解空间的几何。判断风险变化，需要把谱、目标投影与解选择一起计算。这个基线也会继续进入随机特征、核动力学和深度因子化模型，但每次替换都需要核对新的假设。

---

[系列目录](/notes/neural_network_theory/) · [Note 02](/notes/neural_network_theory/kernel_methods/) →
