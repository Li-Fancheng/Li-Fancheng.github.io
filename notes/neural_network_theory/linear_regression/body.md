> **I · 泛化与高维模型 · Note 01**
>
> [系列目录](/notes/neural_network_theory/) · [Note 02](/notes/neural_network_theory/double_descent/) →

**可解性层次：** 有限维精确可解。

## 1. 现象：训练得更久，为什么测试误差可能回升

一个模型可以不断降低训练损失，同时逐渐拟合标签里的随机噪声。要理解这个分离，不必从深层网络开始。线性回归已经包含三个关键对象：数据允许观察的方向、每个方向的恢复速度，以及测试分布怎样给这些方向加权。

本篇只用平方损失和连续时间梯度流，先把每个方向精确求解。它是后面 NTK、double descent 和时间幂律共同使用的动力学基线。这里的“时间”是吸收了学习率的优化时间，没有物理秒的含义。

## 2. 最小模型与三个不同的误差

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

### 风险分解中究竟平均了什么

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

### 损失为何是一个二次型

把经验损失逐项展开，得到

$$
\mathcal L_\lambda(w)
=\frac12w^\top(H+\lambda I)w-b^\top w
+\frac{y^\top y}{2D}.
$$

因此 $H+\lambda I$ 同时是 Hessian 和梯度流的线性恢复算子。$H$ 的大本征值表示沿该方向移动一点就会显著改变训练预测，损失曲率大，梯度也大；小本征值表示训练数据对该方向不敏感。将“小曲率”理解成“容易沿它快速走远”并不适用于未经预条件的梯度流，因为速度还由梯度本身决定。这里的小曲率恰恰意味着慢恢复。

## 3. 从矩阵方程到独立松弛模态

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

### 积分因子与 SVD 的同一解释

单模态方程乘上 $e^{(s_i+\lambda)t}$ 后，左边成为 $d[e^{(s_i+\lambda)t}z_i]/dt$。从 0 积分到 $t$，再除回指数，才得到上面的解。这个过程也说明：初值和标签驱动分别产生齐次解与受迫解，不能只保留长期固定点而丢掉初始化。

令 $X=U\operatorname{diag}(\sqrt{Ds_i})V^\top$ 为薄 SVD。对非零模态，

$$
b_i=\sqrt{\frac{s_i}{D}}\,u_i^\top y,
\qquad
z_i(\infty)=\frac{\sqrt{s_i/D}}{s_i+\lambda}u_i^\top y.
$$

从参数坐标看，是用 $1/(s_i+\lambda)$ 反演曲率；从标签坐标看，是先把 $y$ 投影到左奇异向量，再映射回右奇异向量。二者是同一线性算子在两个空间中的表示。参数零空间和标签中不可拟合的子空间，也分别由右、左奇异向量补空间给出。

## 4. 早停怎样形成谱滤波

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

### 噪声协方差怎样穿过学习算子

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

### 单模态最优停止时间

对于一个正特征值，记 $A_i=(z_i^*)^2$、$B_i=\sigma_\varepsilon^2/(Ds_i)$，无 ridge 风险是 $E_i=A_iu^2+B_i(1-u)^2$，其中 $u=e^{-s_it}$。直接对 $u$ 求导，得到

$$
u_* =\frac{B_i}{A_i+B_i},\qquad
\boxed{t_i^*=\frac1{s_i}\log\left(1+
\frac{Ds_i(z_i^*)^2}{\sigma_\varepsilon^2}\right).}
$$

这个表达同时包含学习速度和该模态的信噪比。没有教师信号时 $t_i^*=0$；噪声趋于零时最优时间趋于无穷。多模态模型只有一个全局停止时间，一般无法同时取到所有 $t_i^*$，因此早停是一种受限的谱正则化，而不是逐模最优去噪。

若 $s_i$ 很小，同时 $Ds_iA_i\ll\sigma_\varepsilon^2$，对数展开给出 $t_i^*\simeq DA_i/\sigma_\varepsilon^2$。这一极限中，最优时刻本身未必随 $1/s_i$ 发散，但该时刻的信号恢复比例仍很小。慢模态可能根本不值得完整学习，这与“所有方向最终都要学完”的优化直觉不同。

### Gaussian 先验、ridge 与 Bayes 最优滤波

若教师先验是 $w_*\sim\mathcal N(0,\tau^2I)$，标签噪声是方差 $\sigma_\varepsilon^2$ 的 Gaussian，负对数后验为

$$
\frac{\|Xw-y\|^2}{2\sigma_\varepsilon^2}
+\frac{\|w\|^2}{2\tau^2}.
$$

与本文平均平方损失比较，后验均值对应
$\lambda=\sigma_\varepsilon^2/(D\tau^2)$，后验协方差是
$\sigma_\varepsilon^2(X^\top X+\sigma_\varepsilon^2I/\tau^2)^{-1}$。
同一个 $\lambda$ 对所有模态都是合理的，来自教师各向同性先验这一额外信息；若教师方差依方向变化，Bayes 滤波也会变成方向依赖。频率或谱方向上的偏好，可以被理解为对教师结构的先验，而不只是一个技术性的稳定化参数。

## 5. 训练残差、零空间与离散梯度下降

令 $r=Xw-y$，在无正则化梯度流下

$$
\dot r=-\frac{XX^\top}{D}r,
\qquad r(t)=e^{-Kt}r(0),\qquad K=XX^\top/D.
$$

无法落入 $\operatorname{im}X$ 的标签分量不会消失。若 $X$ 满行秩、$d\geq D$，则训练集可以插值，但参数零空间里的教师分量依然没有被数据确定。

离散更新 $w_{n+1}=w_n-\eta(Hw_n-b)$ 中，每个齐次模态乘以 $1-\eta s_i$。收敛条件是 $0\lt \eta\lt 2/s_{\max}$；要让所有模态单调松弛还需 $\eta\leq1/s_{\max}$。连续时间解不能用于解释越过稳定边界后的振荡或发散。

## 6. 数值结果

![线性回归各模态的恢复与早停风险。](figures/experiment.svg)

*图 1 · 精确谱模型，自制计算。* 脚本使用 $s=(1,0.1,0.01)$、$z_i^*=1$、$D=100$、$\sigma_\varepsilon=0.5$。左图为信号滤波 $1-e^{-s_it}$；右图将条件期望中的偏差和方差分别求和，测试二阶矩取单位阵。右图没有重新抽训练集，也不是深度网络的实测曲线。

复现命令为 `python experiments.py --note 1`，脚本见[系列目录](/notes/neural_network_theory/#reproduce)。最慢方向是否含有教师信号，决定了晚期拟合是在恢复有用结构，还是主要吸收噪声。

作为数值核对，取 $s=0.01$、$t=100$：信号恢复比例约为 $0.632$，噪声项为 $0.25(1-e^{-1})^2\simeq0.100$。损失是否除以样本数会使噪声项相差 $D$ 倍。

## 7. 解释范围与相关文献

这里精确解释的是固定设计矩阵、平方损失下的优化与噪声拟合。它没有产生新特征，没有解释非线性网络的表征重组。随机训练集的风险平均，还需要下一篇的随机矩阵工具。

**相关文献。** [Advani 与 Saxe（2017），High-dimensional dynamics of generalization error in neural networks](https://arxiv.org/abs/1710.03667)把线性可解动力学与高维谱联系起来。本文的有限维计算直接从上述损失推导；该预印本的作者为两人，后续扩展版本的书目信息不与之混用。

---

[系列目录](/notes/neural_network_theory/) · [Note 02](/notes/neural_network_theory/double_descent/) →
