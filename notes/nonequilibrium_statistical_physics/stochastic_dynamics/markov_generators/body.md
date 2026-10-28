若想知道一个粒子多久以后到达目标，直接求出所有时刻的概率密度，再从边界流提取事件统计，往往做了很多不必要的工作。Markov 结构提供另一条路线：从给定起点预测未来，利用一个很短时间步的条件平均，直接为目标量写出方程。

生成元把局部随机演化规则压缩成一个算子。它既控制概率密度，又控制观测量、时间关联和首次到达。本篇会在这些具体计算中建立统一结构，并解释为什么同一个双稳过程会同时出现很快和极慢的模式。

## 4.1 Markov 性：当前状态何时足以预测未来

对当前状态 $X_t$，Markov 性要求：给定 $X_t$ 后，未来的条件分布不再依赖更早的完整历史。时间齐次时，转移核只依赖时间差，记为 $P_t(dy\mid x)$。若有密度则写成 $P_t(y\mid x)dy$。

这是关于状态选择的性质。只记录惯性粒子的位置，短时间未来还取决于当前速度，位置本身通常不是 Markov 过程；保留 $(X,V)$ 才能闭合。粗粒化后若速度足够快，位置在较长时间尺度上才可近似成为过阻尼 Markov 扩散。观测隐藏了慢变量时，也会留下类似记忆。

时间齐次 Markov 核满足 Chapman–Kolmogorov 关系

$$
P_{t+s}(dy\mid x)=\int P_t(dz\mid x)P_s(dy\mid z).
$$

它先对中间状态求和，再接上下一段未来。对可观测量 $f$ 定义

$$
(P_tf)(x)=\mathbb E_x[f(X_t)],
$$

便有 $P_{t+s}=P_tP_s$、$P_0=\mathrm{Id}$。这一族算子称为 Markov 半群。它把任意有限时间的预测拆成了局部演化的连续复合。

## 4.2 生成元怎样编码局部随机规律

对适当定义域中的函数，生成元为

$$
\boxed{(\mathcal Lf)(x)
=\lim_{h\downarrow0}\frac{\mathbb E_x[f(X_h)]-f(x)}{h}.}
$$

它读取从 $x$ 出发的一小段时间内，$f$ 的平均变化率。对于

$$
dX=b(X)dt+\sigma(X)dW,
$$

由 Itô 公式立即得到

$$
\mathcal Lf=b f'+\frac12\sigma^2f''.
$$

漂移贡献一阶方向变化，随机增量的二次变差贡献二阶扩散。取 $f=1$ 得到 $\mathcal L1=0$，对应守恒过程的总概率保持；取 $f=x,x^2$ 则立即返回矩方程。

生成元的微分表达式、作用域和边界条件共同定义过程。相同内部算子配上吸收或反射边界，就会有不同的未来统计；因此“给出 $b,D$ 就完全给出过程”还默认了边界、初始准备与适定性。对于这里的 OU 和光滑受限双稳扩散，这些条件可按前文自然边界或指定有限区间边界来建立。[Pavliotis 教材](https://www.ma.imperial.ac.uk/~pavl/PavliotisBook.pdf)第 2 章系统讨论了半群、生成元与 Kolmogorov 方程。

### 前向与后向分别在演化什么

固定初始密度 $p_0$，终点密度满足

$$
\partial_t p=\mathcal L^\dagger p,
\qquad
\mathcal L^\dagger p=-\partial_x(bp)+\frac12\partial_x^2(\sigma^2p).
$$

这就是 Fokker–Planck 方程。若改为固定终点观测量 $f$，把每个起点的未来平均记为 $u(x,t)=\mathbb E_x[f(X_t)]$，则

$$
\partial_tu=\mathcal Lu,\qquad u(x,0)=f(x).
$$

后一个方程的空间导数作用于起点 $x$，因此叫后向方程；这里的 $t$ 是预测时长，仍然向正方向增加。若改用“日历时间到固定终点 $T$”的记法 $v(x,s)=\mathbb E[f(X_T)\mid X_s=x]$，则写成 $\partial_sv+\mathcal Lv=0$，终值为 $v(x,T)=f(x)$。

二者都计算同一个平均：

$$
\int f(x)p(x,t)dx=\int (P_tf)(x)p_0(x)dx.
$$

前向把初始概率传播到各终点，后向把终点观测量传播回各起点。选择哪一个取决于所求问题，而不是时间是否真的倒流。

## 4.3 首次到达时间：直接为一个事件写方程

设目标位置为 $r$，粒子在区间 $[\ell,r]$ 中运动，左端反射、右端吸收。定义

$$
\tau_r=\inf\{t\ge0:X_t=r\},\qquad T(x)=\mathbb E_x\tau_r.
$$

在内部点先走很短的 $h$，已经花掉这段时间，接下来平均还需 $T(X_h)$。利用 Markov 性，

$$
T(x)=h+\mathbb E_x[T(X_h)]+o(h).
$$

展开条件平均 $\mathbb E_xT(X_h)=T(x)+h\mathcal LT(x)+o(h)$，消去 $T(x)$，得到

$$
\boxed{\mathcal LT=-1,\qquad T(r)=0,\qquad T'(\ell)=0.}
$$

吸收端的剩余时间为零；反射端对后向算子给出 Neumann 条件。这个推导也可由 Dynkin 公式与停止时间论证严格化，条件包括所需的正则性和首达时间可积性。

### 一维梯度扩散的积分解

对 $b=-U'$、扩散常数 $D$，方程为

$$
DT''-U'T'=-1.
$$

乘以积分因子 $e^{-U/D}$，有

$$
\frac d{dx}\left(e^{-U(x)/D}T'(x)\right)
=-\frac1D e^{-U(x)/D}.
$$

先从反射端 $\ell$ 积分到 $y$，再利用 $T(r)=0$ 积分，得到

$$
\boxed{
T(x)=\frac1D\int_x^rdy\,e^{U(y)/D}
\int_{\ell}^y dz\,e^{-U(z)/D}.}
$$

这个表达式在指定有限区间和边界下是精确的，不需要弱噪声。内层积分累积出发阱里的统计质量，外层积分对高势垒附近给予大权重。势阱占据和穿越障碍在同一公式中共同决定等待时间。

### 从精确积分得到 Kramers 前因子

令 $x_a$ 为出发阱底，$x_s$ 为相邻势垒顶。左端反射位置足够远离阱底，右端目标位于势垒之后、另一个阱内；两端与有关驻点保持固定距离。在 $D\to0$ 且极值非退化时，

$$
U(z)\simeq U(x_a)+\frac12U''(x_a)(z-x_a)^2,
$$

$$
U(y)\simeq U(x_s)-\frac12|U''(x_s)|(y-x_s)^2.
$$

内、外积分分别由一个完整高斯峰主导，得到

$$
\boxed{
T(x_a)\simeq
\frac{2\pi}{\sqrt{U''(x_a)|U''(x_s)|}}
e^{[U(x_s)-U(x_a)]/D}.}
$$

对本文双稳势，$U''(-1)=2$、$U''(0)=-1$，所以

$$
T_{L\to R}\simeq\sqrt2\pi\,e^{1/(4D)},\qquad
k\simeq\frac{1}{\sqrt2\pi}e^{-1/(4D)}.
$$

这补上了路径大偏差未给出的曲率前因子：阱底的回复强度与势垒顶的不稳定性共同影响尝试越垒的时间尺度。[Berglund 的综述](https://arxiv.org/abs/1106.5799)给出了这一公式及其多维推广的条件。

如果目标恰好就是势垒顶 $r=x_s$，外层 Gaussian 只积分到峰顶，前因子变成 $\pi/\sqrt{U''(x_a)|U''(x_s)|}$。第一次触及分界面与成功进入另一势阱是两个事件；在对称局部极限中，到顶后的运动仍可能返回，因而它们的平均时间不能使用同一个前因子。

![首达时间与双稳扩散的弛豫谱](../figures/05_first_passage_and_spectrum.svg)

*图 5｜左图比较双稳势从 $-1$ 到 $0.8$ 的精确积分首达时间和弱噪声 Kramers 近似，左端 $-2$ 反射。指数尺度一致，前因子近似在较弱噪声时更准确。右图比较 OU 的整数间隔弛豫率与 $D=0.05$ 双稳扩散的前几个非零弛豫率；双稳最慢率远小于其余模式。*

## 4.4 到达哪个目标：分裂概率

现在换一个边界问题：$\ell,r$ 两端都吸收。令

$$
q(x)=\Pr_x(\tau_r\lt\tau_{\ell}).
$$

这次短时间步不累积等待时间，因而 $q(x)=\mathbb E_x[q(X_h)]+o(h)$，得到

$$
\mathcal Lq=0,\qquad q(\ell)=0,\qquad q(r)=1.
$$

对于梯度扩散，解为

$$
\boxed{q(x)=
\frac{\int_{\ell}^{x}e^{U(y)/D}dy}
{\int_{\ell}^{r}e^{U(y)/D}dy}.}
$$

在双稳系统中，它衡量当前状态对两个去向的承诺程度，常称 committor。弱噪声时，$q$ 通常主要在势垒附近从接近 $0$ 变到接近 $1$：一旦已经深入某个阱，接下来的去向较确定；势垒附近则对微小噪声很敏感。

同一个 $\mathcal L$，右端为 $-1$ 时计算平均时间，右端为 $0$ 且边界值不同则计算事件概率。局部规律没有改变，目标观测量改变了方程与边界条件。

## 4.5 弛豫谱：不同时间尺度从哪里出现

若生成元具有合适的离散谱，可将后向演化的某个模式写成

$$
\mathcal L\phi_n=-\lambda_n\phi_n,\qquad
P_t\phi_n=e^{-\lambda_nt}\phi_n.
$$

零模对应常数观测量与总概率守恒，非零 $\lambda_n$ 决定不同形状的偏离衰减多快。并非所有 Markov 过程都有可直接使用的纯离散谱；无穷空间的自由扩散等问题还需要连续谱。这里先看受限的可逆模型。

### OU 的整数谱与可观测量

把 OU 坐标缩放为 $y=\sqrt{\gamma/D}\,x$，则

$$
\mathcal L=\gamma(\partial_y^2-y\partial_y).
$$

概率论约定的 Hermite 多项式 $\mathrm{He}_n$ 满足

$$
(\partial_y^2-y\partial_y)\mathrm{He}_n=-n\mathrm{He}_n,
\qquad \lambda_n=n\gamma.
$$

前两个非平凡模式可直接在原坐标验证：

$$
\mathcal Lx=-\gamma x,\qquad
\mathcal L(x^2-D/\gamma)=-2\gamma(x^2-D/\gamma).
$$

所以均值和方差偏离的两种弛豫率来自不同的谱模式。在稳态加权内积中展开 $f-\langle f\rangle_{\mathrm{ss}}=\sum_{n\ge1}c_n\phi_n$，若模式正交归一，则

$$
\langle f(X_t)f(X_0)\rangle_c
=\sum_{n\ge1}c_n^2e^{-\lambda_nt}.
$$

时间关联既取决于系统有哪些模式，也取决于观测量与哪些模式重叠。OU 的位置只投影到第一模，因此关联是单指数；$x^2$ 的连通关联读取第二模，衰减率为 $2\gamma$。

### 双稳系统的慢模是阱间质量不平衡

弱噪声双稳扩散中，第一非零率近似为

$$
\lambda_1\simeq k_{LR}+k_{RL},
$$

其余主要阱内模式的率保持在局部回复尺度。最慢的后向特征函数在两个阱内分别近似常数、符号相反，读取“粒子属于哪个阱”；相应前向模式则把概率从一个阱搬到另一个阱，而保持总质量不变。

因此阱内初始形状很快被遗忘，阱间占据差却长期保留。第二篇看到的概率输运、第三篇计算的势垒代价和这里的谱隙，描述了同一慢过程的不同方面。

对于常数 $D$ 的梯度扩散，还可以令 $p=e^{-U/(2D)}\psi$，把前向方程变换为

$$
\partial_t\psi=-H\psi,\qquad
H=-D\partial_x^2+\frac{(U')^2}{4D}-\frac12U''.
$$

在对应的边界条件下，$H$ 是自伴的 Schrödinger 型算子，基态为 $e^{-U/(2D)}$、本征值为零。它提供了计算弛豫谱的另一种方式。这里的时间演化是指数衰减，物理对象仍然是概率分布。

## 4.6 离散状态空间：主方程自然来自同一生成元

连续扩散不是 Markov 过程的唯一形式。若系统处于离散状态 $i$，以速率 $k_{ij}$ 从 $i$ 跳到 $j$，则

$$
(\mathcal Lf)_i=\sum_{j\ne i}k_{ij}(f_j-f_i).
$$

每项都是“发生速率乘以观测量的变化”。对偶的概率方程为

$$
\dot p_i=\sum_{j\ne i}(k_{ji}p_j-k_{ij}p_i).
$$

这是主方程。若把 $p$ 写成列向量，矩阵的第 $i,j$ 个非对角元应为从 $j$ 到 $i$ 的速率 $k_{ji}$，每列和为零；后向生成元则作用于状态函数。说明向量和指标约定可以避免把速率方向写反。

双稳系统的两态约化正是这个形式，但它只在阱内混合显著快于阱间跳跃时合理。若在到达势垒附近后仍保留长时间的历史依赖，仅按左右位置归类未必已经产生精确 Markov 链。

离散生成元还可以导出连续近似。若跳跃幅度为 $\delta$，左右速率为 $r_\pm(x)$，Taylor 展开给出

$$
\mathcal Lf\simeq
\delta(r_+-r_-)f'
+\frac{\delta^2}{2}(r_++r_-)f''.
$$

漂移来自两个方向事件率的差，扩散来自它们的和。扩散极限要求跳幅缩小、速率相应缩放，高阶项可忽略；不能只因主方程和 Fokker–Planck 形式相似，就认为任意少量大跳跃都可替换为高斯噪声。

## 4.7 Feynman–Kac：路径上的累计量怎样变成算子演化

最后考虑一个真正依赖整段历史的量：

$$
G(x,t)=\mathbb E_x\left[
e^{-\int_0^t V(X_s)ds}f(X_t)
\right].
$$

这里 $V$ 是给定的路径权重率，与前文的势函数 $U$ 无关；例如 $V=\alpha\,\mathbf1_R$ 可以记录粒子在区域 $R$ 中累计停留时间的 Laplace 变换。

先拆出最初一段 $h$，用 Markov 性接上剩余未来：

$$
G(x,t+h)\simeq\mathbb E_x[(1-hV(x))G(X_h,t)].
$$

展开条件平均，得到

$$
\boxed{\partial_tG=(\mathcal L-V)G,\qquad G(x,0)=f(x).}
$$

这就是这一形式的 Feynman–Kac 关系。它在适当的正则性、可积性和边界条件下，把路径累计量转成一个确定的算子问题。[Goodman 的随机微积分讲义](https://math.nyu.edu/faculty/goodman/teaching/StochCalc2004/notes/l6.pdf)用短时间条件平均发展了这种对应。

令 $f=1$、$V=\alpha\mathbf1_R$，对 $\alpha$ 在零点求导，即得

$$
-\left.\partial_\alpha G(x,t)\right|_{\alpha=0}
=\mathbb E_x\int_0^t\mathbf1_R(X_s)ds.
$$

所以相同的局部生成元还能计算一段历史里积累的驻留时间。由于 $(\mathcal L-V)1=-V$，加权演化通常不再守恒概率；它在这里用于统计路径权重，而不是直接替代原过程的 Markov 生成元。

## 4.8 四种视角怎样共同完成一个物理计算

给定漂移、噪声与边界，Langevin 描述先告诉我们局部运动及涨落注入；Fokker–Planck 描述追踪概率的输运与积累；路径统计回答约束下哪些历史较可能发生；生成元则把局部规则变成未来观测量与事件问题的方程。

OU 使这几种计算都能精确对照。双稳扩散说明，即使没有完整的闭式解，也能分别计算局部宽度、慢交换率、最优爬坡代价和首达时间，并检查它们是否彼此一致。实际使用时，先确定所问的是分布、关联、响应、去向概率还是等待时间，再选择最便于保留这些信息的表示。

**上一篇：**[Path Probabilities and Path Integrals](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/path_integrals/) · **[返回专题目录](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/)**。
