前三篇从不同的观察方式出发，研究了同一种随机运动：Langevin 方程描述粒子怎样一步步运动，Fokker–Planck 方程描述一群粒子的概率怎样重新分配，路径积分则保留一整段历史的统计权重。这一章把它们放在一起，重新问一次：给定其中一种描述，如何得到另外两种？每一步又需要哪些信息？

我们先在最熟悉的情形下工作：位置已经是闭合的 Markov 状态，漂移为 $b(x)$，扩散系数 $D>0$ 为常数，空间为整条实轴。这样可以看清主要关系；边界、惯性和状态依赖噪声会在最后单独讨论。

## 4.1 从一小步看三种描述 {#one-step}

想象在时刻 $t$ 找到粒子位于 $x$。再过很短时间 $h$，它向哪里走尚不确定，但下一步的**条件分布**可以确定。对加性噪声的 Langevin 方程

$$
dx=b(x)dt+\sqrt{2D}\,dW_t,
$$

短时间增量的条件均值和方差分别为

$$
\langle\Delta x\mid x\rangle=b(x)h+o(h),\qquad
\langle(\Delta x)^2\mid x\rangle=2Dh+o(h).
$$

平均位移按 $h$ 缩放，典型随机位移却按 $\sqrt h$ 缩放。前者决定概率云的搬运，后者的平方与时间同阶，决定概率云的扩散。把这一步从所有可能的出发位置累加，得到密度演化；把连续各步的条件概率相乘，得到历史的权重。

| 描述 | 要保留的对象 | 回答的问题 |
| --- | --- | --- |
| Langevin | 当前状态与随机增量 | 一次运动如何随噪声发展？ |
| Fokker–Planck | 每个位置的概率密度及其流 | 初始概率云如何移动、展宽或弛豫？ |
| 路径积分 | 一系列中间位置及其联合权重 | 满足条件的历史有哪些，各占多大权重？ |

共同的桥梁是转移密度 $G(x,t\mid x_0,0)$。它给定起点以后终点的分布。对一个 Markov 模型，知道所有时间间隔的转移核，便可把多个时刻的条件分布依次拼接；而仅有某一个时刻的密度，通常不足以恢复时间顺序。

接下来沿三组关系各走一个来回：[Langevin 与 Fokker–Planck](#langevin-to-fp)之间，从随机增量得到概率输运，再从演化规律提取局部运动；[Langevin 与路径权重](#langevin-to-path)之间，从噪声构造历史，再从条件历史恢复运动方程；[Fokker–Planck 与路径积分](#fp-to-path)之间，从传播核拼接历史，再对历史求和恢复密度演化。沿途会遇到到达概率、固定终点后的有效动力学，以及累计观测量的生成函数。

## 4.2 从 Langevin 到 Fokker–Planck {#langevin-to-fp}

先从一条轨道计算可观测量 $f(x)$ 的变化。虽然 $\Delta x$ 很小，它的平方却不能忽略：$\langle(\Delta x)^2\rangle$ 与 $h$ 同阶。对光滑的 $f$ 展开并取条件平均，

$$
\langle f(x+\Delta x)-f(x)\mid x\rangle
=\left[b(x)\frac{df}{dx}+D\frac{d^2f}{dx^2}\right]h+o(h).
$$

一次导数来自平均位移，二次导数来自随机位移的方差。这正是 Itô 公式中二次变差的物理作用：位移的正负在平均时抵消，位移平方却持续积累。再对所有当前位置按 $p(x,t)$ 平均，得到

$$
\frac{d}{dt}\int f(x)p(x,t)dx
=\int\left(b\frac{df}{dx}+D\frac{d^2f}{dx^2}\right)p(x,t)dx.
$$

在无穷远处边界项消失的条件下，将导数通过分部积分移到密度上。由于上式对任意合适的测试函数成立，

$$
\boxed{\frac{\partial p}{\partial t}=-\frac{\partial}{\partial x}(bp)+D\frac{\partial^2p}{\partial x^2}}.
$$

它也可以写成概率守恒 $\partial_t p=-\partial_xJ$，其中 $J=bp-D\partial_xp$。第一项将概率随漂移搬运，第二项描述随机位移在密度不均匀时造成的净流。轨迹上的随机方程因此给出确定的密度方程；这里的确定性指给定初始分布和动力学以后，分布的演化确定，而不是每次实验的结果相同。

反过来看，如果已知上式对**所有**初始分布成立，且已规定常数 $D$、状态空间及边界行为，就读出了局部漂移 $b$ 和扩散强度 $D$，可以构造相应的 Langevin 过程。只有某一次测得的 $p(x,t)$ 则一般不够：不同动力学可能产生相同的一时刻分布，多维情形还可能有不改变密度的环流。要恢复动力学，需要更多关于时间演化和流的信息。

### 从概率方程反推随机运动

现在把方向倒过来。假设一维、非退化扩散的密度对任意初态都服从

$$
\partial_t p=-\partial_x[b(x)p]+\partial_x^2[D(x)p],\qquad D(x)>0.
$$

这里给定的是整个演化规律，而不只是一张概率密度的照片。取初态 $\delta(y-x)$，用其短时传播核 $G(y,h\mid x,0)$ 测量下一步，便可以从条件增量提取

$$
b(x)=\lim_{h\to0}\frac1h\int(y-x)G(y,h\mid x,0)dy,\qquad
D(x)=\lim_{h\to0}\frac1{2h}\int(y-x)^2G(y,h\mid x,0)dy.
$$

在满足通常的正则性和扩散过程假设时，可以选取 Itô 表示 $dx=b(x)dt+\sqrt{2D(x)}\,dW_t$。两式中的二阶矩没有减去均值平方，因为漂移引起的均值平方只有 $h^2$ 阶，在极限中消失。反向构造依赖连续扩散、状态空间和边界条件；仅给出两个矩而不限制跳跃，不能排除其他运动机制。高维时噪声矩阵的不同平方根还可给出相同扩散张量，因此驱动噪声的分解并非唯一。

为什么一张稳态图不够？在周长 $L$ 的圆环上，$dx=vdt+\sqrt{2D}\,dW_t\pmod L$ 对任意常数 $v$ 都有 $p_{\mathrm{st}}=1/L$，但稳态流为 $J_{\mathrm{st}}=v/L$。同一张均匀直方图可以对应顺时针、逆时针或无净流运动；记录相邻时刻的位移才能区分它们。静态占据反映粒子在哪里，短时条件转移则告诉我们它接下来怎样走。

### 自然引出的工具：后向生成元与到达概率

前向方程追踪概率云；若问题改成“从不同起点出发，最终会发生什么”，让起点成为变量更自然。记

$$
\mathcal L f=b(x)\frac{df}{dx}+D\frac{d^2f}{dx^2},\qquad
\mathcal L^\dagger p=-\partial_x(bp)+D\partial_x^2p.
$$

$\mathcal L$ 是后向生成元：它给出从位置 $x$ 出发时，可观测量平均值的初始变化率；$\mathcal L^\dagger$ 则在相应边界条件下作用于密度。4.2 的分部积分正是两者的对应，不是给系统增添一条新动力学。

例如粒子在区间 $(a,c)$ 中运动，两端都吸收。令 $q(x)$ 为从 $x$ 出发先到右端 $c$、而非先到左端 $a$ 的概率。先走一小步，再从新位置问同一个问题：$q(x)=\mathbb E_x[q(x_h)]$。展开后得到

$$
\mathcal Lq=0,\qquad q(a)=0,\qquad q(c)=1.
$$

这个 $q$ 称为到达概率或 committor。它和平均首次到达时间回答不同问题：前者问先去哪里，后者问要等多久。若力由势能产生、$b=-\mu\,dV/dx$，且 $D=\mu k_{\mathrm B}T$，积分上式得

$$
q(x)=\frac{\displaystyle\int_a^x e^{V(y)/(k_{\mathrm B}T)}dy}
{\displaystyle\int_a^c e^{V(y)/(k_{\mathrm B}T)}dy}.
$$

在双稳势中，可以把 $a,c$ 选在两个势阱附近，问一条从中间出发的轨道更可能先落入哪边。$q=1/2$ 是两种结果等可能的位置；对于一维对称势和对称端点，它位于对称中心。它由整段区间的势形共同决定，一般不能只看当前位置的力或某一张稳态密度图。这里给出的是有限区间上的精确到达概率，而不是跨势垒速率；讨论速率还需要时间尺度和初始准备。

## 4.3 从 Langevin 方程到路径权重 {#langevin-to-path}

回到 Langevin 方程，从每一步的 Gaussian 噪声直接构造历史权重。以前点 Euler 规则记 $\Delta x_k=x_{k+1}-x_k$，一步密度为

$$
G_h^{\mathrm E}(x_{k+1}\mid x_k)=\frac{1}{\sqrt{4\pi Dh}}
\exp\left[-\frac{[\Delta x_k-b(x_k)h]^2}{4Dh}\right].
$$

独立的噪声增量使各步密度相乘，指数中的代价便逐步相加：

$$
S_h[x]=\sum_{k=0}^{N-1}\frac{[\Delta x_k-b(x_k)h]^2}{4Dh},\qquad
\mathcal P_h=(4\pi Dh)^{-N/2}e^{-S_h[x]}.
$$

这里 $\mathcal P_h$ 是给定 $x_0$ 时关于后续 $N$ 个位置的密度。对一般漂移，这是一份**离散 Euler 模型的精确权重**，不等于原连续过程在有限 $h$ 下的精确核；在适当条件下，缩小步长才得到连续过程。路径积分是这套带归一化的时间切片构造的简写，而不是忽略归一化后的普通曲线求和。

对光滑参考曲线，$S_h$ 形式上对应 $\int_0^T[\dot x-b(x)]^2dt/(4D)$，表达实现该曲线所需的噪声代价。但典型 Brownian 路径不可微，不能把它直接作为普通光滑曲线代入；完整有限噪声的连续权重还取决于离散约定与相应因子。比较不同推导时，先确认它们使用的是同一时间切片和测度。

### 换一种时间切片：多出来的因子从哪里来

同一个连续过程可以用不同的离散规则逼近，但必须同时改变位置求值点和变量变换的因子。以常数 $D$ 为例，若用中点位置 $\bar x_k=(x_k+x_{k+1})/2$ 写每一步，

$$
\Delta x_k=b(\bar x_k)h+\sqrt{2Dh}\,\xi_k,
\qquad \xi_k\sim\mathcal N(0,1),
$$

从噪声换到位置时，$\xi_k$ 还通过中点漂移依赖 $x_{k+1}$。在步长足够小、该隐式更新局部可逆时，单步的 Jacobian 为

$$
\left|\frac{\partial\xi_k}{\partial x_{k+1}}\right|
=\frac{1}{\sqrt{2Dh}}\left|1-\frac h2\frac{db}{dx}(\bar x_k)\right|.
$$

前点规则中漂移只依赖已知的 $x_k$，对应的 Jacobian 是常数；中点规则却多出依赖位置的因子。对光滑参考曲线形式地取连续极限，它产生熟悉的中点作用量

$$
S_{\mathrm{mid}}[x]=\int_0^T\left[
\frac{(\dot x-b(x))^2}{4D}+\frac12\frac{db}{dx}\right]dt.
$$

例如 OU 的 $b=-\gamma x$，修正项为 $-\gamma T/2$。它不是新的物理力；它记录的是从噪声变量到路径变量的变换，以及中点求值方式。前点和中点的平方项也按各自规则求值，不能只拿一个光滑连续表达式比较而忽略 Brownian 路径的二次变差。连续路径积分的测度、归一化及离散约定必须一起指定。对状态依赖噪声，Jacobian 和漂移转换更复杂，不能直接套用这里的常数 $D$ 公式。

## 4.4 从路径权重反推 Langevin 方程 {#path-to-langevin}

若已知的是一整份路径概率，能否重新得到运动规律？先不要对整条路径作形式变分，而是取出某一时刻之后的一小段条件转移。对前点权重，固定当前 $x$ 并对未来其余位置积分，归一化后的一步增量满足

$$
\Delta x=b(x)h+\sqrt{2Dh}\,Z,\qquad Z\sim\mathcal N(0,1).
$$

读出条件均值 $bh$、条件方差 $2Dh$，再把各步独立的新噪声累加，连续极限就回到 Itô 方程 $dx=b(x)dt+\sqrt{2D}\,dW_t$。这里的独立性是给定当前状态后的局部规则，并不意味着轨迹的位置或位移都彼此独立：OU 的回复会使不同时段的位移相关。

反推也有范围：必须知道规范化的路径测度及其时间切片，不能仅凭一个未指定测度的形式作用量，或一条最小作用量曲线就断定原来的随机方程。若历史中有记忆，当前状态不足以确定下一步；需要保留历史或扩大状态空间。

### 固定终点后，轨道方程会怎样改变？

这是路径描述提供的一个新问题。先看所有从 $x_0$ 出发的历史，再只保留终点为 $x_f$ 的条件系综。给定时刻 $t$ 的当前位置 $x$，记接下来抵达指定终点的密度为

$$
H(x,t)=G(x_f,T\mid x,t).
$$

在 $t<T$ 且 $H>0$ 的区域，Bayes 条件化给出下一位置 $y$ 的转移密度：

$$
G_{\mathrm{cond}}(y,t+h\mid x,t)
=G(y,t+h\mid x,t)\frac{H(y,t+h)}{H(x,t)}.
$$

第二个因子提高那些更容易到达终点的下一步的权重。展开这个归一化的核，并用后向方程 $\partial_tH+\mathcal LH=0$，可读出条件过程的局部漂移：

$$
dx=\left[b(x)+2D\partial_x\ln H(x,t)\right]dt+\sqrt{2D}\,dW_t^{\mathrm{cond}}.
$$

额外的项不是原实验中新施加的力，而是终点筛选后样本的有效运动规律；条件过程使用自身的 Wiener 增量。对自由 Brownian，$H$ 是剩余时间 $T-t$ 的 Gaussian 核，因此

$$
dx=\frac{x_f-x}{T-t}dt+\sqrt{2D}\,dW_t^{\mathrm{cond}},\qquad 0\le t<T.
$$

离终点越远、剩余时间越短，条件漂移越强；单条条件轨道仍然粗糙并有随机涨落。端点处的表达式具有奇异形式，应通过 $t\uparrow T$ 的 bridge 极限理解。这就把第三篇“先按终点筛选路径”的做法转换成一条可生成同一条件系综的随机方程：已知原动力学及约束后，路径权重给出了新的局部条件动力学。这种条件转移的构造称为 Doob 变换，扩散桥的运动方程可参见[相关研究](https://www.numdam.org/item/10.1016/j.crma.2004.12.025.pdf)。

## 4.5 从概率传播到路径权重 {#fp-to-path}

Fokker–Planck 方程配上初态 $\delta(x-x_0)$，给出传播核 $G(x,h\mid x_0,0)$。对 Markov 过程，中间时刻 $s$ 的位置若未观测，就要对所有可能的 $y$ 积分：

$$
G(x,T\mid x_0,0)=\int dy\,G(x,T\mid y,s)G(y,s\mid x_0,0).
$$

这是概率的收支：每条从 $x_0$ 到 $x$ 的历史都在时刻 $s$ 经过某个位置，把这些位置的贡献相加便得到终点密度。继续插入更多时刻，记 $t_k=kh$、$Nh=T$，得到

$$
G(x_N,T\mid x_0,0)=\int\prod_{k=1}^{N-1}dx_k\prod_{k=0}^{N-1}G(x_{k+1},h\mid x_k,0).
$$

保留中间位置时，核的乘积就是这组位置的联合条件密度；对中间位置积分，便回到只问终点的传播核。再对初始位置按 $p_0(x_0)$ 平均，即得 $p(x,T)=\int G(x,T\mid x_0,0)p_0(x_0)dx_0$。这就是从路径统计回到概率云的方式。

## 4.6 从路径权重回到 Fokker–Planck {#path-to-fp}

这条联系也可以倒过来检验。先将路径乘积中最后一步保留，其余中间位置全部积分。一步传播满足

$$
p(x,t+h)=\int dy\,G_h^{\mathrm E}(x\mid y)p(y,t).
$$

对这个式子直接展开 Gaussian 核会遇到依赖 $y$ 的漂移；更清楚的做法是再次用测试函数。给定出发点 $y$，Euler 增量 $\Delta x=b(y)h+\sqrt{2Dh}\,Z$，$Z$ 是标准 Gaussian。于是

$$
\int f(x)p(x,t+h)dx
=\int p(y,t)\left[f(y)+h\bigl(b(y)f'(y)+Df''(y)\bigr)\right]dy+o(h).
$$

两侧相减、除以 $h$ 并取极限，再把导数从 $f$ 移到 $p$，便重新得到 4.2 的 Fokker–Planck 方程。这里没有额外假定新的概率流：漂移流和扩散流都来自路径权重的一步均值与方差。

在常漂移 $b=v$ 时，还可以直接完成整段路径的积分。每一步都是均值 $vh$、方差 $2Dh$ 的独立 Gaussian 增量；将 $N$ 步相加，终点的均值为 $x_0+vT$，方差为 $2DT$。因此路径权重对所有中间位置积分的结果是

$$
G(x,T\mid x_0,0)=\frac{1}{\sqrt{4\pi DT}}
\exp\left[-\frac{(x-x_0-vT)^2}{4DT}\right].
$$

它正是漂移扩散方程从点初态出发的解：每一步的偏移累加成云的平移，每一步的方差累加成云的展宽。

因此，路径积分与密度方程之间不是两套互不相干的规则。前者保留中间位置，后者将这些位置积分掉并追踪剩下的单时刻分布。如果想求已经到达终点的样本在时刻 $s$ 的分布，则保留 $s$ 处的位置并除以终点密度：

$$
p(y,s\mid x_0,0;x_f,T)
=\frac{G(x_f,T\mid y,s)G(y,s\mid x_0,0)}{G(x_f,T\mid x_0,0)}.
$$

分子是经过 $y$ 再抵达终点的权重；分母完成条件化。这也说明，未条件化的传播和固定终点后的路径系综虽然使用同一动力学，回答的问题并不相同。

### 从终点密度走向整段历史：Feynman–Kac 公式

有些问题既需要路径，又希望使用微分方程。例如粒子在一段时间里积累了多少位移面积、消耗了多少反应物，或在某个区域停留了多久？设沿途累计量为

$$
A_t=\int_0^t a(x(s))ds,\qquad
Z(\lambda,t)=\left\langle e^{\lambda A_t}\right\rangle.
$$

其中 $a(x)$ 可以是位置本身，也可以是区域指示函数。$Z$ 在存在时生成累计量的矩；例如 $\partial_\lambda Z|_{\lambda=0}=\langle A_t\rangle$。它依赖整段历史，普通的终点密度 $p(x,t)$ 不能单独提供它。

把权重 $e^{\lambda A_t}$ 放入路径和中，并同时保留终点，定义加权密度

$$
\rho_\lambda(x,t)=\left\langle e^{\lambda A_t}\delta(x-x(t))\right\rangle.
$$

一小段时间内，运动仍按原来的核传播，但每段历史额外获得约 $1+\lambda a(x)h$ 的权重。时间切片拼接并取极限，得到

$$
\partial_t\rho_\lambda=\mathcal L^\dagger\rho_\lambda+\lambda a(x)\rho_\lambda,
\qquad \rho_\lambda(x,0)=p_0(x),\qquad
Z(\lambda,t)=\int\rho_\lambda(x,t)dx.
$$

这就是前向的 Feynman–Kac 关系：路径上的累计权重变成密度方程中的局部项。最后一项会改变总权重，所以 $\rho_\lambda$ 一般不是归一化概率密度；当 $\lambda=0$，它才回到普通 Fokker–Planck 演化。吸收或反射问题仍需配合相应边界条件。

如果固定的是起点而不是终点，令 $u(x,t)=\mathbb E_x[e^{\lambda A_t}]$，从第一小步递推则得到后向形式

$$
\partial_tu=\mathcal Lu+\lambda a(x)u,\qquad u(x,0)=1.
$$

前向方程汇集各条历史到达终点的加权密度；后向方程从不同起点计算未来历史的权重。它们都源于同一套短时转移，差别在于把哪个位置作为变量。更一般的路径泛函与前、后向方程可参见 [Feynman–Kac 方程的研究](https://arxiv.org/abs/1806.03663)。

作为具体核对，取自由 Brownian $dx=\sqrt{2D}\,dW_t$，选择 $a(x)=x$。若从 $x$ 出发，路径面积的平均值是 $xt$，方差是 $2Dt^3/3$。由于面积是 Gaussian 变量，

$$
u(x,t)=\exp\left(\lambda xt+\frac{D\lambda^2t^3}{3}\right).
$$

把它代回 $\partial_tu=D\partial_x^2u+\lambda xu$，左右两边相等。对 $\lambda$ 求两次导数，就恢复面积的二阶矩；取对数后再求二次导数，则直接得到方差。这与第三篇从两时刻关联函数计算面积方差的结果一致，但现在还给出了整个 Gaussian 面积分布的生成函数。

## 4.7 用 OU 过程对照三种计算 {#ou-comparison}

把抽象联系落到同一个例子。对 $\gamma>0$，考虑从固定位置 $x_0$ 出发的 OU 过程：

$$
dx=-\gamma x\,dt+\sqrt{2D}\,dW_t.
$$

**轨道视角：**积分因子给出

$$
x(t)=x_0e^{-\gamma t}+\sqrt{2D}\int_0^t e^{-\gamma(t-s)}dW_s.
$$

初始偏离指数衰减，每一次新的随机位移也会被随后的回复逐渐削弱。噪声积分为 Gaussian，其均值为零，方差是各时段独立贡献之和。因此位置的均值和方差为

$$
m(t)=x_0e^{-\gamma t},\qquad
\sigma^2(t)=\frac D\gamma(1-e^{-2\gamma t}).
$$

**概率视角：**同一个过程的密度满足

$$
\partial_t p=\gamma\partial_x(xp)+D\partial_x^2p,\qquad p(x,0)=\delta(x-x_0).
$$

轨道解告诉我们位置服从 Gaussian，所以传播核是

$$
G(x,t\mid x_0,0)=\frac{1}{\sqrt{2\pi\sigma^2(t)}}
\exp\left[-\frac{(x-m(t))^2}{2\sigma^2(t)}\right].
$$

将它代入密度方程即可检验：均值以速率 $\gamma$ 回到原点，方差先增长、最终停在 $D/\gamma$。前者体现回复，后者体现回复对扩散展宽的限制。对一般初始密度，用这个核按初始位置积分；有限时间的结果不一定仍是 Gaussian。

也可以直接从密度方程读取这两个结果：分别乘以 $x$ 和 $x^2$ 积分，边界项消失时，

$$
\frac{d\langle x\rangle}{dt}=-\gamma\langle x\rangle,\qquad
\frac{d\sigma^2}{dt}=-2\gamma\sigma^2+2D.
$$

点初态给出 $\langle x\rangle(0)=x_0$、$\sigma^2(0)=0$；解这两个收支式，就得到轨道计算中的均值和方差。密度方程中的 $-2\gamma\sigma^2$ 是回复造成的压缩，$2D$ 是噪声持续注入的宽度。

下面重看第二篇的概率云：三幅图从同一个有有限宽度的初态出发，分别展示单独的回复、单独的扩散和二者共同作用。它与上面的点初态不同，但竞争机制是相同的。

![确定性回复、扩散和 OU 过程中的概率云](../figures/06_transport_mechanisms.svg)

*图｜初始均值为 $2$、方差为 $0.04$。左图只有回复，中图只有扩散，右图将回复与扩散结合；三图的纵轴标度各自独立。*

**历史视角：**将同一个 OU 核按时间相乘，就得到任意有限组中间位置的联合密度。对其中未观测的位置积分，可回到上面的终点核；若固定终点，再按上一节的条件密度保留中间位置，就能讨论 OU bridge。单时刻的稳态 Gaussian 不能单独告诉我们粒子多快探索它；时间序列中的关联，以及由转移核构造的路径统计，保留了这份信息。

## 4.8 这些对应何时需要修改 {#scope}

首先要选对状态。第一篇的惯性方程 $m\ddot x=-dV/dx-\zeta\dot x+\eta(t)$ 中，相同位置的粒子可能有不同速度，因此单独的 $x$ 通常不能闭合为 Markov 过程。完整状态是 $(x,v)$；相应的密度在相空间演化。只有在速度弛豫足够快、满足过阻尼近似时，才可以使用上面的位置方程和位置传播核。

例如在平衡白噪声热浴中，$\langle\eta(t)\eta(s)\rangle=2\zeta k_{\mathrm B}T\delta(t-s)$，相空间密度 $P(x,v,t)$ 满足 Kramers 方程

$$
\partial_tP=-v\partial_xP+\frac1m\partial_v\!\left[\left(\frac{dV}{dx}+\zeta v\right)P\right]
+\frac{\zeta k_{\mathrm B}T}{m^2}\partial_v^2P.
$$

位置流由速度产生，摩擦和热噪声则分别在速度方向带来漂移与扩散。相空间的转移核同样可以拼成历史权重；消去速度得到位置描述时，还必须检查时间尺度和边界条件。

其次，局部方程必须配合初始条件与边界。反射边界限制概率流，吸收边界允许概率流出；同一个内部漂移扩散方程会因此得到不同的传播核和路径集合。对吸收问题，核的积分可小于一，缺少的概率对应已被吸收的历史。仅检查离散时刻是否越界，可能漏掉两次采样之间的穿越。

最后，如果噪声幅度依赖位置，必须明确随机积分约定。例如 Itô 方程 $dx=b(x)dt+\sqrt{2D(x)}dW_t$ 对应

$$
\partial_t p=-\partial_x(bp)+\partial_x^2[D(x)p],\qquad
J=bp-\partial_x[D(x)p].
$$

这时不能简单把常数扩散流中的 $D$ 换成 $D(x)$；不同积分约定下相同符号的漂移项也可能代表不同过程。[第五篇的随机微积分](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/stochastic_processes/#ito-calculus)会解释这些规则。如果存在长记忆、跳跃或重尾等待，也需要重新选择完整状态或传播规律；第五篇将沿运动机制继续讨论。

回看整条主线：运动方程规定局部随机增量；增量的条件统计决定概率如何输运；条件转移的拼接则组织起整段历史。转换时明确状态、初态、边界和离散约定，三种描述便可以彼此检验，也各自保留最适合回答的问题。

**上一篇：**[Path Integrals：从终点分布到随机历史](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/path_integrals/) · **下一篇：**[Stochastic Processes：随机过程与随机微积分导览](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/stochastic_processes/)。
