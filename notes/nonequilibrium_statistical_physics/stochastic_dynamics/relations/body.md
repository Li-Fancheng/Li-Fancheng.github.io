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

接下来沿三组关系各走一个来回：[Langevin 与 Fokker–Planck](#langevin-to-fp)之间，从随机增量得到概率输运，再从演化规律提取局部运动；[Langevin 与路径权重](#langevin-to-path)之间，从噪声构造历史，再从条件历史恢复运动方程；[Fokker–Planck 与路径积分](#fp-to-path)之间，从传播核拼接历史，再对历史求和恢复密度演化。随后用矩与关联、弛豫模、趋近平衡和双稳势逃逸，把这些联系变成处理具体物理问题的方法。

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

<span id="自然引出的工具后向生成元与到达概率"></span>

### 前向与后向生成元

前向方程追踪概率云；若问题改成“从不同起点出发，最终会发生什么”，让起点成为变量更自然。记

$$
\mathcal L f=b(x)\frac{df}{dx}+D\frac{d^2f}{dx^2},\qquad
\mathcal L^\dagger p=-\partial_x(bp)+D\partial_x^2p.
$$

$\mathcal L$ 是后向生成元：它给出从位置 $x$ 出发时，可观测量平均值的初始变化率；$\mathcal L^\dagger$ 则在相应边界条件下作用于密度。4.2 的分部积分正是两者的对应，不是给系统增添一条新动力学。

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

## 4.8 矩与关联：不求完整分布能知道什么 {#moments-correlations}

实验中，我们有时并不需要整个概率密度。要判断粒子平均向哪里移动、分布如何展宽，先计算均值和方差就够了。这些量既可以从单条 Langevin 轨道的随机方程求出，也可以从 Fokker–Planck 方程对密度积分得到。

对 $dx=b(x)dt+\sqrt{2D}\,dW_t$，4.2 中的公式对光滑函数 $f$ 给出

$$
\frac{d}{dt}\langle f(x)\rangle
=\left\langle b(x)\frac{df}{dx}+D\frac{d^2f}{dx^2}\right\rangle.
$$

取 $f=x$ 和 $f=x^2$，得到

$$
\frac{d\langle x\rangle}{dt}=\langle b(x)\rangle,\qquad
\frac{d\langle x^2\rangle}{dt}=2\langle xb(x)\rangle+2D.
$$

第一式中无偏噪声不直接改变均值；第二式中它不断提供 $2D$。减去均值平方的变化，还得到

$$
\frac{d\operatorname{Var}(x)}{dt}
=2\bigl(\langle xb(x)\rangle-\langle x\rangle\langle b(x)\rangle\bigr)+2D.
$$

这一式把分布宽度的变化写成两份收支：漂移把当前位置的差异放大或压缩，随机增量则持续注入方差。对 OU，$b=-\gamma x$，它恰好回到 $d\sigma^2/dt=-2\gamma\sigma^2+2D$。

### 为什么非线性系统通常不能只算均值

对线性漂移 $b=-\gamma x$，均值方程只含均值，二阶矩方程只含二阶矩，因而可以逐一求解。但若回复力含有非线性，例如 $b=-\gamma x-gx^3$，则

$$
\frac{d\langle x\rangle}{dt}=-\gamma\langle x\rangle-g\langle x^3\rangle.
$$

均值依赖三阶矩；三阶矩的方程又会引入更高阶矩。这是一串**矩方程层级**，并不会因为写出了第一条就自动闭合。将 $\langle x^3\rangle$ 换成某个均值和方差的函数，需要额外假设或近似；对于一般非线性动力学，这不是精确恒等式。因此，矩方程往往能揭示机制，但不总能代替完整分布。

### 两个时刻之间还需要什么信息

单时刻的均值和方差不会告诉我们一条轨迹多久失去记忆。假设系统以稳态密度 $\pi(x)$ 准备，两个可观测量的时间关联可由同一个传播核写成

$$
\langle A(x(0))B(x(t))\rangle_{\mathrm{st}}
=\int dx\,\pi(x)A(x)\int dy\,G(y,t\mid x,0)B(y).
$$

内层积分问：已知起点 $x$，过了 $t$ 以后 $B$ 的平均值是多少？外层再对稳态起点加权，并与初始的 $A$ 相乘。因此，关联函数不仅需要稳态分布，还需要连接两个时刻的转移规律。

对 OU，已知起点后的平均位置是 $xe^{-\gamma t}$，稳态位置方差为 $D/\gamma$，立刻得到 $\langle x(0)x(t)\rangle_{\mathrm{st}}=(D/\gamma)e^{-\gamma t}$。这正是第一篇从轨道解算出的关联，也说明单时刻直方图为何不能确定关联时间：它只提供上式的 $\pi$，没有给出传播核。


## 4.9 趋于稳态有多快：弛豫模与时间尺度 {#relaxation-modes}

OU 的平均位置按 $e^{-\gamma t}$ 返回原点，方差与稳态值的差却按 $e^{-2\gamma t}$ 衰减。为什么同一个系统会出现两个速率？因为一份概率分布可以偏离平衡的方式不止一种：它的中心可以偏移，宽度也可以偏大或偏小；不同形状的偏离不必以同样速度消失。

对 OU 的 Fokker–Planck 算符，令 $\mathcal L^\dagger p=\gamma\partial_x(xp)+D\partial_x^2p$，平衡密度为

$$
\pi(x)=\sqrt{\frac{\gamma}{2\pi D}}
\exp\left(-\frac{\gamma x^2}{2D}\right).
$$

直接代入可以验证以下三个模式：

$$
\mathcal L^\dagger\pi=0,\qquad
\mathcal L^\dagger(x\pi)=-\gamma x\pi,\qquad
\mathcal L^\dagger\!\left[\left(x^2-\frac D\gamma\right)\pi\right]
=-2\gamma\left(x^2-\frac D\gamma\right)\pi.
$$

第一个模式是保持不变的稳态；第二个在原点两侧符号相反，描述左右不对称的偏离；第三个描述中心与两侧之间的概率重分配。把适当的初始密度展开为这些模式及更高模式后，每一项按自己的速率衰减。这里 OU 模型的这些速率是 $n\gamma$，$n=0,1,2,\ldots$。具体哪一种衰减在观测中出现，还取决于初态和所测量的量是否对该模式有非零投影。

例如，测量 $\langle x\rangle$ 对左右反对称的偏离敏感，得到速率 $\gamma$；测量关于均值的宽度则反映二阶模式，得到 $2\gamma$。当相应模式存在时，最慢的非零衰减率给出长时间弛豫尺度；如果某个初态或观测量没有这一模式，它所呈现的衰减可以更快。

### 和虚时间 Schrödinger 方程的联系

平衡扩散还有一个有用的变换。设粒子在势能 $V(x)$ 中运动，$D=\mu k_{\mathrm B}T$，并令 $U(x)=V(x)/(k_{\mathrm B}T)$。对无净流的可归一化平衡分布，写

$$
p(x,t)=e^{-U(x)/2}\psi(x,t).
$$

代入 Fokker–Planck 方程，可以整理为 $\partial_t\psi=-\hat H\psi$，其中

$$
\hat H=-D\frac{d^2}{dx^2}
+D\left[\frac14\left(\frac{dU}{dx}\right)^2
-\frac12\frac{d^2U}{dx^2}\right].
$$

它在合适的边界条件下具有类似量子力学 Hamiltonian 的自伴形式，但这里演化变量是实时间中的概率，指数 $e^{-\hat Ht}$ 描述衰减，而不是实时间量子振幅的幺正演化。零模 $\psi_0\propto e^{-U/2}$ 对应平衡密度；其余本征值给出弛豫速率。谐势得到类似谐振子的算符，因此 OU 的等间隔弛豫率并非偶然。[UCSD 随机动力学讲义](https://courses.physics.ucsd.edu/2015/Fall/physics210b/LECTURES/STOCHASTIC.pdf)对 OU 的模式展开有详细计算。

这份变换依赖平衡势能结构和相应边界条件；有稳态环流的一般非平衡过程不能不加区分地套用同一个自伴构造。


## 4.10 从稳态到平衡：概率流与自由能 {#equilibrium-relaxation}

已知一份不随时间变化的分布，并不意味着系统一定处于平衡。想象圆环上有一群粒子沿同一方向不断绕行：每一小段的流入与流出相等，密度可以始终均匀，但运动记录仍有明确方向。要判断动力学是否满足平衡，需要连同概率流和正反运动一起看。

先看与热浴接触的过阻尼粒子。若 $F=-dV/dx$ 且 $D=\mu k_{\mathrm B}T$，平衡密度和概率流为

$$
\pi(x)=Z^{-1}e^{-V(x)/(k_{\mathrm B}T)},\qquad
J=\mu Fp-D\partial_xp
=-Dp\,\partial_x\ln\frac{p}{\pi}.
$$

这里假设 $\pi$ 可以归一化。平衡时 $p=\pi$，流在每一点为零：向低势能方向的漂移与热扩散恰好相抵。第二篇已经从概率流讨论过这一平衡；现在可以进一步问，如果初始分布不是 $\pi$，它为什么会向平衡靠近？

定义与平衡分布的相对熵，以及相应的非平衡自由能差：

$$
\mathcal D[p\|\pi]=\int p\ln\frac{p}{\pi}dx,\qquad
\mathcal F[p]-\mathcal F[\pi]=k_{\mathrm B}T\,\mathcal D[p\|\pi].
$$

对归一化密度，相对熵非负且在 $p=\pi$ 时为零。利用 $\partial_t p=-\partial_xJ$ 并分部积分，在实线远处边界项消失或反射边界 $J=0$ 的条件下，

$$
\frac{d\mathcal F}{dt}=-k_{\mathrm B}T\int\frac{J^2}{Dp}\,dx\le 0.
$$

因此，这个固定势能、单一平衡热浴模型中的自由能差不会自发增加。分布离平衡越远并不意味着每一点的密度都单调变化：概率可以先在某一处聚集再离开；单调的是这里定义的整体泛函。它把 Fokker–Planck 的概率输运与趋向平衡的热力学方向联系起来。

对于驱动圆环，稳态流可以非零。一般稳态密度本身仍不随时间改变，但正反方向的轨迹权重可以不同；第二篇已经以绕圈运动解释了这一点。上面的自由能耗散恒等式使用了**零流平衡分布**的具体关系，不能直接当作所有非平衡稳态的熵产生公式。


## 4.11 双稳势：热涨落怎样带来跨势垒运动 {#barrier-escape}

OU 描述粒子在一个稳定位置附近的涨落。若势能有两个势阱，粒子还可能在一个阱里徘徊很久，再因热涨落跨过中间的势垒进入另一个阱。单次轨迹上，这看起来像两种不同运动交替出现：阱内频繁的小幅涨落，以及少见的跨阱事件。

当势垒远高于 $k_{\mathrm B}T$ 时，阱内局部弛豫通常远快于逃逸。因而可以把长时间运动看作粒子大部分时间待在势阱附近、偶尔发生一次转换。我们先算平均等待时间，而不是仅凭势能图猜测它。

![双稳势中概率密度的阱内弛豫与跨阱转移](../figures/03_bistable_transport.svg)

*图｜无量纲势能为 $(x^2-1)^2/4$，取迁移率 $1$、扩散系数 $D=0.07$，初态集中在左阱附近。上图给出势能，中图比较不同时刻的概率密度，下图比较两阱所占概率；虚线是由 Kramers 速率构造的近似两态模型。图中 $U$ 表示势能，与正文的 $V$ 对应。*

先看中图：左阱内的峰较快形成；随后概率缓慢流向右阱。下图将左右阱的概率分别加总，可以更清楚地看到跨阱尺度比阱内弛豫慢。虚线是高势垒下的近似，不应将它当成任意参数下的精确结果。

### 从后向方程算平均逃逸时间

令左侧边界 $a$ 反射、右侧边界 $c$ 吸收；粒子位于 $x\in(a,c)$，第一次到达 $c$ 的时间为 $\tau$。固定起点的平均时间 $T(x)=\langle\tau\rangle_x$ 满足

$$
D\frac{d^2T}{dx^2}-\mu\frac{dV}{dx}\frac{dT}{dx}=-1,
\qquad \left.\frac{dT}{dx}\right|_{a}=0,\quad T(c)=0.
$$

一小步要先花费一段时间，然后从新的位置继续等待，这便是右侧 $-1$ 的来源。取 $U=V/(k_{\mathrm B}T)$、$D=\mu k_{\mathrm B}T$，用积分因子 $e^{-U}$ 求解，得到

$$
T(x)=\frac1D\int_x^c dy\,e^{U(y)}
\int_a^y dz\,e^{-U(z)}.
$$

内层积分衡量势阱中可以积累的概率权重，外层积分则强调通过高势能区域的困难。这个表达式适用于所述一维区间和边界；改变反射或吸收位置，就改变了要计算的逃逸实验。

### 高势垒下的 Kramers 结果

设粒子从左侧势阱底附近 $x_a$ 出发，势垒顶为 $x_b$；令 $\Delta V=V(x_b)-V(x_a)>0$，并记两处曲率

$$
\kappa_a=\left.\frac{d^2V}{dx^2}\right|_{x_a}>0,\qquad
\kappa_b=-\left.\frac{d^2V}{dx^2}\right|_{x_b}>0.
$$

若势垒远高于热能、两处可用二次势近似，反射端足够远且吸收端位于势垒另一侧足够远，则内层积分主要来自阱底，外层积分主要来自势垒顶。分别用 Gaussian 积分近似，可得过阻尼逃逸率

$$
k\simeq\frac{1}{T(x_a)}\simeq
\frac{\mu\sqrt{\kappa_a\kappa_b}}{2\pi}
\exp\left(-\frac{\Delta V}{k_{\mathrm B}T}\right).
$$

指数是 Arrhenius 因子：势垒每增加一个热能单位，逃逸率就大约乘上 $e^{-1}$。前因子同时依赖摩擦和势阱、势垒附近的形状，所以仅知道势垒高度还不足以算出完整速率。这个公式有高势垒、过阻尼以及指定边界等前提，不能直接推广到任意摩擦或任意初始准备。[MIT 的扩散课程](https://ocw.mit.edu/courses/18-366-random-walks-and-diffusion-fall-2006/resources/lec15/)也从势阱逃逸引出 Kramers 公式。

### 同一个指数，如何出现在路径权重里

从轨道看，向上爬势垒需要噪声持续对抗确定性的回复。对一条从 $x_a$ 到 $x_b$ 的光滑参考路径，其弱噪声代价为

$$
S_0[x]=\frac1{4D}\int\left(\dot x+\mu\frac{dV}{dx}\right)^2dt.
$$

利用 $D=\mu k_{\mathrm B}T$ 并展开平方，得到

$$
S_0[x]=\frac{\Delta V}{k_{\mathrm B}T}
+\frac1{4D}\int\left(\dot x-\mu\frac{dV}{dx}\right)^2dt.
$$

第二项非负，因此上坡所需的最小指数代价是 $\Delta V/(k_{\mathrm B}T)$。对应的参考上坡运动满足 $\dot x=+\mu\,dV/dx$，方向恰好与无噪声的下坡弛豫相反；从精确的阱底到精确势垒顶通常需要取长时间极限。概率方程算出逃逸时间，路径代价解释了其中的指数：二者都反映克服同一个势垒所需的热涨落。仅靠最小代价不能确定上面的曲率前因子，还需要计算附近路径的涨落与概率归一化。


### 先进入哪一个势阱：到达概率

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


## 4.12 多维系统与惯性：状态空间中发生了什么 {#multidimensional-phase-space}

一维位置模型使概率流只有向左和向右两个方向，但实际系统可能包含多个位置、速度或内部变量。首先仍要问：哪些变量放在一起，才能使当前状态足以预测下一步的条件分布？选定状态后，前面的三种描述都可以在这个更大的空间中使用。

令状态为 $\boldsymbol x=(x_1,\ldots,x_n)$，取 Itô 方程

$$
dx_i=b_i(\boldsymbol x)dt+\sum_\alpha B_{i\alpha}(\boldsymbol x)dW_\alpha,
\qquad D_{ij}(\boldsymbol x)=\frac12\sum_\alpha B_{i\alpha}B_{j\alpha}.
$$

不同分量的噪声可以相关；扩散矩阵 $D_{ij}$ 记录短时随机位移之间的协方差。相应的概率方程是

$$
\partial_t p=-\sum_i\partial_i(b_ip)+\sum_{i,j}\partial_i\partial_j(D_{ij}p),
\qquad J_i=b_ip-\sum_j\partial_j(D_{ij}p).
$$

因此 $\partial_t p=-\nabla\cdot\boldsymbol J$。一维中，稳态要求流为常数；在二维及更高维，稳态只要求流的散度为零，概率可以在状态空间绕圈而不改变密度。

一个具体例子是在二维位置上令

$$
\boldsymbol b=(-\gamma x-\Omega y,\;\Omega x-\gamma y),\qquad
D_{ij}=D\delta_{ij},\quad \gamma,D>0.
$$

回复部分 $-\gamma(x,y)$ 将概率拉向原点，旋转部分 $\Omega(-y,x)$ 使粒子绕原点运动。其稳态密度和流为

$$
\pi(x,y)=\frac{\gamma}{2\pi D}e^{-\gamma(x^2+y^2)/(2D)},\qquad
\boldsymbol J_{\mathrm{st}}=\Omega(-y,x)\pi(x,y).
$$

密度只取决于到原点的距离，圆周方向的流无散度，所以密度不变；改变 $\Omega$ 的符号，静态分布仍相同，环流方向却反转。这让“一张分布图不够确定动力学”成为可以看见的几何事实。对于时间反演下为偶的过阻尼位置变量，这个非零环流也表明该稳态不满足详细平衡。

当 $\boldsymbol D$ 为常数且正定时，光滑参考路径的弱噪声代价推广为

$$
S_0[\boldsymbol x]=\frac14\int_0^T
(\dot{\boldsymbol x}-\boldsymbol b)^{\mathsf T}
\boldsymbol D^{-1}(\dot{\boldsymbol x}-\boldsymbol b)\,dt.
$$

它衡量噪声需要补上的多方向运动；若某些方向涨落更强，同样大小的偏离在这些方向上的代价就更低。这里仍采用前面说明的光滑参考路径和离散化限制；退化或状态依赖的扩散不能直接套用这个形式。

### 保留速度时：从牛顿方程到相空间概率

现在回到第一篇的惯性 Langevin 方程，取 $v=\dot x$，热浴噪声满足 $\langle\eta(t)\eta(s)\rangle=2\zeta k_{\mathrm B}T\delta(t-s)$。完整状态 $(x,v)$ 的概率密度 $P(x,v,t)$ 满足 Kramers 方程：

$$
\partial_tP=-v\partial_xP+\frac1m\partial_v\!\left[\left(\frac{dV}{dx}+\zeta v\right)P\right]
+\frac{\zeta k_{\mathrm B}T}{m^2}\partial_v^2P.
$$

第一项是位置方向的输运；确定性力、阻力和热噪声改变的是速度。对于束缚势，平衡密度为 $P_{\mathrm{eq}}\propto\exp[-(mv^2/2+V)/(k_{\mathrm B}T)]$。它在位置方向仍有 $J_x=vP_{\mathrm{eq}}$：向右运动的粒子穿过同一位置，向左运动的粒子也穿过同一位置。

因此不能仅凭相空间概率流的某个分量非零就断定存在非平衡耗散。时间反演会将 $v$ 变号，需要区分 Hamiltonian 运动的可逆流与摩擦、噪声对应的不可逆部分。在速度充分迅速弛豫后，才能在合适条件下消去 $v$，恢复前面使用的位置扩散。这个例子说明，同一个物理粒子在不同时间尺度上可以有不同的有效状态空间。


## 4.13 条件历史与路径累计量 {#feynman-kac}

前面已经用概率密度、弛豫和势垒逃逸串起三种描述。接下来回到整段历史，看看对路径施加终点条件或为累计量赋权，会怎样改变我们使用的运动方程。

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

<span id="从终点密度走向整段历史feynmankac-公式"></span>

### 路径累计量与 Feynman–Kac 公式

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

## 4.14 这些对应何时需要修改 {#scope}

局部方程必须配合初始条件与边界。反射边界限制概率流，吸收边界允许概率流出；同一个内部漂移扩散方程会因此得到不同的传播核和路径集合。对吸收问题，核的积分可小于一，缺少的概率对应已被吸收的历史。仅检查离散时刻是否越界，可能漏掉两次采样之间的穿越。

最后，如果噪声幅度依赖位置，必须明确随机积分约定。例如 Itô 方程 $dx=b(x)dt+\sqrt{2D(x)}dW_t$ 对应

$$
\partial_t p=-\partial_x(bp)+\partial_x^2[D(x)p],\qquad
J=bp-\partial_x[D(x)p].
$$

这时不能简单把常数扩散流中的 $D$ 换成 $D(x)$；不同积分约定下相同符号的漂移项也可能代表不同过程。[第五篇的随机微积分](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/stochastic_processes/#ito-calculus)会解释这些规则。如果存在长记忆、跳跃或重尾等待，也需要重新选择完整状态或传播规律；第五篇将沿运动机制继续讨论。

回看整条主线：运动方程规定局部随机增量；增量的条件统计决定概率如何输运；条件转移的拼接则组织起整段历史。转换时明确状态、初态、边界和离散约定，三种描述便可以彼此检验，也各自保留最适合回答的问题。

**上一篇：**[Path Integrals：从终点分布到随机历史](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/path_integrals/) · **下一篇：**[Stochastic Processes：随机过程与随机微积分导览](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/stochastic_processes/)。
