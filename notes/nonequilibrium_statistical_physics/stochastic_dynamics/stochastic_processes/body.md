Brownian 的每一小步都重新抽取一个 Gaussian 位移；OU 又把过去的偏离逐渐拉回。如果运动改成一次次离散事件、有限速度的直行、偶尔的大跳，或很久不动的等待，前面的哪些结论还会成立？这一篇从运动规则出发，提供可独立查阅的模型和微积分导览。

Wiener、OU、Poisson 是过程；Itô/Stratonovich 是积分规则；Fokker–Planck 是密度演化方程；路径积分是组织历史统计的表示。先分清这些层次，再决定某个物理问题适合用哪一种工具。

## 4.1 先看运动机制，再看数学标签 {#process-map}

一个随机变量 $X$ 描述一次不确定读数；随机过程 $\{X_t\}_{t\ge0}$ 描述随时间变化的一族读数；将一次实现的所有读数连起来，才是一条样本路径。一个时刻的 $p(x,t)$ 不单独包含两个时刻的联合统计，更不单独规定整条路径。

比较模型时，至少有几条彼此独立的分类轴：路径连续还是跳跃，增量是否 Gaussian，二阶矩是否存在，更新是否独立，方向是否持续，是否有长等待，统计是否依赖准备后的年龄。它们不能排成一条由简单到复杂的单轴。

| 标准或典型模型 | 路径与运动机制 | Markov 状态 | 代表性统计 |
| --- | --- | --- | --- |
| Wiener / Brownian | 连续粗糙；独立 Gaussian 增量 | 位置本身 | MSD 正比于 $t$；位置非平稳 |
| OU | 连续粗糙；线性回复与 Gaussian 扰动 | 位置本身 | 可有平稳分布，时间关联非零 |
| 齐次 Poisson / 复合 Poisson | 计数或位置瞬时跳跃 | 计数或位置 | 指数等待；独立平稳增量 |
| 一维 run-and-tumble | 有限速度；随机翻向 | $(X,\sigma)$；单独 $X$ 通常不够 | 短时弹道、长时扩散 |
| 对称稳定 Lévy flight | 重尾增量；含瞬时跳跃 | 位置本身 | $0\lt \alpha\lt 2$ 时二阶矩发散 |
| 经典恒速 Lévy walk | 长飞行；位移与时长耦合 | 单独位置通常不够 | 有限传播范围，可超扩散 |
| 重尾等待 CTRW | 平台与跳跃；长时间陷获 | 单独位置通常不够 | 次扩散、年龄依赖 |
| fBm，$H\ne1/2$ | 连续 Gaussian；相关增量 | 一般非 Markov | 方差 $t^{2H}$，可持续或反持续 |

**Markov 的问题是当前完整状态是否已包含预测未来所需的信息。**惯性粒子只记录位置会漏掉速度；加入速度以后才可能闭合。OU 有时间相关却仍是 Markov，fBm 是 Gaussian 却一般非 Markov；Poisson 有独立增量却存在跳跃。Gaussian、独立、平稳、连续和 Markov 分别回答不同问题。

这里的“平稳”还须区分对象：位置本身的联合分布不随整体平移时间改变，才叫过程平稳；只要求等时长增量的分布相同，则是增量平稳。Wiener 属于后者，稳态准备的 OU 才属于前者。

## 4.2 Wiener：连续扩散怎样由小步组成 {#wiener}

要具体生成一条自由 Brownian 记录，选择步长 $h$，独立抽取 $Z_k\sim\mathcal N(0,1)$，令 $\Delta W_k=\sqrt h Z_k$，再累加这些增量。缩小步长时，这个构造趋向标准 Wiener 过程：

$$
W_0=0,\qquad W_t-W_s\sim\mathcal N(0,t-s)\quad(t>s),
$$

不相交区间上的增量独立，样本路径连续。对 $t\ge s$，写 $W_t=W_s+(W_t-W_s)$，后面的增量与 $W_s$ 独立且均值为零，因而

$$
\boxed{\mathbb E[W_tW_s]=\min(t,s),\qquad
\mathbb E[(W_t-W_s)^2]=|t-s|.}
$$

$\operatorname{Var}(W_t)=t$ 持续增大，所以位置并不平稳；任意等时长增量却有相同分布。乘上 $\sqrt{2D}$ 并加上初态，就得到[第一篇的自由位置过程](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/langevin/)。

典型增量按 $\sqrt h$ 缩小，差商却按 $h^{-1/2}$ 增大。连续不保证可微，Brownian 路径几乎处处没有普通切线；形式白噪声 $\xi=\dot W$ 必须通过时间积分理解。这也解释了下一节为什么需要改变微积分规则。

## 4.3 Itô：随机增量的平方为什么不能丢掉 {#ito-calculus}

Itô 不是另一种运动模型，而是用什么信息定义随机积分、怎样对该积分运算的规则。最短的例子是平方：每一步有

$$
W_{k+1}^2-W_k^2=2W_k\Delta W_k+(\Delta W_k)^2.
$$

对整段时间相加，最后一项不会像光滑运动那样消失：在分割变细时，$\sum_k(\Delta W_k)^2\to t$（依概率收敛）。所以

$$
\boxed{d(W_t^2)=2W_t\,dW_t+dt.}
$$

这是二次变差的作用。把固定时长切成 $N$ 份，每个随机位移平方约为 $1/N$，总共 $N$ 份，贡献仍有限；光滑位移的平方则约为 $1/N^2$，总和趋于零。

考虑 Itô 扩散

$$
dX_t=b(X_t,t)\,dt+\sigma(X_t,t)\,dW_t.
$$

先对 $f(X,t)$ 作 Taylor 展开。一次项保留 $f_x\Delta X$，二次项 $f_{xx}(\Delta X)^2/2$ 的主要部分为 $f_{xx}\sigma^2\Delta t/2$；$\Delta t\Delta W$ 和更高阶项在适当正则条件下消失。因此

$$
\boxed{
df(X_t,t)=
\left[f_t+b f_x+\frac{\sigma^2}{2}f_{xx}\right]dt
+\sigma f_x\,dW_t.}
$$

这就是 Itô 公式。额外的二阶导数项把随机轨道的粗糙性转化成了可观测量的系统性变化。

也可以从测量函数的曲率理解这个项。把一个点左右对称地随机移动一点，位置的平均值没有改变；但若探测器读数是位置平方，两边移动都会让读数增大。对于一般弯曲的 $f(x)$，正负位移对读数的影响也不完全抵消，留下的量由局部曲率 $f_{xx}$ 决定。Itô 修正记录了这种“零平均扰动经过非线性测量后产生平均变化”的效应。

例如自由扩散 $dX=\sqrt{2D}\,dW$ 的平均位置不动，但

$$
d(X^2)=2D\,dt+2\sqrt{2D}X\,dW,
\qquad \frac{d}{dt}\langle X^2\rangle=2D.
$$

噪声平均为零，并不意味着它对二阶量没有平均效应。扩散展宽正来自独立随机增量的平方累积。

对带回复作用的运动，这项持续注入的方差还要与确定性收缩竞争。如果遗漏二次变差，OU 模型的方差方程就只剩下衰减，错误地预言热浴中的粒子最终完全停止涨落。

### 随机积分约定与乘性噪声

对于普通积分，区间左端或中点取值在极限中通常没有差别。对于随机积分，被积函数也会随同一个噪声变化，这种差别能够留下有限贡献。

物理上的关键在于：这一小步的随机扰动，是作用在扰动到来之前的状态上，还是用已经被这一步改变过的状态来计算它的作用强度？如果噪声振幅处处相同，这个差别无关紧要；如果振幅依赖位置，同一次扰动既改变位置，又改变它自身被计入的幅度，相关性就可能累积成可观测的漂移。

Itô 积分使用每步左端的状态：

$$
\int_0^t X_s\,dW_s
=\lim\sum_n X_{t_n}(W_{t_{n+1}}-W_{t_n}).
$$

在适应性和平方可积条件下，$X_{t_n}$ 不使用未来的增量，积分均值为零，并满足 Itô 等距

$$
\mathbb E\left[\left(\int_0^t H_s\,dW_s\right)^2\right]
=\mathbb E\int_0^t H_s^2\,ds.
$$

取 $X=W$，离散恒等式 $W_{n+1}^2-W_n^2=2W_n\Delta W_n+(\Delta W_n)^2$ 给出

$$
\int_0^t W_s\,dW_s=\frac12(W_t^2-t).
$$

Stratonovich 积分以对称取值定义，记作 $\circ dW$，对应

$$
\int_0^t W_s\circ dW_s=\frac12W_t^2.
$$

它保留普通形式的链式法则，但积分中的噪声与区间内状态已有关联，不能直接断言其均值为零。

这里可以看到“同一个物理过程有不同写法”与“换了一套物理规律”的区别。左端求值没有把本步噪声反馈进本步系数，中点求值则包含了部分反馈。若想让二者描述相同的轨迹统计，就必须在平均漂移里补偿这项差别。这正是下面转换公式的作用。

对于一维乘性噪声，相同过程的两种表达为

$$
dX=b_{\mathrm S}(X)\,dt+\sigma(X)\circ dW
$$

和

$$
dX=\left[b_{\mathrm S}(X)+\frac12\sigma(X)\sigma'(X)\right]dt
+\sigma(X)\,dW.
$$

选择约定时要同时转换漂移。若 $\sigma$ 为常数，两者的方程相同；若 $\sigma$ 依赖状态，保持同一漂移再改变积分约定就改变了模型。由平滑有色外力取某类短关联时间极限时常得到 Stratonovich 形式；由离散跳跃过程的条件增量得到扩散近似时常自然得到 Itô 形式。具体粗粒化和同时消去哪些变量，决定正确的极限，不能仅靠“物理噪声”四个字选择。



因此，一维转换可以概括为 $b_{\mathrm I}=b_{\mathrm S}+\sigma\sigma'/2$。使用 Itô 漂移时，一般密度方程是

$$
\frac{\partial p(x,t)}{\partial t}
=-\frac{\partial}{\partial x}[b_{\mathrm I}(x)p(x,t)]
+\frac12\frac{\partial^2}{\partial x^2}[\sigma^2(x)p(x,t)].
$$

导数作用于整个 $\sigma^2p$，不能将常数扩散流机械替换为 $-D(x)p_x$。当 $\sigma=\sqrt{2D}$ 为常数时，转换修正为零，二阶项就是第二篇的 $Dp_{xx}$。以上规则使同一过程的不同写法相容；如果保持同一漂移却改换积分约定，就改变了模型。


## 4.4 OU：有回复、有关联，为什么仍然 Markov {#ou-model}

一个稳定位置附近的涨落可以持续很久，但这并不要求运动记住整段过去。OU 的模型卡是

$$
dX_t=-\gamma(X_t-\bar x)dt+\sqrt{2D}\,dW_t,
\qquad\gamma>0.
$$

[第一篇](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/langevin/)已求过解，[第二篇](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/fokker_planck/)已解释概率云。这里只读出分类所需的性质：稳态准备时，均值为 $\bar x$、方差为 $D/\gamma$，连通关联为 $(D/\gamma)e^{-\gamma|t-s|}$；固定初态的有限时间过程则还未平稳。

从当前时刻向后看一小段，

$$
X_{t+h}-\bar x=e^{-\gamma h}(X_t-\bar x)+Z_{t,h},
\qquad
\operatorname{Var}(Z_{t,h})=\frac D\gamma(1-e^{-2\gamma h}),
$$

其中 $Z_{t,h}$ 是与此前历史独立的 Gaussian 新增量。知道当前值以后，不必再给更早记录，就能预测下一段条件分布。这正是 OU 的 Markov 性；当前位置与此前位置仍然相关，二者并不矛盾。

### 将 OU 当作有色噪声，又会怎样

若环境扰动本身要经过时间 $\tau_c$ 才遗忘，可以令

$$
d\eta_t=-\frac{\eta_t}{\tau_c}dt+\sigma_\eta dW_t,
\qquad \dot X_t=b(X_t)+\eta_t.
$$

稳态环境的关联为 $(\sigma_\eta^2\tau_c/2)e^{-|t-s|/\tau_c}$，不是 Delta 函数。联合状态 $(X_t,\eta_t)$ 可闭合为 Markov 过程；只知道 $X_t$，一般还不知道当前的驱动偏向，因而会显出记忆。这个例子回到第一篇的变量选择：隐藏了慢环境变量，常会把局部规律变成有记忆的约化运动。

## 4.5 Poisson 与跳跃：随机性也可以是一件件事件 {#poisson-jumps}

如果一次反应或跳跃在短时间 $h$ 内发生的概率为 $rh+o(h)$，多次事件的概率为 $o(h)$，而不相交时间窗的事件独立，就得到齐次 Poisson 计数。令 $P_n(t)=\Pr(N_t=n)$，流入第 $n$ 个计数的速率来自 $n-1$，流出来自下一次事件：

$$
\frac{dP_n(t)}{dt}=rP_{n-1}(t)-rP_n(t),
\qquad P_n(0)=\mathbf1_{n=0}.
$$

先解 $P_0=e^{-rt}$，再逐个积分，得到

$$
\boxed{\Pr(N_t=n)=e^{-rt}\frac{(rt)^n}{n!},
\qquad\mathbb EN_t=\operatorname{Var}(N_t)=rt.}
$$

从 $P_0(t)$ 还可读出首次等待的密度 $\psi(\tau)=re^{-r\tau}$。等了多久并不改变接下来的事件率，这种无记忆等待允许只用当前计数预测未来。$N_t$ 本身随时间增加，非平稳；增量的分布却只依赖时间窗长度。

### 计数发生一次，与位置跳一步是不同对象

若第 $k$ 次事件带来独立位移 $Y_k$，且跳幅与事件时间独立，则

$$
X_t=X_0+\sum_{k=1}^{N_t}Y_k
$$

是复合 Poisson 过程。对有限二阶矩的跳幅，$\mathbb E(X_t-X_0)=rt\mathbb EY$，$\operatorname{Var}(X_t-X_0)=rt\mathbb E(Y^2)$。后一式同时包含跳幅随机和事件数随机；即使每次同向跳固定距离，计数起伏仍会造成位置方差。

对离散状态 $i\to j$ 的事件率 $w_{i\to j}$，概率收支推广为主方程

$$
\boxed{\frac{dp_i(t)}{dt}
=\sum_{j\ne i}\left[w_{j\to i}p_j(t)-w_{i\to j}p_i(t)\right].}
$$

每一项都是流入减流出，与第二篇的守恒思路相同。连续 Gaussian 小步不是唯一机制：有限的大跳需要保持非局部转移；只有跳幅变小、事件加快并有适当矩条件时，才可取漂移扩散极限。不能把重尾大跳无条件截成局部二阶 Fokker–Planck 方程。

## 4.6 持续游走：方向不必每一瞬间重新选择 {#persistent-walk}

对一维 run-and-tumble 模型，粒子以固定速率 $v$ 向左或右运动，方向 $\sigma_t=\pm1$ 以速率 $\alpha$ **翻转**。从左右等概率的方向初态出发，均值位移为零，路径却先是一段段直行：

$$
\dot X_t=v\sigma_t.
$$

若到时刻 $t$ 翻转 $N_t$ 次，$\sigma_t\sigma_0=(-1)^{N_t}$。Poisson 求和给出 $\mathbb E[(-1)^{N_t}]=e^{-2\alpha t}$，因此方向的关联时间为 $1/(2\alpha)$。这里的 $\alpha$ 指翻转率；若文献规定“重新随机选方向”的事件率，数值中的因子 2 会随约定改变。

像第一篇积分速度关联一样，

$$
\begin{aligned}
\langle[X_t-X_0]^2\rangle
&=2v^2\int_0^t(t-s)e^{-2\alpha s}ds\\
&=\boxed{\frac{v^2}{\alpha}t
-\frac{v^2}{2\alpha^2}(1-e^{-2\alpha t})}.
\end{aligned}
$$

短时为 $v^2t^2$，长时为 $2D_{\mathrm{eff}}t$，其中 $D_{\mathrm{eff}}=v^2/(2\alpha)$。这与惯性粒子有相同形式的跨越：把[图 1a](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/langevin/#inertial-msd)中的记忆时间换成 $1/(2\alpha)$ 即可。但一个来自热速度弛豫，一个来自随机换向，长时的扩散相同并不固定短时机制。

单独的位置一般不足以预测下一小步，因为还需要知道向哪边走；联合状态 $(X,\sigma)$ 才闭合。最简模拟也不必每步抽 Gaussian 数：独立抽指数等待 $-\log U/\alpha$，沿当前方向走完这一段，再翻转方向。持续游走和扩散近似的联系可对照 [MIT 随机游走讲义](https://ocw.mit.edu/courses/18-366-random-walks-and-diffusion-fall-2006/pages/lecture-notes/)。

## 4.7 Lévy flight 与 walk：长跳不等于长时间飞行 {#levy-flight-walk}

如果运动中出现很长的一次位移，需要先问它怎样完成：在模型里瞬间跳到远处，还是以有限速度走了很久？两个机制会给出不同的路径和可测统计量。

### Lévy flight：重尾增量和瞬时跳跃

选择典型的一维对称稳定 Lévy 过程，其独立平稳增量满足

$$
\boxed{\mathbb E\,e^{ik(X_t-X_0)}=e^{-Kt|k|^\alpha},
\qquad K>0,\quad0\lt \alpha\lt 2.}
$$

$K$ 的单位是长度的 $\alpha$ 次方除以时间。特征函数依赖 $k(Kt)^{1/\alpha}$，故分位宽度等典型尺度正比于 $t^{1/\alpha}$。但这不是均方位移：稳定分布尾部使二阶矩发散。有限样本当然能算出一个平方平均，却不是一个有限理论 MSD 的可靠估计；可以改用分位宽度或 $0\lt q\lt \alpha$ 的绝对矩。

这个过程是非 Gaussian、含跳跃的 Markov 过程，有独立增量；它没有规定最大传播速度。不能因平均图像像“更快扩散”，就继续使用 Brownian 的局部小步解释。

### Lévy walk：位移与飞行时长耦合

经典恒速模型每次选择方向和飞行时长 $\tau$，以速度大小 $v$ 持续前进，飞行结束再更新。若 $\psi(\tau)\sim\tau^{-1-\alpha}$，会出现长直行；位移大小 $|\ell|=v\tau$ 与耗时相连，因此任意有限时间都有传播范围 $|X_t-X_0|\le vt$。

在常见的一维恒速、独立对称方向更新模型中，当 $1\lt \alpha\lt 2$，长时可有 $\langle\Delta X^2\rangle\propto t^{3-\alpha}$。这要求上述模型和指数区间；不是所有 Lévy walks 的定义。$\alpha$ 改变或加入停顿、不同速度和不同初始年龄，都需要重新检查统计结论。

![稳定 Lévy flight 与恒速 Lévy walk 的空间时间记录](../figures/11_levy_flight_walk.svg)

*图 4a｜两图取尾指数 $\alpha=1.3$、$X_0=0$、观察时长 20。左图为 $K=1$ 的对称稳定 flight，在网格上采样，以阶梯连线标出跳跃；它没有画出每格内部更小的跳跃，阶梯不代表额外陷获。右图从新飞行开始，速度大小 $v=1$，飞行时长取尺度 $\tau_0=0.3$ 的 Lomax 分布、方向独立等概率；虚线为有限速度前沿。这里展示路径机制，不用 flight 的样本 MSD 推断不存在的理论二阶矩。*

flight 的长距离来自重尾跳幅，walk 的长距离来自持续飞行，两者的传播前沿与路径连续性不同。[Zaburdaev、Denisov 与 Klafter 的综述](https://doi.org/10.1103/RevModPhys.87.483)系统区分了这些模型。

## 4.8 CTRW：有时次扩散是因为很久没有迈出一步 {#ctrw}

Lévy walk 的长时段仍在走；重尾等待 CTRW 却可以很久不动。每次先等一个独立时间 $\tau_k$，再跳一个独立、零均值、有限方差 $\sigma_Y^2$ 的位移 $Y_k$，两组随机量彼此独立。从时刻零开始一次新的等待，令完成的步数为 $N(t)$：

$$
X_t-X_0=\sum_{k=1}^{N(t)}Y_k,
\qquad \psi(\tau)\sim c\tau^{-1-\beta},\quad0\lt \beta\lt 1.
$$

平均等待发散，但每步并不一定很小。真正变慢的是时钟：实际步数随物理时间增长得比线性更慢。

### 从等待分布算出步数标度

令 $m(t)=\mathbb EN(t)$。至少走完 $n$ 步，需要前 $n$ 次等待之和不超过 $t$。将这些概率对 $n\ge1$ 求和，并作时间 Laplace 变换，独立等待的卷积变为乘积：

$$
\widetilde m(s)
=\frac1s\sum_{n=1}^\infty[\widetilde\psi(s)]^n
=\frac{\widetilde\psi(s)}{s[1-\widetilde\psi(s)]}.
$$

对上述重尾，低频行为为 $1-\widetilde\psi(s)\sim a_\beta s^\beta$，$a_\beta>0$ 的单位是时间的 $\beta$ 次方。于是 $\widetilde m(s)\sim1/(a_\beta s^{1+\beta})$，利用 $t^\beta$ 的 Laplace 变换，得到长时

$$
\boxed{m(t)\sim\frac{t^\beta}{a_\beta\Gamma(1+\beta)},
\qquad\langle[X_t-X_0]^2\rangle
=\sigma_Y^2m(t)\propto t^\beta.}
$$

最后一步由给定 $N$ 后的独立零均值跳幅方差相加得到。若跳幅有偏、方差无限，或与等待耦合，这个推导会改变。这里不需要引入分数阶 Fokker–Planck；少走的步数已经解释了次扩散。[MIT 的 CTRW 讲义](https://ocw.mit.edu/courses/18-366-random-walks-and-diffusion-fall-2006/f3713d068fc4f5de471d599badfa9691_lec23_chris_06.pdf)从同一更新结构展开。

### 已经等了多久，可能影响接下来多久会动

指数等待没有年龄；重尾等待中，一个已持续很久的平台，接下来再继续等待的统计会不同。系统从制备开始已经过 $t_a$，未来观察窗长为 $T$ 时，平均新增步数是 $m(t_a+T)-m(t_a)$，依赖年龄 $t_a$。在上述长时区间，相应位移平方平均正比于 $(t_a+T)^\beta-t_a^\beta$。

因此只记录位置通常不是 Markov。若再保留当前等待的年龄，可以建立扩充状态的更新模型。位置平台、等待尾部和观察年龄，提供了比一个 MSD 指数更多的物理信息。

## 4.9 fBm：同样的幂律，也可能来自相关增量 {#fractional-brownian}

次扩散并不一定需要停滞。标准化的 fractional Brownian motion $B_H(t)$ 是零均值 Gaussian 过程，$B_H(0)=0$，具有连续路径与协方差

$$
\boxed{\mathbb E[B_H(t)B_H(s)]
=\frac12\left[t^{2H}+s^{2H}-|t-s|^{2H}\right],
\qquad0\lt H\lt 1.}
$$

它给出 $\operatorname{Var}(B_H(t))=t^{2H}$；位置非平稳，增量平稳。恢复物理单位可用 $X_t=X_0+\sqrt{2K_H}B_H(t)$，其中 $K_H$ 的单位是长度平方除以时间的 $2H$ 次方。$H=1/2$ 恰好恢复 Wiener。

等时长增量 $\Delta B_k=B_H((k+1)h)-B_H(kh)$ 的相关性，可以直接从协方差相减得到。相隔 $n\ge1$ 步时，标准化关联为

$$
\rho(n)=\frac12\left[(n+1)^{2H}+|n-1|^{2H}-2n^{2H}\right].
$$

$H>1/2$ 时正相关，前一步的方向较易延续；$H\lt 1/2$ 时负相关，较易反向抵消。这是持续或反持续的 Gaussian 增量机制，不是随机等待造成的少走几步。除 $H=1/2$ 外，它一般非 Markov，也不能直接照搬普通 Wiener 驱动的 Itô 规则或用普通局部扩散方程代表全部历史。

### 相同次扩散指数的两种运动

取 $\beta=2H=0.6$。重尾等待 CTRW 和 $H=0.3$ 的 fBm 都可有系综 MSD 正比于 $t^{0.6}$，却有不同的路径与时间统计：

![相同次扩散指数下 CTRW 与 fBm 的路径、MSD 与增量相关](../figures/12_ctrw_fbm.svg)

*图 4b｜CTRW 从新等待开始，$\beta=0.6$、Lomax 等待尺度 $\tau_0=0.05$、独立对称单位跳幅；fBm 取 $H=0.3$，用相关 Gaussian 增量的协方差矩阵生成。上方为各一条样本；下方使用 1024 条实现，分别比较归一化 MSD 的趋势与等间隔增量的关联。MSD 各自除以终值以比较形状，CTRW 在有限时间仍有更新瞬态；平台与反持续关联已能区分机制。*

这个比较不靠两条样本恰好长得怎样来判定模型。CTRW 的长等待形成平台，离散更新时钟随年龄改变；fBm 的连续 Gaussian 样本没有这套等待事件，反持续性直接出现在增量协方差中。

## 4.10 怎样从记录中诊断运动，而非只拟合一个指数 {#diagnosis}

有限采样会把很快的运动压成一个跳跃，也会把微小位移藏成平台。先明确时间分辨率、定位误差和事件识别规则，再比较多种统计：

| 观察对象 | 能帮助识别什么 | 不能单独断言什么 |
| --- | --- | --- |
| 轨迹形状 | 长直行、平台、跳跃或连续粗糙 | 有限分辨率不能证明真正瞬时跳跃 |
| MSD 与有效斜率 | 展宽随时间增长的快慢 | 相同指数不固定机制；稳定 flight 无有限二阶矩 |
| 位移尾部、分位宽度 | Gaussian、重尾、有限速度前沿 | 有限样本尾部拟合可能很不稳定 |
| 速度或增量关联 | 持续、反持续与记忆时间 | 有关联不等于非 Markov |
| 等待和停留时间 | 陷获与年龄依赖 | 需要先定义一次事件或一个区域 |
| 条件转移 | 当前状态是否足够 | 隐藏速度或内部态可能制造表观记忆 |
| 时间平均与系综平均 | 准备、年龄和测量时长的作用 | 单条短记录不足以确立遍历性质 |

一个有判别力的补充是固定延迟 $\Delta$ 的时间平均位移平方

$$
\overline{\delta^2}(\Delta)
=\frac1{\mathcal T-\Delta}\int_0^{\mathcal T-\Delta}
[X_{t+\Delta}-X_t]^2dt.
$$

对标准 fBm，其期望为 $\Delta^{2H}$，不依赖测量总时长 $\mathcal T$，因为增量平稳。对从新等待开始、无偏且有限跳幅方差的重尾 CTRW，在 $\mathcal T\gg\Delta$ 的相应长时极限，期望却可按 $\Delta\mathcal T^{\beta-1}$ 缩小：记录越长，后期长平台占的作用越明显。两者的系综指数可以相同，时间平均行为却不同；这一对照可参见 [Deng 与 Barkai 的原始研究](https://arxiv.org/abs/0809.2430)。

图 4b 的复现脚本还输出两种过程在前后等长观察窗内的平均新增步数或位移平方，以及改变测量总时长后的时间平均量。可以先对照平台与增量关联，再看年龄和时间平均；不应只对两条有限轨迹各拟合一个斜率便宣告模型相同。

在图 4b 的同一组无量纲记录里，取延迟 $\Delta=0.2$，1024 条实现给出下表。比较同一行随测量时长的变化，而非将不同模型的绝对幅度相互等同：

| 过程 | 测量时长 $\mathcal T=10$ | 测量时长 $\mathcal T=20$ |
| --- | --- | --- |
| CTRW 的平均时间平均量 | 0.265 | 0.199 |
| fBm 的平均时间平均量 | 0.381 | 0.380 |

fBm 理论值为 $0.2^{0.6}\simeq0.381$；CTRW 则随记录延长而减小。有限样本值用于区分机制，不能替代对任意时长的渐近证明。完整参数、代码与输出见[复现脚本](../generate_intro_figures.py)和[数值记录](../intro_numerical_checks.json)。

Wiener、OU、Poisson、持续游走、Lévy、CTRW 和 fBm 分别改变了回复、事件更新、方向持续、跳幅尾部、等待时钟或增量相关。先判断这些物理机制，再选择过程模型和数学工具。


## 选读附录 A：乘性涨落的平均增长与典型增长

一个能直接看出差异的例子是几何 Brownian 运动：

$$
dX=aX\,dt+sX\,dW,\qquad X_0>0.
$$

对 $\log X$ 用 Itô 公式，得到

$$
X_t=X_0\exp\left[(a-s^2/2)t+sW_t\right].
$$

其均值为 $X_0e^{at}$，而中位数为 $X_0e^{(a-s^2/2)t}$。少量增长很大的轨迹可以支撑均值，即使典型轨迹正在减小。这说明对随机运动，“平均增长”和“典型增长”也需要分开。[Pavliotis 的教材](https://www.ma.imperial.ac.uk/~pavl/PavliotisBook.pdf)第 3 章给出了随机积分和转换公式的系统处理。

一个直观原因是，乘性起伏的增减不对称地影响长期结果。一个量先乘以 $1+\varepsilon$，再乘以 $1-\varepsilon$，最后得到的是原值的 $1-\varepsilon^2$ 倍；相同幅度的一次相对增长和一次相对下降，不能相互复原。连续极限中，这种二阶积累正对应对数增长率里的修正。少数极大的结果又会抬高算术平均，所以平均轨迹和多数样本的经历可以明显不同。

## 选读附录 B：多变量状态与密度方程的核对

进一步把所保留的状态统记为 $\boldsymbol Z$，许多白噪声驱动的 Markov 扩散都可以写成

$$
\boxed{d\boldsymbol Z_t=\boldsymbol a(\boldsymbol Z_t,t)dt
+B(\boldsymbol Z_t,t)d\boldsymbol W_t.}
$$

这是 Itô 随机微分方程的一般形式。漂移 $\boldsymbol a$ 规定给定当前状态时的平均增量，矩阵 $B$ 规定各变量的随机增量如何由同一组独立 Wiener 驱动产生。小时间窗内的条件协方差为 $BB^{\mathsf T}h+o(h)$，所以噪声怎样耦合不同自由度也包含在这个表达式中。

它的准确含义是积分关系

$$
\boldsymbol Z_t=\boldsymbol Z_0
+\int_0^t\boldsymbol a(\boldsymbol Z_s,s)ds
+\int_0^t B(\boldsymbol Z_s,s)d\boldsymbol W_s.
$$

第一项积累确定性漂移，第二项积累随机增量。对 Itô 积分，每一小步的系数使用该步新噪声到来前的信息。状态依赖的系数仍然可以是随机的，但不会提前使用尚未发生的噪声；4.3 节已经解释这个先后关系如何影响微积分规则。

这样，Langevin 与 Itô 的联系就清楚了：Langevin 组织物理上的受力、耗散和涨落；Itô 方程为其中的白噪声扩散提供统一且明确的状态空间表述。有记忆或有色噪声的广义 Langevin 方程通常还需要扩充状态或保留历史，不能无条件地压成同一有限维 Itô 方程。

### 用 Itô 公式核对概率输运

物理收支给出了方程的结构，Itô 公式则将它与轨道模型精确对照。对一般一维 Itô 过程 $dX=b(X)dt+\sigma(X)dW$，令 $a(x)=\sigma^2(x)/2$。对足够光滑的函数 $f$，

$$
\frac d{dt}\langle f(X_t)\rangle
=\int p(x,t)[b(x)f'(x)+a(x)f''(x)]dx.
$$

同一个平均也等于 $\int f(x)p(x,t)dx$。在边界项消失或与边界条件相容时，对右侧分部积分一次和两次，得到

$$
\int f\,\partial_t p\,dx
=\int f[-\partial_x(bp)+\partial_x^2(ap)]dx.
$$

等式对足够丰富的 $f$ 成立，因而

$$
\partial_t p=-\partial_x(bp)+\partial_x^2(ap),
\qquad J=bp-\partial_x(ap).
$$

常数 $a=D$ 时恢复前面的漂移扩散流；状态依赖时，导数必须作用于整个 $ap$。改变噪声强度的空间分布，会改变各处随机增量的统计，不能只在常数 $D$ 的 Fick 形式里替换一个符号。若原方程采用 Stratonovich 约定，还需先换成相应的 Itô 漂移。

**上一篇：**[Path Integrals](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/path_integrals/) · **[返回专题首页](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/)**。

需要继续计算事件与响应时，可分别查阅[生成元与首达选读](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/markov_generators/)、[罕见事件选读](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/rare_events/)和[MSRJD 选读](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/path_integral_methods/)。
