一个粒子长期停留在左势阱，随后很快进入右势阱。单次转换的等待时间很随机，但如果收集很多次已经完成转换的历史，爬升阶段会不会呈现某种共同形状？在确定性力一直向下拉的地方，噪声需要怎样配合才能把粒子送上势垒？

路径概率把这个问题变成可计算的比较：给定相同的起点与终点，哪些历史需要较小的随机驱动代价？本篇先从离散增量建立路径权重，再讨论弱噪声下的稀有路径，最后用泛函同时组织时间关联与响应。

## 3.1 概率的对象从一个终点变成整条历史

考虑加性噪声方程

$$
dX_t=b(X_t)dt+\sqrt{2D}\,dW_t.
$$

把时间离散成 $N$ 步，步长 $h=T/N$。Itô 的 Euler 离散给出

$$
x_{n+1}-x_n=b(x_n)h+\sqrt{2Dh}\,Z_n,
\qquad Z_n\sim\mathcal N(0,1).
$$

固定初值 $x_0$ 后，各步的条件密度相乘，得到离散模型的联合路径密度

$$
\mathcal P_h(x_1,\ldots,x_N\mid x_0)
=\prod_{n=0}^{N-1}\frac{1}{\sqrt{4\pi Dh}}
\exp\left[-\frac{(x_{n+1}-x_n-b(x_n)h)^2}{4Dh}\right].
$$

每一项衡量实际增量偏离局部确定性增量多少。把所有中间位置积分掉，就得到终点转移密度：

$$
p(x_N,T\mid x_0)
=\lim_{h\to0}\int dx_1\cdots dx_{N-1}\,\mathcal P_h.
$$

有限步长下右边是 Euler 模型的传播核；在适当条件下取极限，才恢复连续扩散。路径积分记号是这个有明确测度和离散规则的极限的简写。

单条精确连续轨道通常具有零概率，因而不能把路径权重当作“恰好发生这条曲线”的普通概率。它可以表示离散观测值的密度、路径邻域的概率，或相对于某个参考路径测度的权重；使用哪一种含义，需要同时给出约定。

## 3.2 为什么有限噪声的作用量需要指定离散化

把离散平方和形式上写成积分，会出现

$$
\frac{1}{4D}\int_0^T[\dot x-b(x)]^2dt.
$$

但典型 Brownian 轨道没有普通导数，单靠这个表达式并没有定义完整路径测度。离散方式、Jacobian、归一化和端点约定都可能参与连续表示。

例如在加性噪声的中点离散中，设 $\bar x_n=(x_n+x_{n+1})/2$，将增量写成 $x_{n+1}-x_n=b(\bar x_n)h+\sqrt{2Dh}Z_n$。从 $Z_n$ 换元到 $x_{n+1}$，每步会产生

$$
\left|1-\frac h2b'(\bar x_n)\right|
\simeq\exp[-hb'(\bar x_n)/2]
$$

这一 Jacobian 因子。相应的常用 Onsager–Machlup 表示写成

$$
S_{\mathrm{OM}}[x]
=\int_0^T\left\{
\frac{[\dot x-b(x)]^2}{4D}+\frac12b'(x)
\right\}dt,
$$

并配合该中点规则的路径积分测度。在左端 Itô 离散中，噪声到增量的 Jacobian 对加性噪声为常数，相应的信息由不同的离散求值规则承载。不能把一种规则下的作用量和另一种规则下的测度任意拼接。

对于光滑参考曲线附近的小管概率，Onsager–Machlup 泛函也有相应解释，但管的定义和参考测度必须固定。这里需要记住的是：**在有限噪声下比较整条历史，平方偏差项之外的因素可能影响结果**。[Chow 与 Buice 的路径积分综述](https://link.springer.com/article/10.1186/s13408-015-0018-5)详细讨论了随机离散约定及响应场表示。

## 3.3 弱噪声极限：把历史理解为控制噪声的代价

如果目标是 $D\to0$ 时罕见事件的主导指数，描述会更简单。在适当的漂移正则性与非爆炸条件下，固定时间区间上的弱噪声路径大偏差泛函为

$$
\boxed{
I_T[x]=\frac14\int_0^T[\dot x-b(x)]^2dt.}
$$

这里 $x$ 是从指定初值出发的绝对连续候选路径；不满足相应条件的路径取无穷代价。围绕一类路径的事件 $A$，其概率的指数尺度由 $\inf_{x\in A}I_T[x]$ 控制，形式上记为

$$
\Pr(X\in A)\asymp
\exp\left[-\frac{\inf_{x\in A}I_T[x]}{D}\right].
$$

严格结论对闭集和开集分别给出上、下界；只有极小代价相容的事件才可直接写这种对数等价。对一条光滑曲线说“$X\approx x$”时，也是在指足够窄的路径邻域，并采用相应极限。

为了让粒子沿候选历史运动，需要提供额外速度

$$
u(t)=\dot x(t)-b(x(t)).
$$

把它看作受控动力学 $\dot x=b(x)+u$，则 $I_T=\int u^2dt/4$。高斯噪声对持续、大幅偏离的支撑代价更高。沿确定性方程 $\dot x=b(x)$ 的轨道有零作用量；需要噪声长期抵抗漂移的历史具有正作用量。

这里的“代价”首先是负对数概率的指数尺度。只有另行指定热浴与物理功的关系，才能给它更多热力学解释。

## 3.4 双稳势中的爬坡历史

令 $b=-U'$。对于从 $x_0$ 到 $x_T$ 的绝对连续曲线，恒等式

$$
(\dot x+U')^2=(\dot x-U')^2+4\dot xU'
$$

给出

$$
\boxed{
I_T[x]=\frac14\int_0^T[\dot x-U'(x)]^2dt
+U(x_T)-U(x_0).}
$$

若从阱底 $x_a$ 上升到相邻势垒顶 $x_s$，则

$$
I_T[x]\ge U(x_s)-U(x_a)=\Delta U.
$$

使非负积分消失的条件是 $\dot x=+U'(x)$：沿与确定性下坡相反的方向爬升。对 $U=(x^2-1)^2/4$，左阱到势垒之间有 $U'(x)=x(x^2-1)>0$，最优爬坡方向正好向右。

### 为什么要允许足够长的时间

精确从稳定点出发，反向方程也会永远停留在该固定点。因此，从阱底到势垒顶的最小代价 $\Delta U$ 一般是在允许时间趋长、端点邻域趋小的意义下逼近；它不意味着有限时间内存在一条从固定点启动并精确到达另一个固定点的反向确定性轨道。

实际转换轨迹会先经历长时间的阱内涨落，随后沿接近最优爬坡的通道进入势垒附近，再沿确定性下坡落入另一势阱。下坡段的弱噪声作用量为零。因而离开一个阱所需的指数代价由向势垒上升的那一段决定，等待时间与跨越本身的持续时间可以显著不同。

对弱噪声、光滑非退化的单势垒问题，这给出

$$
k_{\mathrm{escape}}\asymp e^{-\Delta U/D},\qquad
\mathbb E\tau_{\mathrm{escape}}\asymp e^{\Delta U/D}.
$$

这是速率和平均时间的主导指数，尚未给出前因子。若已达到准稳态逃逸状态，且 $kT\ll1$，固定观察窗内的事件概率才近似为 $kT$；不能把 $e^{-\Delta U/D}$ 直接当成任意观察时间的完整逃逸概率。[Berglund 的综述](https://arxiv.org/abs/1106.5799)将这种大偏差观点与首达时间、Kramers 公式联系起来。

![双稳势中的确定性下坡与最优爬坡](../figures/04_optimal_paths.svg)

*图 4｜左图在势函数上标出确定性下坡与弱噪声最优爬坡的方向。右图比较从 $-0.95$ 到 $-0.05$、用时相同的三条光滑候选历史：反向漂移轨道以及它的两种扰动。反向轨道使作用量达到这两个端点之间的势差；端点避开固定点，使有限时间比较有明确含义。图中候选曲线不是典型 Brownian 样本。*

### 最优路径怎样推广

一般漂移下，极小化作用量可引入共轭变量

$$
q=\frac12(\dot x-b),\qquad
H(x,q)=qb(x)+q^2.
$$

Euler–Lagrange 方程等价于 Hamilton 方程

$$
\dot x=b(x)+2q,\qquad \dot q=-q\,b'(x).
$$

它们将“寻找最省噪声的历史”化为带端点条件的动力学问题。多维、常数各向同性扩散时也有 $I=\int|\dot{\boldsymbol x}-\boldsymbol b|^2dt/4$。如果漂移不再是梯度，最优上升路径一般不能通过把原来的确定性运动简单反向得到；这个结论依赖前面使用的可逆梯度结构。

## 3.5 OU 桥：路径极值怎样接回一个精确概率

稀有路径不只用于跨势垒。给 OU 指定 $X_0=x_0$、$X_T=x_f$，即使 $x_f$ 远离自由均值，也可以问哪条光滑历史最省代价。

令控制满足 $\dot x=-\gamma x+u$，则终点偏离为

$$
\delta=x_f-x_0e^{-\gamma T}
=\int_0^T e^{-\gamma(T-s)}u(s)ds.
$$

由 Cauchy–Schwarz 不等式，

$$
\delta^2\le
\left[\int_0^T e^{-2\gamma(T-s)}ds\right]
\left[\int_0^T u^2(s)ds\right].
$$

因此

$$
\boxed{I_T^{\min}
=\frac{\gamma(x_f-x_0e^{-\gamma T})^2}{2(1-e^{-2\gamma T})}.}
$$

等号对应 $u(s)$ 正比于 $e^{-\gamma(T-s)}$，积分得到

$$
x_*(t)=x_0e^{-\gamma t}
+(x_f-x_0e^{-\gamma T})\frac{\sinh(\gamma t)}{\sinh(\gamma T)}.
$$

第一篇的精确转移密度，其指数恰好是 $-I_T^{\min}/D$。高斯前因子则来自围绕极值历史的涨落积分。线性高斯模型中，这条光滑曲线也等于给定两端点后的条件均值轨道；实际桥轨道仍然粗糙，在它周围涨落。

这个结果展示了 saddle-point 方法的意义：先找出满足约束的主导历史，再积分它附近的涨落。对一般非线性模型，弱噪声允许系统地组织这种近似；对 OU，二次结构使计算可以精确完成。

## 3.6 用泛函同时组织关联与响应

路径方法还可以研究典型涨落。加入外源 $h(t)$，并用另一个源 $j(t)$ 记录位置统计：

$$
Z[j,h]=\mathbb E_h\exp\left[\int j(t)X_tdt\right],
\qquad W[j,h]=\log Z[j,h].
$$

对 $j$ 的一次导数给出均值，二次导数给出连通关联：

$$
\left.\frac{\delta W}{\delta j(t)}\right|_{j=0}=\langle X_t\rangle_h,
\qquad
\left.\frac{\delta^2W}{\delta j(t)\delta j(s)}\right|_{j=h=0}
=\langle X_tX_s\rangle_c.
$$

再对 $h(s)$ 求导，则得到对物理扰动的响应。$j$ 是提取统计的记账工具，$h$ 实际改变漂移，二者有不同作用。

### 从噪声积分引入响应场

形式上把方程约束 $\dot x-b(x)-h=\eta$ 用 Fourier 场 $\hat x$ 表示。噪声满足 $\langle\eta(t)\eta(s)\rangle=2D\delta(t-s)$，积分掉它以后，Itô 离散的加性噪声路径表示为

$$
Z[j,h]=\int\mathcal Dx\,\mathcal D\hat x\;p_0(x_0)
\exp\left\{\int dt\left[
i\hat x(\dot x-b-h)-D\hat x^2+jx
\right]\right\}.
$$

这里的测度仍按左端离散定义，固定初始准备，终点自由积分，归一化使 $Z[0,h]=1$。Jacobian 的处理也与这一规则绑定。

令 $\widetilde x=-i\hat x$，常用的 MSRJD 作用量记为

$$
S[x,\widetilde x]=\int dt\left[
\widetilde x(\dot x-b(x))-D\widetilde x^2
\right],
$$

从而指数为 $-S+\int(jx+h\widetilde x)dt$。响应场沿相应的虚轴轮廓积分；含 $-D\widetilde x^2$ 的形式不是一个可以沿实轴直接当作概率权重积分的正作用量。

响应于是写成

$$
R(t,s)=\left.\frac{\delta^2W}{\delta j(t)\delta h(s)}\right|_0
=\langle x(t)\widetilde x(s)\rangle_c.
$$

它描述时刻 $s$ 的驱动怎样影响未来位置，因果性要求 $t\lt s$ 时为零。左端 Itô 规则还规定了等时响应的处理，不能脱离离散规则任意选择 $\Theta(0)$。

### OU 的二次泛函

取 $b=-\gamma x$，作用量为

$$
S_{\mathrm{OU}}=\int dt\left[
\widetilde x(\partial_t+\gamma)x-D\widetilde x^2
\right].
$$

其迟致 Green 函数满足 $(\partial_t+\gamma)R(t,s)=\delta(t-s)$，给出

$$
R(t,s)=\Theta(t-s)e^{-\gamma(t-s)}.
$$

在稳态准备下，噪声经两个响应核传播，得到

$$
C(t,s)=2D\int_{-\infty}^{\min(t,s)}R(t,u)R(s,u)du
=\frac D\gamma e^{-\gamma|t-s|}.
$$

这与第一篇的直接求解完全一致。泛函表示的价值是把更多时刻、更高阶关联和非线性修正放在同一套源导数中组织；对非线性漂移，可以围绕可解的高斯部分展开，而无需为每个观测量重新设计一次计算。

随机历史的视角同时覆盖了两种问题：典型涨落由路径集合的关联和响应描述，稀有事件由约束下的最小噪声代价组织。下一篇将回到局部演化规则，直接从生成元计算事件时间，并与这里的势垒指数对照。

**上一篇：**[Fokker–Planck Dynamics](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/fokker_planck/) · **下一篇：**[Markov Processes and Generators](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/markov_generators/)。
