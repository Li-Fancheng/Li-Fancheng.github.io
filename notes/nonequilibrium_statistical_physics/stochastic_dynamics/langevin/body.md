在统计物理中，我们通常不需要追踪系统的全部微观自由度，在非平衡统计物理中我们将延续这个想法。我们真正关心的往往只是其中一小部分，例如粒子的位置和速度、局部密度变化，或某些集体行为等。但这些变量并不是孤立演化的，它们始终受到其余大量自由度的影响。即使完整系统遵循确定性的微观动力学，仅凭我们所关心变量的当前状态，通常也不足以确定它们未来的运动，因为那些未被显式描述的自由度仍然参与着系统的演化。如何在不追踪全部微观细节的情况下，建立少数变量的有效动力学，是统计物理中的一个基本问题。

一种自然的思路是将未解析自由度的影响进行统计处理。环境的作用往往包含相对稳定、可以重复的平均部分，也包含因微观状态不同而产生的涨落。例如，流体会对运动的粒子产生阻力，同时又不断通过分子碰撞向它传递动量。但仅仅区分平均作用与涨落还不够：这些涨落有多强？不同时刻的扰动是否相关？环境对粒子过去运动的影响能够持续多久？在什么条件下，复杂的环境作用可以近似为简单的阻尼和随机力？这些都是建立有效描述时需要回答的物理问题。

随机动力学为研究这类问题提供了一种方法。我们不再试图仅凭少数变量预测唯一的运动轨迹，而是描述它们可能经历的不同运动，以及这些运动所服从的统计规律。虽然单次运动可能难以预测，但大量轨迹的统计性质，例如关联函数、弛豫时间和扩散系数，却可以具有明确而可重复的规律。

本篇将从最经典的 Brownian 运动出发，以流体中的微粒为例，考察如何从牛顿方程得到 Langevin 方程。我们将看到，环境引起的阻尼与涨落并不是彼此无关的两种作用；它们共同决定粒子的热运动，并使粒子从短时间的弹道运动逐渐过渡到长时间的扩散。

<span id="12-从牛顿方程到开放系统"></span>
<span id="13-消去环境以后记忆阻尼与随机力"></span>
<span id="环境为什么会记住粒子的过去"></span>

## 1.1 从牛顿方程到 Langevin 方程 {#bath-memory}
考虑一个与环境相互作用的粒子，将环境自由度统记为 $\boldsymbol q$。粒子的牛顿方程为

$$
m\ddot x=-\frac{dV}{dx}
+F_{\mathrm{bath}}(x,\boldsymbol q,\dot{\boldsymbol q}).
$$

环境本身也满足牛顿方程，因此

$$
m_i\ddot q_i=F_i(x,\boldsymbol q),\qquad i=1,\ldots,N.
$$

给定全部初始位置和速度，这套方程的演化仍然是确定的。所谓消去环境自由度，是先在给定粒子运动历史和环境初态的条件下求解 $\boldsymbol q(t)$，再代回粒子方程。环境曾被粒子扰动，之后又会反过来影响粒子，因此代回后的力通常依赖此前的运动历史。对未指定的环境初态作统计平均，才得到平均作用及其涨落；随机性来自我们没有逐一保留的环境状态。

在接近平衡、空间近似均匀的热浴中，若粒子运动足够缓慢，可以将平均阻力对速度历史的响应保留到一阶。将平均保守作用吸收到 $V(x)$，并采用与平衡热浴相容的初始准备，可写出一种广义 Langevin 形式：

\[
m\dot v(t)=-\left.\frac{dV}{dx}\right|_{x=x(t)}
-\int_0^t K(t-s)v(s)\,ds+\eta(t),
\qquad\dot x=v.
\]

$K(u)$ 是因果的摩擦记忆核，描述时刻 $s$ 的运动在延迟 $u=t-s$ 后仍引起多强的环境反作用；$\eta(t)$ 是均值为零、一般具有时间关联的剩余涨落力。这个卷积已经采用了线性响应和热浴平稳性的简化，更一般的约化还可能包含依赖位置的记忆和噪声。广义 Langevin 描述及其平衡约束可参阅 [Kubo（1966），第 4 节](https://www.zib.de/userpage/donati/stochastics2023/06/references/Kubo1966.pdf)。

现在假设 $K(u)$ 主要集中在很短的环境弛豫时间 $\tau_{\mathrm{bath}}$ 内，而粒子速度显著变化的时间 $\tau_v$ 远长于它。在初始环境瞬态过去后，记忆区间内可用 $v(t-u)\simeq v(t)$，于是

\[
\begin{aligned}
\int_0^t K(t-s)v(s)\,ds
&=\int_0^t K(u)v(t-u)\,du\\
&\simeq v(t)\int_0^\infty K(u)\,du
=\zeta v(t),
\end{aligned}
\]

其中 $\zeta=\int_0^\infty K(u)du$。例如，$K(u)=(\zeta/\tau_{\mathrm{bath}})e^{-u/\tau_{\mathrm{bath}}}$ 的宽度是 $\tau_{\mathrm{bath}}$、面积是 $\zeta$；短记忆近似保留这个面积，忽略反作用在时间上的展宽。这一步才把非局域摩擦化为瞬时摩擦。时间尺度分离对这一步误差的作用，见 [Di Pasquale 等关于 Markov 近似的分析](https://bura.brunel.ac.uk/handle/2438/24744)。

对涨落还需要进一步的统计近似：选取 $\tau_{\mathrm{bath}}\ll\Delta t\ll\tau_v$，使一段时间内包含许多微小碰撞，而粒子的慢变量变化很小。在环境充分混合、单次冲量方差有限且没有少数碰撞占主导等条件下，累积冲量可近似为 Gaussian，不同粗粒化时间段的冲量近似独立，从而得到 Gaussian 白噪声描述。短关联本身并不保证 Gaussian 性；两者分别依赖时间尺度分离和许多弱扰动的统计极限。

最后，还需假设在所研究的空间范围内，摩擦系数与涨落强度对位置的依赖可以忽略，才能使用常数摩擦和加性噪声。于是得到本文采用的模型

\[
m\ddot x=-\frac{dV}{dx}-\zeta\dot x+\eta(t).
\]

这套近似保留了 $(x,v)$ 及环境的有效阻尼和随机冲量，舍去了环境的详细状态与快速记忆，使 $(x,v)$ 构成封闭的 Markov 描述。如果环境记忆不能忽略，或噪声强度明显依赖位置，就需要回到带记忆或乘性噪声的模型。下一节再利用热平衡确定白噪声的强度。

<span id="14-白噪声极限与-langevin-方程"></span>
<span id="从有限时间的随机冲量理解噪声"></span>
<span id="得到带随机力的牛顿方程"></span>

## 1.2 热平衡情况下的涨落-耗散关系 {#thermal-noise}

采用白噪声近似，研究阻尼与随机力如何共同决定粒子的热运动可以得到

$$
m\ddot x=-\frac{dV}{dx}-\zeta\dot x+\eta(t),
\qquad \langle\eta(t)\rangle=0,
\qquad \langle\eta(t)\eta(s)\rangle=A\delta(t-s).
$$

这里的 \(A\) 衡量随机力的强度。虽然随机力的平均值为零，但在短时间 \(h\) 内，它仍会给粒子带来随机冲量，其均方为 \(\langle[\int_t^{t+h}\eta(s)\,ds]^2\rangle=Ah\)，类似随机游走。因此，随机冲量的典型大小随 \(\sqrt{h}\) 增长，而恒定外力产生的冲量随 \(h\) 增长。

为了简化问题，先不考虑外部势能，并记速度弛豫时间 $\tau_m=m/\zeta$。引入速度得到一个一阶线性微分方程，使用积分因子 $e^{t/\tau_m}$ 求解得到（逐步推导见[附录 A.1](#integrating-factor)）

$$
v(t)=v_0e^{-t/\tau_m}
+\frac1m\int_0^t e^{-(t-s)/\tau_m}\eta(s)ds.
$$

这个解将速度分为初值的衰减和后续随机冲量的累积。指数因子 $e^{-(t-s)/\tau_m}$ 记录粒子自身的速度弛豫，前一节的 $K$ 则记录环境反作用的延迟。即使环境的快速记忆已被忽略，粒子速度仍可在有限时间内保持关联，这与 $(x,v)$ 的 Markov 性相容。如果初速度与后续随机力独立，两部分的方差便可以相加。将随机积分的平方写成重积分，再利用白噪声关联函数（见[附录 A.2](#squared-integral)），得到

$$
\begin{aligned}
\operatorname{Var}(v(t))
&=\operatorname{Var}(v_0)e^{-2t/\tau_m}
+\frac A{m^2}\int_0^t e^{-2(t-s)/\tau_m}ds\\
&=\operatorname{Var}(v_0)e^{-2t/\tau_m}
+\frac{A}{2\zeta m}(1-e^{-2t/\tau_m}).
\end{aligned}
$$

前面的计算表明，无论粒子最初具有怎样的速度分布，初始条件的影响都会随着阻尼逐渐消失。与此同时，环境不断通过随机碰撞向粒子传递动量，使速度始终保持一定程度的涨落。长时间后，粒子的平均速度趋于零，而速度方差趋于一个稳定值 \(A/(2\zeta m)\)。这个结果说明，阻尼与随机力可以共同维持稳定的速度涨落，但我们还没有确定这个稳态方差究竟应该有多大，因为随机力的强度 \(A\) 仍然是未知的。

为了确定 \(A\)，进一步假设周围的液体是一个温度为 \(T\) 的平衡热浴，这个假设的适用范围非常广泛，粒子浸泡在温度固定的液体里。粒子不断与环境交换能量，经过足够长的时间后与环境达到热平衡。此时，粒子的速度涨落不能是任意大小，而必须符合我们在平衡态统计物理里知道的要求。根据能量均分定理，粒子一个自由度对应的的平均动能应当为 \(k_{\mathrm B}T/2\)，即

\[
\frac{1}{2}m\langle v^2\rangle_{\mathrm{st}}
=\frac{1}{2}k_{\mathrm B}T.
\]

由于稳态下的平均速度为零，这里的 \(\langle v^2\rangle_{\mathrm{st}}\) 正是我们刚才计算的长期速度方差。因此，将两个结果比较，就可以确定平衡热浴所产生的随机力究竟应该有多强。将这个强度代回受力式，便得到惯性 Langevin 方程：

$$
\boxed{m\ddot x=-\frac{dV}{dx}-\zeta\dot x
+\sqrt{2\zeta k_{\mathrm B}T}\,\xi(t)}.
$$

$\xi$ 表示单位强度 Gaussian 白噪声，$\langle\xi(t)\xi(s)\rangle=\delta(t-s)$。这样得到的 Langevin 方程仍然保留了熟悉的牛顿第二定律的形式，只是粒子除了受到保守力，还受到环境产生的摩擦力和随机力。对于平衡情况，这两种作用的强度不是相互独立的，摩擦不断耗散粒子的动能，随机碰撞又不断激发粒子的热运动，为了达到最终温度的热平衡态，二者必须满足刚才得到的关系\(A=2\zeta k_{\mathrm B}T\)。这是涨落-耗散关系（fluctuation–dissipation relation）的一个特殊情况。如果环境本身处于非平衡态，例如存在主动运动或外部随机驱动，我们仍然可以写出类似的 Langevin 方程，但随机力的强度一般不再由这个平衡关系决定。

## 1.3 速度弛豫与位置扩散 {#inertial-msd}

同一个粒子的运动，在不同时间尺度上会表现出很不一样的性质。在很短的时间内，阻尼和随机碰撞还来不及明显改变粒子的速度，因此粒子近似以当前速度做匀速运动；但随着时间增长，粒子原来的运动方向会逐渐被遗忘，最终表现为扩散运动。决定这一变化的特征时间，就是前面得到的速度弛豫时间 \(\tau_m=m/\zeta\)。

为了更清楚地描述粒子如何逐渐失去速度记忆，我们来计算速度的自关联函数。假设粒子的初速度已经服从温度为 \(T\) 的平衡分布，即 \(\langle v_0\rangle=0\)、\(\langle v_0^2\rangle=k_{\mathrm B}T/m\)，并且与之后的随机力相互独立。将前面求得的速度解 \(v(t)\) 与初速度 \(v_0\) 相乘，再对不同的随机轨迹取平均。由于随机力的平均值为零，且与初速度独立，最终得到

\[
C_v(t)=\langle v(t)v_0\rangle_{\mathrm{st}}
=\frac{k_{\mathrm B}T}{m}e^{-|t|/\tau_m}.
\]

这是一个指数衰减的速度自关联函数。一个容易误解的地方是：平衡态下粒子的平均速度为零，但这并不意味着不同时刻的速度之间没有关联。例如，如果粒子此刻正在向右运动，那么在很短的时间后，它仍然更有可能向右运动。随着时间增长，这种关联逐渐减弱，并在 \(t\gg\tau_m\) 时趋于零。因此，\(\tau_m\) 也可以理解为粒子保持速度记忆的特征时间。

除了时间尺度，我们还可以估计一个相应的空间尺度。平衡态下，粒子的典型热速度为 \(v_{\mathrm{th}}=\sqrt{k_{\mathrm B}T/m}\)。将它乘以速度弛豫时间，可以定义一个特征长度

\[
\ell_p=v_{\mathrm{th}}\tau_m
=\frac{\sqrt{mk_{\mathrm B}T}}{\zeta}.
\]

这个长度称为速度持续长度（persistence length），描述粒子在逐渐失去原来速度记忆的过程中，典型能够运动多远。因此，\(\tau_m\) 和 \(\ell_p\) 分别给出了从短时间弹道运动向长时间扩散运动过渡的时间尺度和空间尺度。

接下来，我们计算粒子的均方位移。由于位移是速度在时间上的积分，\(\Delta x(t)=\int_0^t v(s)\,ds\)，将其平方并取平均后，自然会出现不同时刻的速度乘积 \(\langle v(t_1)v(t_2)\rangle\)，也就是速度关联函数。因此，粒子位移的增长不仅取决于速度的大小，还取决于不同时刻的速度能够保持多强的关联。平方展开、积分区域对称性及换序的完整步骤见[附录 A.2](#squared-integral)。

前面我们已经假设初速度服从热平衡分布，因此速度过程从一开始就是平稳的（stationary）。平稳性意味着速度的统计性质不随时间改变，两个时刻的速度关联只取决于它们之间的时间差，而与具体的观测时刻无关，即 \(\langle v(t_1)v(t_2)\rangle=C_v(|t_1-t_2|)\)。利用这一性质，我们可以将位移的积分表达式代入均方位移，得到两个时间积分。由于交换 \(t_1\) 和 \(t_2\) 的先后顺序不会改变速度乘积的平均值，因此只需要计算 \(t_1\geq t_2\) 的部分。于是


\[
\begin{aligned}
\langle[\Delta x(t)]^2\rangle
&=\int_0^t dt_1\int_0^t dt_2\,
\langle v(t_1)v(t_2)\rangle\\
&=2\int_0^t dt_1\int_0^{t_1}dt_2\,C_v(t_1-t_2)\\
&=2\int_0^t d\tau\,(t-\tau)C_v(\tau).
\end{aligned}
\]

代入 $C_v(s)=(k_{\mathrm B}T/m)e^{-s/\tau_m}$ 后，只需计算 $\int_0^t(t-s)e^{-s/\tau_m}ds$。分部积分，取 $u=t-s$、$dv=e^{-s/\tau_m}ds$，得到

\[
\begin{aligned}
\int_0^t(t-s)e^{-s/\tau_m}\,ds
&=\left[-\tau_m(t-s)e^{-s/\tau_m}\right]_0^t
-\tau_m\int_0^t e^{-s/\tau_m}\,ds\\
&=\tau_m t-\tau_m^2(1-e^{-t/\tau_m}).
\end{aligned}
\]

乘回 $2k_{\mathrm B}T/m$，并用 $D=k_{\mathrm B}T/\zeta=(k_{\mathrm B}T/m)\tau_m$，得到

\[
\boxed{\langle[x(t)-x_0]^2\rangle
=2Dt-2D\tau_m(1-e^{-t/\tau_m}).}
\]

边界项 $\tau_m t$ 给出线性增长的 $2Dt$；余下的积分给出有限速度弛豫带来的暂态修正 $-2D\tau_m(1-e^{-t/\tau_m})$。这两项共同描述弹道到扩散的过渡。

在白噪声描述适用的短时窗口 $\tau_{\mathrm{bath}}\ll t\ll\tau_m$ 内，展开指数后，两项中的线性部分恰好相消，留下 $\langle\Delta x^2(t)\rangle\simeq(k_{\mathrm B}T/m)t^2$。此时速度尚未明显改变，粒子近似匀速运动，因此出现弹道增长。

当 $t\gg\tau_m$ 时，修正项趋于常数 $-2D\tau_m$，长期增长率为 $2D$，所以均方位移的主导项为 $\langle\Delta x^2(t)\rangle\sim2Dt$，对应正常扩散。

![惯性 Brownian 粒子的速度相关与弹道到扩散跨越](../figures/09_inertial_crossover.svg)

*Fig. 1. 左图展示均方位移从短时 \(t^2\) 增长到长时 \(t\) 增长的过渡；右图展示速度关联函数在时间尺度 \(\tau_m\) 上的指数衰减。时间以 \(\tau_m\) 为单位，均方位移以 \(D\tau_m\) 为单位。*

除了得到完整的均方位移，我们还可以从速度关联函数直接确定长时间的扩散系数。对均方位移关于时间求导，有 \(\frac{d}{dt}\langle\Delta x^2(t)\rangle=2\int_0^t C_v(s)\,ds\)。当观察时间远大于速度关联时间时，积分已经接近一个常数，而此时均方位移的增长率正是 \(2D\)。因此，

\[
D=\int_0^\infty C_v(s)\,ds
=\frac{k_{\mathrm B}T}{m}\tau_m
=\frac{k_{\mathrm B}T}{\zeta}.
\]

这就是扩散系数的 Green–Kubo 关系。它给出了一个很直观的物理图像：粒子的速度涨落越强，或者速度能够在同一方向上保持得越久，位置扩散就越快。在当前模型中，质量较大的粒子具有较小的热速度，但其速度弛豫时间 \(\tau_m=m/\zeta\) 也更长。这两种质量依赖恰好相互抵消，因此长时间的扩散系数 \(D\) 与粒子质量无关。不过，质量仍然决定了从弹道运动过渡到扩散运动所需的时间。[附录 C](#green-kubo-extension)进一步将这个积分与迁移率、多维扩散及其他输运系数联系起来。

最后，我们还可以将扩散与粒子受到外力时的响应联系起来。假设对粒子施加一个恒定的外力 \(f\)，那么平均速度满足 \(m\,d\langle v\rangle/dt=-\zeta\langle v\rangle+f\)。经过足够长的时间，粒子达到稳定的漂移速度 \(\langle v\rangle_{\mathrm{st}}=f/\zeta\)。因此，定义迁移率（mobility）为单位外力引起的稳态速度响应，就有 \(\mu=\partial\langle v\rangle_{\mathrm{st}}/\partial f=1/\zeta\)。

将这一结果与前面得到的扩散系数比较，便得到

\[
D=\mu k_{\mathrm B}T.
\]

这就是著名的 Einstein 关系。它将粒子在没有外力时的自发扩散，与受到外力时的运动响应联系了起来。两者看起来是不同的现象，但实际上都受到同一个环境的影响：液体一方面通过摩擦阻碍粒子运动，另一方面又通过热涨落不断推动粒子。对于平衡热浴，涨落与耗散之间的联系最终保证了 Einstein 关系成立。

<span id="15-从-langevin-到-itô在状态空间中统一描述"></span>
<span id="二阶运动方程怎样变成一阶状态方程"></span>
<span id="再消去速度过阻尼动力学"></span>
<span id="brownian-位置为什么不否定原来的牛顿轨道"></span>

## 1.4 从惯性 Langevin 到过阻尼动力学 {#overdamped}

要进一步只保留位置，需要让速度的弛豫快于位置的慢演化。回到惯性方程

\[
m\ddot x=F(x)-\zeta\dot x+\sqrt{2\zeta k_{\mathrm B}T}\,\xi(t),
\]

其中 $F(x)=-\frac{dV}{dx}$，并继续假设 $\zeta,T$ 为常数。记 $\tau_x$ 为位置在所关心的空间尺度上明显演化、或外力明显变化所需的较短时间。过阻尼约化要求存在中间观测窗口

\[
\tau_{\mathrm{bath}}\ll\tau_m=\frac m\zeta
\ll\Delta t\ll\tau_x.
\]

第一层分离使环境可近似为白噪声热浴；第二层分离使速度能够在位置显著改变前弛豫。例如，谐势中的慢回复时间为 $\tau_x=\zeta/\kappa$，要求 $m\kappa/\zeta^2\ll1$。对特征速度为 $U$ 的缓变平均运动，惯性力与摩擦力的量级比为

\[
\frac{mU/\tau_x}{\zeta U}=\frac{\tau_m}{\tau_x}\ll1.
\]

随机速度含有快速涨落，还需比较有限时间内的冲量。将惯性方程积分，得到精确收支

\[
\zeta\Delta x=\int_t^{t+\Delta t}F(x(s))\,ds
+\int_t^{t+\Delta t}\eta(s)\,ds-m\Delta v.
\]

对速度已弛豫后的热涨落，热速度的量级为 $v_{\mathrm{th}}=\sqrt{k_{\mathrm B}T/m}$。当 $\Delta t\gg\tau_m$ 时，端点动量差的均方根量级为 $\sqrt{mk_{\mathrm B}T}$，而随机冲量的均方根为 $\sqrt{2\zeta k_{\mathrm B}T\Delta t}$，因此

\[
\frac{|m\Delta v|_{\mathrm{rms}}}
{\left|\int_t^{t+\Delta t}\eta(s)\,ds\right|_{\mathrm{rms}}}
\sim\sqrt{\frac{\tau_m}{\Delta t}}\ll1.
\]

这就说明，在粗粒化位移中可以忽略端点惯性修正。同时 $\Delta t\ll\tau_x$ 允许用 $F(x(t))\Delta t$ 近似外力积分。除以 $\zeta$，并利用噪声冲量的方差 $2\zeta k_{\mathrm B}T\Delta t$，得到

\[
\Delta x\simeq\mu F(x)\Delta t+\sqrt{2D\Delta t}\,Z,
\qquad Z\sim\mathcal N(0,1).
\]

其中 $\mu=1/\zeta$、$D=k_{\mathrm B}T/\zeta$。近似保留了平均漂移 $\mu F\Delta t$ 和位移方差 $2D\Delta t$，消去了速度的快弛豫及其在相邻时间段之间留下的短时关联。引入标准 Wiener 过程 $W_t$，其不相交区间的增量独立，且 $\Delta W\sim\mathcal N(0,\Delta t)$，便可写成

\[
dx=\mu F(x)\,dt+\sqrt{2D}\,dW_t.
\]

这就是位置本身构成 Markov 状态的过阻尼 Langevin 方程。这里的 $dt$ 是慢时间模型中的连续时间记号；推导要求存在上述中间窗口，不能把它用于分辨固定质量粒子在 $t\lesssim\tau_m$ 时的弹道段。有效 Brownian 位置轨迹的不可微性，也来自对速度快过程的消去。关于这一极限的物理说明，可参阅 [Cugliandolo 讲义，第 2.3.7 节](https://giamarchi.unige.ch/wp-content/php_code/notes/cuso/Cugliandolo/DynamicsComplexSystems-lectures1-4.pdf)。

<span id="17-ou-过程涨落如何被产生又如何被消除"></span>
<span id="精确解与有限时间传播"></span>
<span id="矩方程中的机制"></span>

## 1.5 从自由扩散到束缚涨落：Ornstein–Uhlenbeck 过程 {#ou-trajectories}

上一节我们得到了粒子位置的过阻尼 Langevin 方程。对于自由粒子，由于没有外力限制其运动范围，均方位移会随时间不断增长。但如果粒子受到一个指向稳定平衡位置的回复力，情况就会不同：随机碰撞仍然使粒子不断运动，而回复力则倾向于将它拉回平衡位置。两种作用共同决定了粒子能够偏离平衡位置多远，以及这种偏离能够持续多长时间。

我们先考虑最简单的情况：粒子处于一个谐势中，

\[
V(x)=\frac{1}{2}\kappa x^2,
\]

其中 $\kappa>0$ 是势能的刚度，平衡位置取在 $x=0$。相应的回复力为 $F(x)=-\kappa x$，即粒子偏离平衡位置越远，受到的回复力越大，方向始终指向原点。将这一回复力代入上一节得到的过阻尼 Langevin 方程，便有

\[
dx=-\gamma x\,dt+\sqrt{2D}\,dW_t,
\qquad
\gamma=\frac{\kappa}{\zeta}.
\]

这就是著名的 **Ornstein–Uhlenbeck（OU）过程**。其中，$\gamma$ 描述位置偏离平衡点后恢复的快慢，$D$ 则决定随机扰动的强度。我们仍然采用上一节的过阻尼近似，即速度的弛豫远快于位置的回复，$\tau_m\ll\gamma^{-1}$。

虽然我们在这个体系中得到了 OU 过程，但其重要性并不局限于谐势。对于一个更一般的随机动力学系统，只要它存在稳定固定点 (stable fixed point) $x_*$，就可以在固定点附近将确定性漂移展开为

\[
b(x)\simeq b'(x_*)(x-x_*),
\]

其中 $b(x_*)=0$，稳定性要求 $b'(x_*)<0$。如果噪声强度在这一小范围内也可以近似看作常数，那么将坐标原点移到 $x_*$，便重新得到 OU 方程。OU 过程是描述稳定状态附近小涨落的一个基本模型，线性恢复力模型的广泛性也正体现于此。

为了理解 OU 过程的运动，我们先暂时忽略随机扰动。此时方程退化为普通的确定性回复运动，$\dot x=-\gamma x$，其解为 $x(t)=x_0e^{-\gamma t}$。这表明，如果粒子最初偏离平衡位置，它会逐渐返回原点，偏离量以指数形式衰减。相应的回复时间为 $\tau=1/\gamma$。

加入随机扰动后，情况变得更加有趣。粒子在返回平衡位置的过程中，仍然不断受到环境的随机推动。每一次扰动都会使位置发生微小变化，而这些变化在后续运动中又受到回复力的作用，逐渐衰减。因此，粒子在某一时刻的位置，既包含初始位置尚未衰减的部分，也包含此后所有随机扰动留下的影响。由于 OU 方程是线性的，用[附录 A.1](#integrating-factor)中的积分因子法可以直接求出它的精确解：

\[
x(t)=x_0e^{-\gamma t}
+\sqrt{2D}\int_0^t e^{-\gamma(t-s)}\,dW_s.
\]

这个解很好地体现了回复力与随机扰动的共同作用。第一项表示初始位置的影响，它会随着时间逐渐消失。第二项则是从初始时刻到现在的所有随机扰动的累积。

特别值得注意的是积分中的指数因子 $e^{-\gamma(t-s)}$。它表明，在较早时刻 $s$ 产生的一次扰动，到时刻 $t$ 时只保留了原来的一部分。时间间隔越长，留下的影响就越小。也就是说，OU 过程会逐渐忘记过去的扰动，记忆衰减的特征时间正是 $\gamma^{-1}$。这与自由 Brownian 运动有一个重要区别。自由粒子的一次随机位移会一直保留在后续位置中，而在 OU 过程中，同样的位移会被回复力逐渐消除。尽管环境施加的仍然是没有时间关联的白噪声，粒子的位置却会因为自身的回复动力学而表现出有限的时间关联。

接下来考虑从同一个初始位置 $x_0$ 出发的许多条独立轨迹。由于随机增量的平均值为零，粒子的平均位置只保留初始条件的贡献。另一方面，不同时刻的噪声增量相互独立，因此可以将它们对位置方差的贡献相加。利用上面的精确解，得到

\[
\begin{aligned}
\langle x(t)\rangle
&=x_0e^{-\gamma t},\\[6pt]
\left\langle [x(t)-\langle x(t)\rangle]^2\right\rangle
&=2D\int_0^t e^{-2\gamma(t-s)}\,ds\\
&=\frac{D}{\gamma}\left(1-e^{-2\gamma t}\right).
\end{aligned}
\]

这里可以看到两个不同的过程。平均位置随着时间逐渐接近平衡点，而位置方差则从零开始增长，最终趋于一个有限值 $D/\gamma$。

在短时间 $t\ll\gamma^{-1}$ 内，回复力还来不及明显改变粒子的运动。此时将指数函数展开，得到位置方差约为 $2Dt$，与自由 Brownian 扩散相同。这里的短时间仍然是在过阻尼描述适用的尺度上，即 $\tau_m\ll t\ll\gamma^{-1}$。

而在长时间 $t\gg\gamma^{-1}$ 后，粒子已经基本失去初始位置的信息，平均位置趋于零，方差则达到稳定值

\[
\sigma_{\mathrm{st}}^2=\frac{D}{\gamma}.
\]

方差达到稳定并不意味着粒子停止运动。它仍然不断受到随机碰撞，有时靠近中心，有时远离中心。只是回复力阻止了这些随机位移无限积累，使粒子的位置始终分布在平衡点附近的一个有限范围内。因此，这个问题上达到稳态的是位置的概率分布，而不是某一条轨迹。

对于前面考虑的平衡热浴，$D$ 和 $\gamma$ 并不是两个完全独立的参数。利用 Einstein 关系 $D=k_{\mathrm B}T/\zeta$，以及 $\gamma=\kappa/\zeta$，可以得到

\[
\boxed{
\sigma_{\mathrm{st}}^2=\frac{k_{\mathrm B}T}{\kappa},
\qquad
\tau=\frac{\zeta}{\kappa}.
}
\]

这两个结果具有非常清楚的物理意义。

首先，平衡位置涨落的大小由温度和势能刚度决定。温度越高，热运动越强，粒子可以偏离平衡位置更远；势能越陡，将粒子推离平衡位置所需的能量越大，因此位置涨落越小。事实上，$\sigma_{\mathrm{st}}^2=k_{\mathrm B}T/\kappa$ 正是平衡态能量均分定理的结果（振动自由度）：

\[
\frac{1}{2}\kappa\langle x^2\rangle_{\mathrm{st}}
=\frac{1}{2}k_{\mathrm B}T.
\]

实际上，我们还可以进一步确定粒子的完整稳态分布。前面求得的 OU 精确解表明，位置由确定性衰减项和 Gaussian 随机积分组成，因此从固定初态出发，任意时刻的位置分布都是 Gaussian 的。随着初始位置的影响逐渐消失，长期分布的均值趋于零，方差趋于 \(k_{\mathrm B}T/\kappa\)。于是

\[
\begin{aligned}
p_{\mathrm{st}}(x)
&=\frac{1}{\sqrt{2\pi\sigma_{\mathrm{st}}^2}}
\exp\left(-\frac{x^2}{2\sigma_{\mathrm{st}}^2}\right)\\
&=\sqrt{\frac{\kappa}{2\pi k_{\mathrm B}T}}
\exp\left(-\frac{\kappa x^2}{2k_{\mathrm B}T}\right).
\end{aligned}
\]

对于谐势 \(V(x)=\kappa x^2/2\)，这恰好就是熟悉的 Boltzmann 分布 \(p_{\mathrm{eq}}(x)\propto e^{-V(x)/(k_{\mathrm B}T)}\)。

其次，摩擦系数 \(\zeta\) 决定的是粒子接近平衡状态的速度，而不是最终的涨落幅度。增大摩擦，一方面会减慢随机扩散，另一方面也会减慢回复运动。对于固定的温度和势能刚度，扩散系数 \(D=k_{\mathrm B}T/\zeta\) 和回复率 \(\gamma=\kappa/\zeta\) 都随摩擦系数以相同的比例减小，因此稳态方差 \(D/\gamma=k_{\mathrm B}T/\kappa\) 保持不变，但弛豫时间 \(\tau=\zeta/\kappa\) 会变长。

平衡位置分布告诉我们粒子通常出现在什么位置，却无法反映粒子运动的快慢。即使两个系统具有完全相同的平衡分布，它们的弛豫时间也可能相差很大。要区分这样的动力学差异，仅仅统计粒子的位置是不够的，还需要研究不同时刻的位置之间存在怎样的关联。

OU 过程还可以用来描述其他具有线性回复机制的涨落。例如，在前面的自由惯性 Langevin 方程中，速度受到线性阻尼和随机力的共同作用，因此速度本身也是一个 OU 过程。此时对应的回复率为 $\gamma=\zeta/m$，噪声强度参数为 $D=\zeta k_{\mathrm B}T/m^2$，稳态方差 $D/\gamma=k_{\mathrm B}T/m$ 正好给出热平衡速度涨落。位置在谐势中的回复与速度在摩擦作用下的弛豫，实际上具有相同的随机动力学结构。

下面的图展示了 OU 过程的典型行为。尽管从固定初态出发的平均位置逐渐返回平衡点，单条轨迹仍然围绕平衡点持续涨落；位置方差则逐渐达到稳定值。

![OU 过程的随机轨迹、方差、时间关联与响应](../figures/01_ou_fluctuations.svg)

*Fig. 2. 取 $\gamma=1,D=0.4$。左上图展示从 $x_0=2$ 出发的随机轨迹及其集合平均，阴影表示平均值两侧一个标准差的范围；右上图展示位置方差逐渐趋于 $D/\gamma$ 的过程。下方分别给出稳态位置关联函数与受外力扰动后的响应，它们将在下一节讨论。轨迹采用 OU 过程的精确转移分布采样。*

下面先计算这些轨迹的时间关联，再比较微小外力引起的平均响应。

<span id="18-时间关联频谱与线性响应"></span>

## 1.6 从自发涨落到线性响应 {#correlation-response}

### 稳态涨落的时间关联 {#stationary-correlations}

对处于热平衡的 OU 过程，平均位置为零，定义位置自关联函数

\[
C(\tau)=\langle x(t+\tau)x(t)\rangle_{\mathrm{st}}.
\]

它衡量此刻的位置偏离在时间 $\tau$ 后还能保留多少；若两个时刻的位置独立，由于均值为零，关联也为零。平稳性使 $C$ 只依赖时间差。利用上一节的精确解，对 $\tau\geq0$ 有

\[
x(t+\tau)=e^{-\gamma\tau}x(t)
+\sqrt{2D}\int_t^{t+\tau}e^{-\gamma(t+\tau-s)}\,dW_s.
\]

其中，第一项来自时刻 \(t\) 的位置，第二项来自之后新产生的随机扰动。由于后续噪声与 \(x(t)\) 相互独立，将两边乘以 \(x(t)\) 并取平均，噪声项便不再贡献关联。再利用稳态方差 \(\langle x^2\rangle_{\mathrm{st}}=D/\gamma\)，得到

\[
C(\tau)
=e^{-\gamma\tau}\langle x^2\rangle_{\mathrm{st}}
=\frac{D}{\gamma}e^{-\gamma\tau},
\qquad \tau\geq0.
\]

对于负的时间差，利用关联函数的对称性，便有 \(C(\tau)=(D/\gamma)e^{-\gamma|\tau|}\)。因此，OU 过程的位置关联按照指数规律衰减，其特征时间正是前面得到的回复时间 \(\tau=\gamma^{-1}\)。在短于这个时间的尺度上，粒子的位置仍然保留着此前运动的信息；经过几个回复时间后，原来的偏离基本被新的随机运动所取代。

在实验和数值模拟中，时间关联还会影响统计量的测量精度。如果相邻两次采样相隔的时间远短于关联时间，那么这些数据并不是独立的。即使记录了很多个数据点，它们也可能主要反映同一段持续的涨落。要获得更多近似独立的样本，就需要延长观测时间，使记录包含足够多个关联时间。

除了直接观察时间关联，也可以研究它在频率空间中的表现。对 \(C(\tau)\) 作 Fourier 变换，定义功率谱（power spectral density）

\[
S(\omega)=\int_{-\infty}^{\infty}C(\tau)e^{i\omega\tau}\,d\tau,
\]

对于 OU 过程可以得到

\[
S(\omega)=\frac{2D}{\gamma^2+\omega^2}.
\]

这就是 Lorentzian 频谱。虽然环境提供的是白噪声，粒子的位置涨落却并不是各个频率都具有相同强度。回复动力学使位置主要表现出低频涨落，而高于特征频率 \(\gamma\) 的快速变化受到抑制。时间域中的指数关联和频率域中的 Lorentzian 频谱，实际上描述的是同一套动力学。

### 自发涨落与外力响应 {#fluctuation-response}

到这里，我们一直在研究没有额外扰动时，粒子自身的热涨落。现在考虑另一个问题：如果对已经处于平衡态的粒子施加一个很小的外力，它的平均位置会怎样变化？这对应着两种不同的测量方式。时间关联函数研究粒子在自然热运动中，不同时刻的位置如何相关；响应函数则研究人为施加一个扰动以后，粒子的平均运动会发生多大的变化。

对于谐势中的过阻尼粒子，施加一个随时间变化的小外力 \(f(t)\)，Langevin 方程变成

\[
dx=[-\gamma x+\mu f(t)]\,dt+\sqrt{2D}\,dW_t,
\]

其中 \(\mu=1/\zeta\) 是迁移率。由于噪声的平均值为零，粒子的平均位置满足

\[
\frac{d\langle x\rangle}{dt}
=-\gamma\langle x\rangle+\mu f(t).
\]

假设系统最初处于平衡态，并从过去某个时刻开始受到外力。求解这个线性方程，可以得到外力引起的平均位置变化：

\[
\delta\langle x(t)\rangle
=\mu\int_{-\infty}^{t}e^{-\gamma(t-s)}f(s)\,ds.
\]

这个结果与之前的 OU 精确解非常相似。发生在时刻 \(s\) 的外力会改变粒子的平均位置，但这次扰动的影响也会随着时间按照 \(e^{-\gamma(t-s)}\) 衰减。因此，我们可以定义线性响应函数

\[
R_f(\tau)=\mu\Theta(\tau)e^{-\gamma\tau},
\]

其中 \(\Theta(\tau)\) 是阶跃函数，表示粒子只能在受到外力之后作出响应。响应函数描述了某一时刻施加的微小外力，会怎样影响之后的平均位置。值得注意的是，刚才求得的自发关联函数与响应函数具有相同的指数衰减形式。对于平衡热浴，它们的关系还不止于此。利用 \(C(\tau)=(D/\gamma)e^{-\gamma\tau}\)，以及前面得到的 Einstein 关系 \(D=\mu k_{\mathrm B}T\)，可以直接验证

\[
\boxed{
R_f(\tau)
=-\frac{1}{k_{\mathrm B}T}\frac{dC(\tau)}{d\tau},
\qquad \tau>0.
}
\]

这就是涨落–耗散定理（fluctuation–dissipation theorem）在当前模型中的一个具体形式。[附录 B](#fdt-extension)给出不依赖谐势的平衡推导，并将它写成频谱与耗散响应之间的关系。

它告诉我们，在热平衡条件下，即使不对系统施加任何外力，也可以通过测量粒子自发涨落的时间关联，推断它受到微小外力时将如何响应。关联函数与响应函数原本对应两种不同的实验，却因为热平衡条件而受到同一个关系约束。

作为一个简单的例子，如果施加恒定外力 \(f\)，粒子的平均位置最终会移动到 \(\langle x\rangle_{\mathrm{st}}=f/\kappa\)。这正是弹簧在恒定外力下的新平衡位置。与此同时，由于势能的曲率没有改变，位置涨落的稳态方差仍然是 \(k_{\mathrm B}T/\kappa\)。外力改变了平均位置，却没有改变围绕平均位置的热涨落宽度。在更一般的非平衡系统中，自发关联与外力响应未必满足上述平衡关系。因此，对关联函数和响应函数进行比较，也是研究系统是否具有平衡动力学性质的一种重要方法。

通过 Langevin 描述，我们已经能够从单条粒子的随机运动出发，计算扩散、弛豫、时间关联和外力响应。接下来如果希望直接研究许多次随机运动所形成的概率分布，就需要转向 Fokker–Planck 描述。

## 附录 A：反复用到的计算工具 {#calculation-tools}

正文中的速度弛豫、OU 解、方差与均方位移，都用到了下面两个计算步骤。先从普通积分说明，再写出白噪声对应的随机积分形式。

### A.1 一阶线性微分方程的积分因子法 {#integrating-factor}

考虑初值问题

\[
\frac{dy}{dt}+a(t)y=b(t),\qquad y(t_0)=y_0.
\]

我们希望把左边写成一个乘积的导数。将方程乘以一个只依赖时间的函数 $\mathcal I(t)$，根据乘积求导法则，

\[
\frac{d}{dt}[\mathcal I(t)y(t)]
=\mathcal I(t)\frac{dy}{dt}+\mathcal I'(t)y(t).
\]

因此只要选取 $\mathcal I'=a\mathcal I$，就能把原方程左边的两项合起来。取 $\mathcal I(t_0)=1$，得到积分因子

\[
\mathcal I(t)=\exp\left[\int_{t_0}^t a(u)\,du\right],
\qquad
\frac{d}{dt}[\mathcal I(t)y(t)]=\mathcal I(t)b(t).
\]

从 $t_0$ 积分到 $t$，再除以 $\mathcal I(t)$，便得到

\[
\boxed{
y(t)=y_0\exp\left[-\int_{t_0}^t a(u)\,du\right]
+\int_{t_0}^t\exp\left[-\int_s^t a(u)\,du\right]b(s)\,ds.
}
\]

这里用 $\mathcal I$ 表示积分因子，以免与正文表示迁移率的 $\mu$ 混淆。若 $a(t)=a$ 为常数，公式简化为

\[
y(t)=y_0e^{-a(t-t_0)}+\int_{t_0}^t e^{-a(t-s)}b(s)\,ds.
\]

对正文的自由粒子，$\dot v+v/\tau_m=\eta/m$，因此取 $a=1/\tau_m$、$b=\eta/m$，便得到速度解。对受外力的平均位置，取 $a=\gamma$、$b=\mu f(t)$，就得到线性响应的卷积表达式；令初始时刻趋于 $-\infty$，衰减的初值项也随之消失。

理想白噪声不是普通的时间函数，更严格的写法是

\[
dy=[-a(t)y+b(t)]\,dt+\sigma(t)\,dW_t,
\]

其中 $a,b,\sigma$ 为确定性函数。由于积分因子没有随机增量，Itô 乘积公式中不会多出二次变差项，因此仍有

\[
d(\mathcal I y)=\mathcal I b\,dt+\mathcal I\sigma\,dW_t.
\]

积分以后，通解只需在确定性解后再加上

\[
\int_{t_0}^t\exp\left[-\int_s^t a(u)\,du\right]\sigma(s)\,dW_s.
\]

令 $a=\gamma$、$b=0$、$\sigma=\sqrt{2D}$，就得到 OU 的精确解。这说明指数核的来源是线性回复项；随机性改变的是被指数核累积的输入。

### A.2 积分的平方为什么变成二重积分 {#squared-integral}

先考虑一个普通的可积函数 $g$。平方就是两个相同积分相乘，但两个积分各有自己的积分变量：

\[
\begin{aligned}
\left[\int_0^t g(s)\,ds\right]^2
&=\left[\int_0^t g(s)\,ds\right]
  \left[\int_0^t g(u)\,du\right]\\
&=\int_0^t ds\int_0^t du\,g(s)g(u).
\end{aligned}
\]

用离散求和看得更直接：$[\sum_i g_i\Delta t]^2=\sum_{i,j}g_i g_j(\Delta t)^2$。平方展开包含每一对时刻的乘积，而不只是 $i=j$ 的项，因此一般不能写成 $\int_0^t g(s)^2\,ds$。

对随机过程 $X$，定义 $Y(t)=\int_0^t K(s)X(s)\,ds$，其中 $K$ 为确定性权重。在二阶矩有限、可交换平均与积分的条件下，

\[
\begin{aligned}
\langle Y(t)^2\rangle
&=\int_0^t ds\int_0^t du\,
K(s)K(u)\langle X(s)X(u)\rangle,\\
\operatorname{Var}Y(t)
&=\int_0^t ds\int_0^t du\,
K(s)K(u)\operatorname{Cov}[X(s),X(u)].
\end{aligned}
\]

第二行通过对 $X$ 去均值得到。若过程均值非零，二阶矩与方差必须区分；正文中自由平衡粒子的平均速度为零，因此均方位移就是位移方差。

**白噪声怎样把重积分缩成单积分？** 对 $\langle\eta(s)\eta(u)\rangle=A\delta(s-u)$，在随机积分意义下有

\[
\begin{aligned}
\left\langle\left[\int_0^t K(s)\eta(s)\,ds\right]^2\right\rangle
&=A\int_0^t ds\int_0^t du\,K(s)K(u)\delta(s-u)\\
&=A\int_0^t K(s)^2\,ds.
\end{aligned}
\]

这里的 $\delta$ 将两个时刻约束为同一时刻。严格地说，$\eta(s)\,ds=\sqrt A\,dW_s$，上式就是确定性被积函数的 Itô 等距关系。它也可以从独立增量推出：$\langle\Delta W_i\Delta W_j\rangle=\delta_{ij}\Delta t$，因此平方求和中只留下同一时间格的方差。

例如，速度解的随机部分权重为 $K(s)=m^{-1}e^{-(t-s)/\tau_m}$，因此

\[
\operatorname{Var}\left[\frac1m\int_0^t e^{-(t-s)/\tau_m}\eta(s)\,ds\right]
=\frac A{m^2}\int_0^t e^{-2(t-s)/\tau_m}\,ds
=\frac A{2\zeta m}(1-e^{-2t/\tau_m}).
\]

如果初速度与后续噪声独立，交叉协方差为零，再加上 $e^{-2t/\tau_m}\operatorname{Var}(v_0)$，就回到正文的速度方差。同理，OU 随机积分的方差为 $2D\int_0^t e^{-2\gamma(t-s)}ds$。

**平稳关联怎样再化成 $2\int_0^t(t-\tau)C(\tau)d\tau$？** 对均值为零的平稳标量速度，$\langle v(s)v(u)\rangle=C_v(s-u)$，而交换两个经典随机变量给出 $C_v(-\tau)=C_v(\tau)$。对普通可积的速度关联函数，正方形积分区域可沿对角线分成两个贡献相等的三角形：

\[
\begin{aligned}
\left\langle\left[\int_0^t v(s)\,ds\right]^2\right\rangle
&=\int_0^t ds\int_0^t du\,C_v(s-u)\\
&=2\int_0^t ds\int_0^s du\,C_v(s-u)\\
&=2\int_0^t ds\int_0^s d\tau\,C_v(\tau)\\
&=2\int_0^t d\tau\int_\tau^t ds\,C_v(\tau)\\
&=2\int_0^t(t-\tau)C_v(\tau)\,d\tau.
\end{aligned}
\]

第三行令 $\tau=s-u$；第四行交换积分次序，区域 $0\leq\tau\leq s\leq t$ 变成先固定 $\tau$，再让 $s$ 从 $\tau$ 取到 $t$。所以权重 $t-\tau$ 表示观测窗口中具有同一时间差的时间对所占的长度。这里速度关联在对角线上没有奇异项；若直接处理理想白噪声，应使用前面的 $\delta$ 积分或 Itô 等距关系。

<span id="附录-a怎样模拟怎样判断结果可信"></span>
<span id="附录-a数值小步怎样保留正确涨落"></span>

### A.3 数值小步怎样保留正确涨落 {#numerical-checks}

对加性噪声，记 $x_n=x(t_n)$、$t_n=nh$。Euler–Maruyama 更新为 $x_{n+1}=x_n+b(x_n)h+\sqrt{2Dh}Z_n$，其中每步的 $Z_n$ 独立服从标准高斯。随机项必须按 $\sqrt h$ 缩放，才能使固定时长内累积的方差保持有限。

OU 可用精确更新作为基准：

$$
x_{n+1}=e^{-\gamma h}x_n
+\sqrt{\frac D\gamma(1-e^{-2\gamma h})}\,Z_n.
$$

Euler 版本为 $x_{n+1}=(1-\gamma h)x_n+\sqrt{2Dh}Z_n$，稳态方差

$$
\sigma_h^2=\frac{D}{\gamma(1-\gamma h/2)},\qquad0\lt\gamma h\lt2.
$$

稳定并不等于精确：接近稳定边界时，方差会显著偏大。比较减半步长后的宽度和关联时间，比判断一条轨迹“看起来会抖”更有判别力。

## 附录 B：拓展——涨落耗散关系 {#fdt-extension}

正文用谐势中的 OU 过程验证了 $R_f(\tau)=-C'(\tau)/(k_{\mathrm B}T)$。这个关系并不要求回复力一定是线性的；它来自平衡分布与微扰方式之间的联系。以下讨论经典热平衡系统，外力通过能量项 $-f(t)x$ 与位置耦合，并假设位置方差有限、关联最终衰减。线性响应的一般背景可参阅 [Kubo（1966）](https://www.zib.de/userpage/donati/stochastics2023/06/references/Kubo1966.pdf)。

### B.1 从撤去小外力后的弛豫推导

记 $\beta=1/(k_{\mathrm B}T)$、$\delta x=x-\langle x\rangle_{\mathrm{eq}}$，以及

\[
C_x(\tau)=\langle\delta x(\tau)\delta x(0)\rangle_{\mathrm{eq}}.
\]

设在 $t<0$ 时施加恒定的小外力 $f_0$，并等待系统在扰动后的能量 $H_0-f_0x$ 下达到平衡。在 $t=0$ 突然撤去外力，此后按原来的无扰动动力学演化。用 $z$ 表示完整状态，撤去外力时的分布为

\[
\begin{aligned}
p_{f_0}(z)
&=\frac{p_{\mathrm{eq}}(z)e^{\beta f_0x(z)}}
{\langle e^{\beta f_0x}\rangle_{\mathrm{eq}}}\\
&=p_{\mathrm{eq}}(z)[1+\beta f_0\delta x(z)]+O(f_0^2).
\end{aligned}
\]

初始分布的一阶偏离正比于初始位置涨落。将这个分布按无扰动动力学推进到时刻 $t>0$，便得到平均位置的弛豫

\[
\delta\langle x(t)\rangle=\beta f_0 C_x(t)+O(f_0^2).
\]

另一方面，按照线性响应的定义，同一次“先施加、后撤去”的实验满足

\[
\begin{aligned}
\delta\langle x(t)\rangle
&=f_0\int_{-\infty}^0 R_f(t-s)\,ds\\
&=f_0\int_t^\infty R_f(\tau)\,d\tau.
\end{aligned}
\]

比较两式并对 $t$ 求导，得到

\[
\boxed{R_f(t)=-\beta\frac{dC_x(t)}{dt},\qquad t>0.}
\]

这个推导说明，自发涨落之所以能预测外力响应，是因为一个小的平衡外力恰好以 $\delta x$ 为权重改变初始状态的概率。只要无扰动系统及其受力扰动属于同一温度下的平衡动力学，这个结论就可用于一般的束缚势。若改用其他观测量，响应所对应的关联函数也要与外场在能量中耦合的变量相匹配。

若 $C_x(\infty)=0$，积分还给出静态易感率

\[
\chi(0)=\int_0^\infty R_f(t)\,dt
=\beta C_x(0)=\frac{\operatorname{Var}_{\mathrm{eq}}x}{k_{\mathrm B}T}.
\]

在谐势中，$C_x(0)=k_{\mathrm B}T/\kappa$，所以 $\chi(0)=1/\kappa$，与恒力导致的位置移动 $f_0/\kappa$ 一致。

### B.2 频率空间中，为什么叫“耗散”

沿用正文的双边 Fourier 变换约定，定义

\[
S_x(\omega)=\int_{-\infty}^{\infty}C_x(t)e^{i\omega t}\,dt,
\qquad
\chi(\omega)=\int_0^\infty R_f(t)e^{i\omega t}\,dt.
\]

利用 $C_x(t)$ 的偶对称性及分部积分，可得

\[
\begin{aligned}
\operatorname{Im}\chi(\omega)
&=-\beta\int_0^\infty C_x'(t)\sin(\omega t)\,dt\\
&=\beta\omega\int_0^\infty C_x(t)\cos(\omega t)\,dt\\
&=\frac{\beta\omega}{2}S_x(\omega).
\end{aligned}
\]

于是经典涨落耗散关系也可以写成

\[
\boxed{S_x(\omega)=\frac{2k_{\mathrm B}T}{\omega}
\operatorname{Im}\chi(\omega),\qquad\omega\ne0.}
\]

若零频极限存在，可取 $\omega\to0$ 理解。系数 $2$ 对应这里的双边功率谱定义。响应的虚部与耗散的联系及其经典极限，见 [Tong 的线性响应讲义，第 4.2.2、4.3.2 节](https://www.damtp.cam.ac.uk/user/tong/kintheory/four.pdf)。

例如，施加 $f(t)=f_0\cos\omega t$，经过瞬态后，外力平均输入功率为 $\overline{f(t)\,d\langle x(t)\rangle/dt}=f_0^2\omega\operatorname{Im}\chi(\omega)/2$。虚部决定每个周期中净输入并耗散掉的能量，因此称为耗散响应。对正文的 OU 模型，

\[
\chi(\omega)=\frac{\mu}{\gamma-i\omega},
\qquad
S_x(\omega)=\frac{2\mu k_{\mathrm B}T}{\gamma^2+\omega^2}
=\frac{2D}{\gamma^2+\omega^2},
\]

恰好回到正文的 Lorentzian 频谱。正文的噪声强度 $A=2\zeta k_{\mathrm B}T$ 约束热浴的随机力与摩擦；这里的关系约束粒子的自发关联与受力响应。两种形式表达的是同一个平衡热浴对涨落与耗散的限制。对主动驱动等非平衡环境，不能仅凭稳态或指数关联就套用这些平衡公式。

## 附录 C：拓展——Green–Kubo 关系 {#green-kubo-extension}

Green–Kubo 关系把输运系数写成平衡态中相应运动或流的时间关联积分。对粒子扩散，先由位移的统计得到扩散系数，再由小外力响应得到迁移率，可以看清它与涨落耗散关系、Einstein 关系之间的联系。

### C.1 从均方位移到扩散系数

考虑零平均速度的平稳过程。由[附录 A.2](#squared-integral)，

\[
\frac{\langle\Delta x(t)^2\rangle}{2t}
=\int_0^t\left(1-\frac{\tau}{t}\right)C_v(\tau)\,d\tau.
\]

若速度二阶矩有限，且 $\int_0^\infty|C_v(\tau)|d\tau<\infty$，便可取长时间极限，得到

\[
\boxed{D=\lim_{t\to\infty}\frac{\langle\Delta x(t)^2\rangle}{2t}
=\int_0^\infty\langle v(\tau)v(0)\rangle_{\mathrm{eq}}\,d\tau.}
\]

对这个积分恒等式，平稳性和可积的关联已经足够；将它进一步与热力学响应联系起来，才需要平衡条件。如果存在平均漂移，应先使用 $\delta v=v-\langle v\rangle$，并从位移中减去平均漂移。若关联尾部不可积，或速度二阶矩不存在，上述有限扩散系数的推导便不能直接使用，可能出现第四篇讨论的反常输运。

### C.2 从冲量响应到迁移率

对正文的惯性粒子，在 $t=0$ 施加一个很短的力脉冲，总冲量为 $I$，速度便瞬间增加 $I/m$。平衡速度密度 $p_{\mathrm{eq}}(v)\propto e^{-\beta mv^2/2}$ 因而变成

\[
p_I(v)=p_{\mathrm{eq}}(v-I/m)
=p_{\mathrm{eq}}(v)[1+\beta I v]+O(I^2).
\]

脉冲之后按无外力动力学演化，将这个初始分布代入平均速度，得到

\[
\delta\langle v(t)\rangle
=\beta I\langle v(t)v(0)\rangle_{\mathrm{eq}}+O(I^2).
\]

所以速度对外力的响应核是 $R_{vf}(t)=\beta C_v(t)$（$t>0$）。对持续施加的小恒力，将历次冲量的响应相加；若响应可积，就有

\[
\boxed{
\mu=\lim_{f\to0}\frac{\langle v\rangle_{\mathrm{st},f}}{f}
=\int_0^\infty R_{vf}(t)\,dt
=\frac1{k_{\mathrm B}T}\int_0^\infty C_v(t)\,dt.
}
\]

与上一节的扩散系数比较，便得到 $D=\mu k_{\mathrm B}T$。在正文的模型中，$C_v(t)=(k_{\mathrm B}T/m)e^{-t/\tau_m}$，因此 $D=(k_{\mathrm B}T/m)\tau_m$、$\mu=\tau_m/m=1/\zeta$。平衡速度关联不仅记录粒子怎样忘记初速度，也决定它受小外力时能产生多大的持续漂移。

### C.3 多维扩散与其他输运系数

在 $d$ 维各向同性系统中，$\langle|\Delta\boldsymbol r(t)|^2\rangle\sim2dDt$，因此

\[
D=\frac1d\int_0^\infty
\langle\boldsymbol v(t)\cdot\boldsymbol v(0)\rangle_{\mathrm{eq}}\,dt.
\]

若系统各向异性，则用 $\langle\Delta r_i\Delta r_j\rangle\sim2D_{ij}t$ 定义对称扩散张量，

\[
D_{ij}=\frac12\int_0^\infty
\left[\langle\delta v_i(t)\delta v_j(0)\rangle
+\langle\delta v_j(t)\delta v_i(0)\rangle\right]_{\mathrm{eq}}dt.
\]

类似结构也出现在其他输运现象中。例如，对体积 $V$ 的平衡流体，若 $P_{xy}$ 表示瞬时的体积平均剪切应力分量，则剪切黏度可以写成

\[
\eta_{\mathrm{sh}}=\frac{V}{k_{\mathrm B}T}
\int_0^\infty
\langle\delta P_{xy}(t)\delta P_{xy}(0)\rangle_{\mathrm{eq}}\,dt.
\]

这里的 $\eta_{\mathrm{sh}}$ 是黏度，不是正文的随机力 $\eta(t)$。应力记录动量输运，正如速度记录粒子位移；黏度的平衡关联计算及应力归一化示例见 [LAMMPS 的黏度说明](https://docs.lammps.org/Howto_viscosity.html)。不同输运问题的流及其共轭驱动具有不同量纲，前面的体积和温度因子需要随定义确定。

在数值计算中，可以检查截断积分 $D(t_c)=\int_0^{t_c}C_v(t)dt$ 是否在关联衰减后出现稳定区间，并与均方位移的长时斜率比较。更长的积分上限也会累积更多关联尾部的采样噪声，因此应同时增加独立轨迹或观测时长来判断收敛。

## 参考文献 {#references}

以下文献对应本文的 Brownian 运动、环境约化、OU 过程及涨落响应推导；每条后附与正文的联系。

1. Einstein, A. (1905). [Über die von der molekularkinetischen Theorie der Wärme geforderte Bewegung von in ruhenden Flüssigkeiten suspendierten Teilchen](https://doi.org/10.1002/andp.19053220806). *Annalen der Physik* **322**, 549–560. Brownian 扩散与迁移率之间的 Einstein 关系，对应第 1.3 节。
2. Uhlenbeck, G. E. & Ornstein, L. S. (1930). [On the Theory of the Brownian Motion](https://doi.org/10.1103/PhysRev.36.823). *Physical Review* **36**, 823–841. 速度弛豫、位移统计与谐势中 Brownian 运动的经典计算，对应第 1.3、1.5–1.6 节。
3. Kubo, R. (1966). [The fluctuation-dissipation theorem](https://doi.org/10.1088/0034-4885/29/1/306). *Reports on Progress in Physics* **29**, 255–284. 广义 Langevin 方程、摩擦记忆、线性响应与平衡涨落之间的联系，对应第 1.1–1.2、1.6 节及附录 B–C；另有[开放全文](https://www.zib.de/userpage/donati/stochastics2023/06/references/Kubo1966.pdf)。
4. Di Pasquale, N., Hudson, T., Icardi, M., Rovigatti, L. & Spinaci, M. (2023). [A systematic analysis of the memory term in coarse-grained models: The case of the Markovian approximation](https://doi.org/10.1017/S0956792522000158). *European Journal of Applied Mathematics* **34**, 326–345. 分析以瞬时作用近似记忆项的条件与误差，并讨论空间依赖，对应第 1.1 节。
5. Cugliandolo, L. F. (2011). [Out of equilibrium dynamics of complex systems](https://giamarchi.unige.ch/wp-content/php_code/notes/cuso/Cugliandolo/DynamicsComplexSystems-lectures1-4.pdf). Cours de troisième cycle de la Suisse Romande 课程讲义，第 2.3–2.4 节。由热浴模型建立 Langevin 描述，并讨论记忆、过阻尼极限与谐势弛豫，可结合第 1.1、1.4–1.5 节阅读。
6. Tong, D. (2012). [Kinetic Theory](https://www.damtp.cam.ac.uk/user/tong/kintheory/). University of Cambridge 课程讲义，第 3–4 章。第 3 章介绍 Langevin、扩散与 Einstein 关系，第 4 章介绍响应函数、耗散和涨落耗散定理，可用于串联正文与附录 B–C。
7. LAMMPS Developers. [Calculate viscosity](https://docs.lammps.org/Howto_viscosity.html). *LAMMPS Documentation*. 通过平衡剪切应力的自关联计算黏度的 Green–Kubo 示例，对应附录 C.3。

**下一篇：**[Fokker–Planck Dynamics：从随机轨迹到概率流](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/fokker_planck/)。
