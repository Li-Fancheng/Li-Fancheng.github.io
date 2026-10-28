把同一个微小粒子在相同条件下释放很多次，每次得到的轨迹都不同。确定性力仍然约束运动方向，液体分子的碰撞却不断改变粒子的瞬时运动。理论要解释的，既包括为什么一次运动不可重复，也包括为什么许多次实验的扩散速度、涨落幅度和记忆时间可以精确预测。

本篇从阻尼和随机驱动的物理来源开始，再讨论随机轨道为何需要不同的微积分。OU 过程会贯穿计算：回复作用不断消除偏离，噪声不断产生偏离，两者共同形成有限的稳态涨落。

## 1.1 环境自由度怎样变成阻尼和噪声

### 从惯性运动开始

设粒子质量为 $m$，位置为 $X$，速度为 $V$，真实势能为 $V_{\mathrm{pot}}(X)$。最简的热浴模型为

$$
dX_t=V_t\,dt,
\qquad
m\,dV_t=[-\zeta V_t-V_{\mathrm{pot}}'(X_t)]\,dt
+\sqrt{2\zeta k_{\mathrm B}T}\,dW_t.
$$

阻力 $-\zeta V$ 描述粒子运动带给环境的平均反作用；噪声描述环境在相同宏观条件下仍然存在的微观涨落。这里的关系 $\text{噪声振幅}^2=2\zeta k_{\mathrm B}T$ 是这个平衡热浴模型的涨落耗散关系。它保证自由速度的稳态方差为 $k_{\mathrm B}T/m$，不能对任意主动驱动或非热环境直接套用。

这种写法把大量环境变量的影响归入一个平均项和一个随机项。它依赖两个物理判断。首先，所保留的粒子变量比被消去的环境变量慢；其次，环境对先前运动的记忆在观察时间分辨率内已经衰减。如果环境关联时间为 $\tau_{\mathrm{env}}$，白噪声描述要求时间步长显著大于它，同时显著小于所研究的弛豫时间。若没有这段可用的尺度窗口，带记忆的摩擦核和有色噪声就可能必须保留。

对于自由粒子，速度弛豫时间为 $\tau_m=m/\zeta$，稳态速度关联为

$$
\langle V_tV_0\rangle=\frac{k_{\mathrm B}T}{m}e^{-|t|/\tau_m}.
$$

对速度积分得到位置的均方位移

$$
\mathbb E[(X_t-X_0)^2]
=2D_x\left[t-\tau_m(1-e^{-t/\tau_m})\right],
\qquad D_x=\frac{k_{\mathrm B}T}{\zeta}.
$$

短时 $t\ll\tau_m$，展开指数得到 $(k_{\mathrm B}T/m)t^2$：粒子保留初始速度，运动近似弹道。长时 $t\gg\tau_m$，则得到 $2D_xt$：速度记忆已经丢失，位置扩散。因而过阻尼 Brownian 轨道的粗糙性质对应消去惯性以后的尺度，不能无限外推到真实粒子的任意短时运动。[Tong 的随机过程讲义](https://www.damtp.cam.ac.uk/user/tong/kintheory/three.pdf)也从速度阻尼与热浴涨落出发讨论这一过渡。

### 过阻尼极限保留了什么

在位置演化比 $\tau_m$ 慢、摩擦系数为常数时，消去快速速度变量得到

$$
dX_t=-\mu V_{\mathrm{pot}}'(X_t)\,dt+\sqrt{2D_x}\,dW_t,
\qquad \mu=\zeta^{-1}.
$$

以后用 $U=\mu V_{\mathrm{pot}}$、$D=D_x$，把它写成 $dX=-U'\,dt+\sqrt{2D}\,dW$。这解释了为什么 $U/D=V_{\mathrm{pot}}/(k_{\mathrm B}T)$。对一般有效随机模型，$b(x)$ 与 $D$ 可以独立指定；只有补上热浴来源，才有相应的温度解释。

白噪声常记为 $\xi(t)=\dot W_t$，满足

$$
\langle\xi(t)\rangle=0,\qquad
\langle\xi(t)\xi(s)\rangle=\delta(t-s).
$$

这里的导数和相关函数是广义函数意义下的记号。$\xi(t)$ 不是每个时刻都有有限方差的普通随机变量；真正可以直接使用的是有限时间积分 $W_{t+h}-W_t$。Langevin 方程的可靠含义也应从它的积分形式读取。

## 1.2 随机轨道的局部性质为什么改变微积分

标准 Wiener 过程从 $W_0=0$ 出发，具有独立的高斯增量：

$$
W_{t+h}-W_t\sim\mathcal N(0,h).
$$

增量的典型大小是 $\sqrt h$。对于可微轨道，位置增量是 $O(h)$；对于 Brownian 轨道，噪声增量在短时反而比漂移增量更大。它的差商典型幅度随 $h^{-1/2}$ 增长，轨道几乎处处不可微。把 Brownian 位置直接代入普通链式法则，会漏掉累积起来仍然有限的二阶贡献。

取时间分割 $0=t_0\lt\cdots\lt t_N=t$，有

$$
\sum_n(W_{t_{n+1}}-W_{t_n})^2\longrightarrow t
$$

（在分割变细时依概率收敛）。这称为二次变差。它不是说每一个随机增量的平方都恰好等于时间步长，而是说大量平方增量的和有确定的极限。随机微积分里的 $(dW)^2=dt$ 是对这个事实的简写。

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

例如自由扩散 $dX=\sqrt{2D}\,dW$ 的平均位置不动，但

$$
d(X^2)=2D\,dt+2\sqrt{2D}X\,dW,
\qquad \frac{d}{dt}\langle X^2\rangle=2D.
$$

噪声平均为零，并不意味着它对二阶量没有平均效应。扩散展宽正来自独立随机增量的平方累积。

## 1.3 随机积分怎样记录噪声与状态之间的先后关系

对于普通积分，区间左端或中点取值在极限中通常没有差别。对于随机积分，被积函数也会随同一个噪声变化，这种差别能够留下有限贡献。

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

对于一维乘性噪声，相同过程的两种表达为

$$
dX=a_{\mathrm S}(X)\,dt+\sigma(X)\circ dW
$$

和

$$
dX=\left[a_{\mathrm S}(X)+\frac12\sigma(X)\sigma'(X)\right]dt
+\sigma(X)\,dW.
$$

选择约定时要同时转换漂移。若 $\sigma$ 为常数，两者的方程相同；若 $\sigma$ 依赖状态，保持同一漂移再改变积分约定就改变了模型。由平滑有色外力取某类短关联时间极限时常得到 Stratonovich 形式；由离散跳跃过程的条件增量得到扩散近似时常自然得到 Itô 形式。具体粗粒化和同时消去哪些变量，决定正确的极限，不能仅靠“物理噪声”四个字选择。

一个能直接看出差异的例子是几何 Brownian 运动：

$$
dX=aX\,dt+sX\,dW,\qquad X_0>0.
$$

对 $\log X$ 用 Itô 公式，得到

$$
X_t=X_0\exp\left[(a-s^2/2)t+sW_t\right].
$$

其均值为 $X_0e^{at}$，而中位数为 $X_0e^{(a-s^2/2)t}$。少量增长很大的轨迹可以支撑均值，即使典型轨迹正在减小。这说明对随机运动，“平均增长”和“典型增长”也需要分开。[Pavliotis 的教材](https://www.ma.imperial.ac.uk/~pavl/PavliotisBook.pdf)第 3 章给出了随机积分和转换公式的系统处理。

## 1.4 OU 过程：涨落如何被产生，又如何被消除

考虑

$$
dX_t=-\gamma X_t\,dt+\sqrt{2D}\,dW_t,
\qquad \gamma>0.
$$

它既可描述谐势中的过阻尼粒子，也可描述一般稳定固定点附近的线性涨落：若 $b(x_*)=0$ 且 $b'(x_*)=-\gamma$，在噪声振幅近似常数的局部区域，令 $X=x-x_*$ 即得到 OU。这个近似只覆盖典型小涨落，不自动适用于远离固定点的跨势垒运动。

### 精确解与有限时间传播

乘以积分因子 $e^{\gamma t}$，得到

$$
X_t=e^{-\gamma t}X_0
+\sqrt{2D}\int_0^t e^{-\gamma(t-s)}\,dW_s.
$$

若 $X_0=x_0$ 固定，噪声积分是高斯变量。由 Itô 等距，

$$
\mathbb E[X_t\mid x_0]=x_0e^{-\gamma t},\qquad
\operatorname{Var}(X_t\mid x_0)=\frac D\gamma(1-e^{-2\gamma t}).
$$

所以转移密度是

$$
p(x,t\mid x_0)=\frac{1}{\sqrt{2\pi v_t}}
\exp\left[-\frac{(x-x_0e^{-\gamma t})^2}{2v_t}\right],
\qquad v_t=\frac D\gamma(1-e^{-2\gamma t}).
$$

初始位置的影响逐渐衰减，后续噪声的影响逐渐积累。稳态并不是停止接受噪声，而是早先扰动被指数遗忘，新扰动持续补充。

### 矩方程中的机制

对 $X^2$ 使用 Itô 公式：

$$
d(X^2)=(-2\gamma X^2+2D)\,dt
+2\sqrt{2D}X\,dW.
$$

若初始条件与未来噪声独立且二阶矩有限，取平均得到

$$
\frac{dm}{dt}=-\gamma m,\qquad
\frac{dv}{dt}=-2\gamma v+2D,
$$

其中 $m=\mathbb E X$、$v=\operatorname{Var}X$。因此对一般初始方差 $v_0$，

$$
v(t)=v_0e^{-2\gamma t}+\frac D\gamma(1-e^{-2\gamma t}),
\qquad v_{\mathrm{ss}}=\frac D\gamma.
$$

$2D$ 是持续的涨落注入，$2\gamma v$ 是回复作用消除涨落的速率。增大 $D$ 会抬高稳态方差，增大 $\gamma$ 会压低它。均值的弛豫率为 $\gamma$，方差偏离稳态的弛豫率为 $2\gamma$；即使同一个线性系统，不同可观测量也可以看到不同时间尺度。

![OU 的单条轨道、方差、时间关联与响应](../figures/01_ou_fluctuations.svg)

*图 1｜取 $\gamma=1,D=0.4$。左上为从 $X_0=2$ 出发的样本轨道和集合均值，阴影为理论均值两侧一个标准差；右上显示方差向 $D/\gamma$ 饱和。下方将稳态关联与因果响应并列，二者有相同衰减率，却有不同定义。轨道用 OU 精确转移采样。*

## 1.5 时间关联、频谱与线性响应

### 系统保留信息多久

稳态初始条件满足 $\mathbb E X_0=0$、$\mathbb E X_0^2=D/\gamma$。对 $t\ge s$，把精确解从 $s$ 开始写出，后续噪声积分与 $X_s$ 独立，于是

$$
\boxed{C(t-s)=\langle X_tX_s\rangle_{\mathrm{ss}}
=\frac D\gamma e^{-\gamma|t-s|}.}
$$

相关时间 $1/\gamma$ 描述偏离的信息多久以后仍然可见。它同时决定时间平均的统计效率。令 $\overline X_{\mathcal T}=\mathcal T^{-1}\int_0^{\mathcal T}X_tdt$，则

$$
\operatorname{Var}(\overline X_{\mathcal T})
=\frac{2}{\mathcal T^2}\int_0^{\mathcal T}(\mathcal T-s)C(s)\,ds
\simeq\frac{2D}{\gamma^2\mathcal T}
\quad(\mathcal T\gg\gamma^{-1}).
$$

连续采样的点很多，不代表拥有同样多的独立样本。积分相关时间把“记录很密”和“获得新信息”区别开来。

采用双边 Fourier 变换 $S(\omega)=\int_{-\infty}^{\infty}e^{i\omega t}C(t)dt$，得到

$$
S(\omega)=\frac{2D}{\gamma^2+\omega^2}.
$$

OU 对低频扰动跟随得较充分，对快于回复时间的驱动响应较弱，形成 Lorentzian 频谱。测量稳态宽度和时间关联，分别给出 $D/\gamma$ 与 $\gamma$；只测稳态直方图则无法把二者分开。

### 响应需要一次受控扰动

在漂移中加入小外源 $h(t)$：

$$
dX_t=[-\gamma X_t+h(t)]dt+\sqrt{2D}\,dW_t.
$$

平均值的变化为

$$
\delta\langle X_t\rangle=\int_{-\infty}^t e^{-\gamma(t-s)}h(s)ds,
\qquad
R(t,s)=\frac{\delta\langle X_t\rangle}{\delta h(s)}
=\Theta(t-s)e^{-\gamma(t-s)}.
$$

关联是在未施加扰动时比较两个时刻；响应比较改变驱动前后的平均运动。二者在这个模型中满足 $R(\tau)=-C'(\tau)/D$（$\tau>0$），可以直接由上面的公式验证。若扰动是物理力 $f$，$h=\mu f$，热浴条件下 $R_f(\tau)=-C'(\tau)/(k_{\mathrm B}T)$。这种对应依赖相应的平衡条件和扰动耦合方式，不能仅凭两个量都指数衰减就推广到一般非平衡过程。

## 1.6 怎样模拟，怎样判断结果可信

对一般 Itô 方程，Euler–Maruyama 离散为

$$
X_{n+1}=X_n+b(X_n)h+\sigma(X_n)\sqrt h\,Z_n,
\qquad Z_n\overset{\mathrm{iid}}\sim\mathcal N(0,1).
$$

每步必须使用新的独立 $Z_n$，并保留 $\sqrt h$ 的缩放。把随机数乘以 $h$ 会使总噪声在步长趋零时消失。

对 OU，这成为 $X_{n+1}=(1-\gamma h)X_n+\sqrt{2Dh}Z_n$。其离散稳态方差为

$$
v_h=\frac{2Dh}{1-(1-\gamma h)^2}
=\frac{D}{\gamma(1-\gamma h/2)},\qquad 0\lt\gamma h\lt2.
$$

稳定并不意味着误差已经小；步长接近稳定性边缘时，方差会被严重高估。OU 本身还有精确一步更新

$$
X_{n+1}=e^{-\gamma h}X_n
+\sqrt{\frac D\gamma(1-e^{-2\gamma h})}\,Z_n.
$$

它适合用作数值基准。对非线性模型，至少应比较减半步长后的方差、关联或首达统计，并明确是在计算典型涨落还是罕见事件。一个能给出合理稳态直方图的步长，仍可能不能准确计算极少发生的跃迁。

这一篇建立了单条随机运动与可预测统计之间的联系：二次变差解释噪声为什么能持续注入涨落，OU 精确解解释回复如何限制涨落并消除记忆。接下来固定同一组漂移和噪声，追踪概率在状态空间里怎样流动。

**下一篇：**[Fokker–Planck Dynamics：概率输运、稳态与亚稳态](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/fokker_planck/)。
