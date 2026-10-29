本页保留路径积分的计算进阶：先用 Gaussian 泛函组织多时刻统计，再从离散动力学约束推导响应场表示，最后处理 OM 换元的 Jacobian。读者可在读完[路径积分正文](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/path_integrals/)以后，按需要选读；前三篇的物理主线不以此为前提。

本页固定加性白噪声、前点 Itô 规则；涉及中点表示时另行说明测度与换元。以下用 $j(t)$ 表示统计源，概率流仍记为 $J(x,t)$。计算普通关联和因果响应时，初始准备固定而终点自由；bridge 则必须使用条件化后的路径系综。

## C.1 路径观测量与 Gaussian 生成泛函 {#gaussian-functional}

[路径积分正文](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/path_integrals/)已经说明为什么需要路径系综。若希望计算整段记录的累计涨落，或者同时提取许多时刻的统计，就需要进一步组织这些平均。下面保留原有系综，改变我们从历史中读取的量。

### 从一个终点到沿途的累计量

两条历史即使有相同两端，也可能在中间长期停留在不同区域，因而给出不同的时间积分

$$
A_T=\int_0^T X_tdt.
$$

两时关联 $\langle X_tX_s\rangle$ 同样读取历史中的多个位置。对一个正确归一化的路径系综，这些量统一写成 $\langle\mathcal O[X]\rangle=\int\mathcal D x\,\mathcal O[x]e^{-S[x]}$，其中测度包含初始准备、终点条件及前面讨论的归一化因素。

要研究通常的关联和外力响应，以下固定初始分布 $p_0$，将终点自由积分，因而未扰动总路径权重为一。若研究 bridge 下的累计量，就改用前面除以终点传播核后的条件系综。观测量可以相同，平均所使用的历史集合必须说明。

### 生成泛函：一次组织所有多时刻统计

对一个随机变量，可以用矩母函数生成各阶矩。对随机历史，需要在不同时间分别读取位置，于是把单个参数推广为时间函数 $j(t)$：

$$
Z[j]=\left\langle e^{\int_0^T j(t)X_tdt}\right\rangle,
\qquad W[j]=\ln Z[j],\qquad Z[0]=1.
$$

在有限时间切片上，指数就是 $\sum_k\varepsilon j(t_k)x_k$；对某个时刻的源求导，会从指数中取出相应的位置因子。连续极限中，

$$
\left.\frac{\delta W}{\delta j(t)}\right|_0=\langle X_t\rangle,
\qquad
\left.\frac{\delta^2W}{\delta j(t)\delta j(s)}\right|_0
=\langle X_tX_s\rangle_c.
$$

这里 $\langle X_tX_s\rangle_c=\langle X_tX_s\rangle-\langle X_t\rangle\langle X_s\rangle$。对 $Z$ 求导生成原始矩，对 $\ln Z$ 求导会扣除均值等可分解部分，得到连通关联或更高阶累积量。生成泛函因而把许多不同的多时刻问题组织成同一个对象的导数。

$j$ 在这里是提取统计的记账工具。它暂时按观测量给历史加权，求导后再设为零；并没有在原运动方程中施加一个真实外力。下一节的驱动 $h(t)$ 才会改变运动规律，两种源的作用需要区分。

### OU 的 Gaussian 生成泛函，以及一次具体计算

对 Gaussian 过程，任意有限组时刻的位置联合分布都是多元高斯。有限维高斯积分的配方给出“源乘均值，加上源的二次型”；取时间切片极限，就得到

$$
W[j]=\int_0^T j(t)m(t)dt
+\frac12\int_0^Tdt\int_0^Tds\,j(t)C(t,s)j(s),
$$

其中 $m(t)$ 为均值、$C$ 为连通协方差。这个结果同时保留所有时刻，而非只描述某个终点。

先选取稳态准备的零均值 OU：$dX=-\gamma Xdt+\sqrt{2D}\,dW$，$X_0$ 从方差 $D/\gamma$ 的稳态分布采样，并与未来噪声独立。于是

$$
\boxed{W[j]=\frac D{2\gamma}
\int_0^Tdt\int_0^Tds\,
j(t)e^{-\gamma|t-s|}j(s).}
$$

两次源导数立即恢复 $C(t,s)=(D/\gamma)e^{-\gamma|t-s|}$。由于 $W$ 只有二次项，三阶及以上连通累积量为零；Gaussian 多时刻统计已由均值与协方差完整确定。

现在令整段时间内 $j(t)=\lambda$，则 $Z[j]$ 就变成累计量 $A_T$ 的矩母函数。将上述核积分，得到 $A_T$ 为零均值高斯变量，方差为

$$
\boxed{\operatorname{Var}(A_T)
=\frac{2D}{\gamma^2}
\left[T-\frac{1-e^{-\gamma T}}\gamma\right].}
$$

当 $T\ll\gamma^{-1}$，位置尚未显著失去记忆，$A_T$ 近似为 $TX_0$，方差按 $(D/\gamma)T^2$ 增长；当 $T\gg\gamma^{-1}$，许多相隔较远的时间段近似失去关联，方差的主导项改为与 $T$ 成正比。这是沿整条历史累计的涨落，从同一个生成泛函即可算出。

初始准备会改变这个计算。若固定 $X_0=x_0$，则 $m(t)=x_0e^{-\gamma t}$，协方差为 $(D/\gamma)[e^{-\gamma|t-s|}-e^{-\gamma(t+s)}]$；把它们代入一般 Gaussian 公式即可。不能只因后续方程相同，就将固定初态与稳态准备的有限时间统计混用。

## C.2 MSRJD：动力学约束与响应场 {#msrjd}

### 保留动力学约束，为什么有助于计算响应

Onsager–Machlup 直接用位置历史表达权重。若还想系统地研究驱动变化，或把非线性漂移组织成可展开的项，先保留运动方程的约束通常更方便。约束中漂移只出现一次，外源怎样进入也更加直接。

以 $\dot x=b(x)+h(t)+\eta(t)$ 为例，$\langle\eta(t)\eta(s)\rangle=2D\delta(t-s)$。以下仍用前点 Itô 离散，并保持初始准备固定、终点自由；因此不会对未来终点另加条件。将本步实际位移记为 $\Delta x_k$，噪声位移 $\eta_k$ 的方差为 $2D\varepsilon$，方程约束为

$$
\delta\left(\Delta x_k-\varepsilon b(x_k)-\varepsilon h_k-\eta_k\right).
$$

这里的 $h$ 是外部驱动，时间步长始终记为 $\varepsilon$。给定一组噪声实现，这份 Delta 约束只允许符合动力学的下一步位置；再按高斯噪声权重平均，就恢复原过程的转移概率。

### 从 Delta 约束完成一次基本推导

令 $q_k=\Delta x_k-\varepsilon b(x_k)-\varepsilon h_k$。对约束使用 Fourier 表示，再积分噪声：

$$
\delta(q_k-\eta_k)
=\int\frac{d\hat x_k}{2\pi}e^{-i\hat x_k(q_k-\eta_k)},
\qquad
\left\langle e^{i\hat x_k\eta_k}\right\rangle
=e^{-D\varepsilon\hat x_k^2}.
$$

各步相乘后的指数因此为

$$
\sum_k\left[-i\hat x_k
(\Delta x_k-\varepsilon b(x_k)-\varepsilon h_k)
-D\varepsilon\hat x_k^2\right].
$$

取 $\widetilde x_k=i\hat x_k$，并转为连续记号，得到 Martin–Siggia–Rose–Janssen–de Dominicis（MSRJD）作用量

$$
\boxed{S_h[x,\widetilde x]
=\int_0^Tdt\left[
\widetilde x(\dot x-b(x)-h)-D\widetilde x^2\right].}
$$

于是可将生成泛函写成

$$
Z[j,h]=\int dx_0\,p_0(x_0)
\int_{x(0)=x_0}\mathcal Dx\,\mathcal D\widetilde x\;
e^{-S_h[x,\widetilde x]+\int_0^T j x\,dt},
$$

终点自由积分，测度按上述离散构造规定，归一化为 $Z[0,h]=1$。对前点加性噪声，方程约束到位置变量的 Jacobian 没有额外路径依赖；换用别的规则时仍需相应调整。

这份作用量的两个部分具有不同来源。$\widetilde x(\dot x-b-h)$ 实施动力学约束，$-D\widetilde x^2$ 来自高斯噪声的特征函数。$\widetilde x$ 是引入约束时出现的辅助变量，并不是新增的物理坐标或独立噪声。[Chow 与 Buice 的综述](https://link.springer.com/article/10.1186/s13408-015-0018-5)也从离散 Delta 约束建立了这类响应场表示。

### 与 Onsager–Machlup 怎样相接

由于 $\hat x$ 沿实轴积分，$\widetilde x=i\hat x$ 的积分轮廓沿虚轴。作用量中 $-D\widetilde x^2$ 的符号必须连同这个轮廓理解，不能将它作为实轴上可归一化的普通概率密度。

如果重新积分辅助 Fourier 变量，每一步恰好给出

$$
\int\frac{d\hat x_k}{2\pi}
e^{-D\varepsilon\hat x_k^2-i\hat x_kq_k}
=\frac1{\sqrt{4\pi D\varepsilon}}
e^{-q_k^2/(4D\varepsilon)}.
$$

这就是前面的高斯一步权重，包括归一化因子。因此 OM 与 MSRJD 是同一随机动力学的不同表示：积分响应场得到位置路径权重，保留响应场则能更直接地提取驱动作用。与中点 OM 比较时，还要同时转换离散规则和 Jacobian。

### 辅助场为什么称为响应场

外力在作用量中只增加一项

$$
S_h=S_0-\int_0^T\widetilde x(s)h(s)ds.
$$

对 $h(s)$ 求导，会在路径平均中取出 $\widetilde x(s)$。令 $W[j,h]=\ln Z[j,h]$，便得到

$$
\boxed{R(t,s)
=\left.\frac{\delta\langle X_t\rangle_h}{\delta h(s)}\right|_{h=0}
=\left.\frac{\delta^2W}{\delta j(t)\delta h(s)}\right|_0
=\langle x(t)\widetilde x(s)\rangle_c.}
$$

归一化 $Z[0,h]=1$ 还给出 $\langle\widetilde x(s)\rangle=0$，因而这里也可以省略连通下标。响应场的物理作用正由这个导数关系定义：改变真实驱动怎样影响位置统计，可以在扩展路径表示中用一次混合关联计算。

因果性要求早于扰动的观测不受其影响，$t\lt s$ 时 $R(t,s)=0$。这里使用的是没有未来终点约束的原过程；bridge 中对未来结果的筛选不能与这种无条件因果响应混淆。前点 Itô 规则还规定了等时求值，本篇按离散约定取 $\Theta(0)=0$。

### OU：在同一泛函中恢复关联与响应

取 $b=-\gamma x$，未扰动作用量为

$$
S_0=\int_0^Tdt\left[
\widetilde x(\partial_t+\gamma)x-D\widetilde x^2\right].
$$

二次作用量中的确定性算子为 $\partial_t+\gamma$。其迟致 Green 函数满足 $(\partial_t+\gamma)R(t,s)=\delta(t-s)$，并在 $t\lt s$ 为零，所以

$$
\boxed{R(t,s)=\Theta(t-s)e^{-\gamma(t-s)}.}
$$

扰动产生的位移会按与普通弛豫相同的速率消失。对固定初始准备的小驱动，平均改变为 $\delta m(t)=\int_0^tR(t,s)h(s)ds$，与第一篇直接解运动方程的结果一致。

对前面选定的稳态初始分布，协方差包含初始涨落的残余，以及后来噪声的共同贡献：

$$
C(t,s)=\frac D\gamma e^{-\gamma(t+s)}
+2D\int_0^{\min(t,s)}R(t,u)R(s,u)du
=\frac D\gamma e^{-\gamma|t-s|}.
$$

一次早于两个观测时刻的噪声，会分别通过两个响应核影响它们；对共同来源累加，就得到两时关联。晚于其中一个观测时刻的噪声不能改变已经发生的位置，因此积分只到两时刻的较小值。关联读取共同涨落，响应读取一次受控扰动的后果，二者在同一个 Gaussian 理论中有清楚的不同作用。

![OU 的时间关联与因果响应核](../figures/08_ou_correlation_response.svg)

*图 C1｜$\gamma=1,D=0.4$，横轴为时刻 $s$、纵轴为时刻 $t$。左图的稳态关联关于 $t=s$ 对称，两侧的时刻都可以共享过去噪声；右图的响应只在 $t>s$ 一侧有权重。色标分别表示各自的数值，关联与响应的量纲和幅度不相同。*

外部 $h$ 对 OU 只改变均值，不改变未扰动的协方差。因此其完整 Gaussian 生成泛函还可写为

$$
W[j,h]=\int_0^Tdt\,j(t)\int_0^t ds\,R(t,s)h(s)
+\frac12\int_0^Tdt\int_0^Tds\,j(t)C(t,s)j(s).
$$

两次 $j$ 导数读取 $C$，一次 $j$ 和一次 $h$ 导数读取 $R$。这将前面的路径统计与真实动力学响应放进了同一个计算结构，而不是在得到一个形式作用量后便结束。

### 非线性项与计算的扩展

若漂移改成 $b(x)=-\gamma x-gx^3$，取 $g>0$，作用量在 Gaussian 部分之外增加 $g\int\widetilde x x^3dt$。这类非线性项会耦合不同阶的涨落，在微扰组织中称为相互作用顶点。生成泛函仍然定义关联与响应，Gaussian 理论提供展开的起点；这里不继续发展图形计算或重整化。

## C.3 中点 OM 表示中的 Jacobian {#om-jacobian}

前点离散中，定义随机位移 $q_k=x_{k+1}-x_k-b(x_k)\varepsilon$。固定 $x_0$ 后，从 $x_1,\ldots,x_N$ 到 $q_0,\ldots,q_{N-1}$ 的 Jacobian 矩阵是三角的，对角元均为一，因此没有额外的路径依赖因子。

若采用中点规则，漂移改为 $b(\bar x_k)$，$\bar x_k=(x_k+x_{k+1})/2$。这时改变下一步位置，也同时改变了本步用来评估漂移的位置。小步局部换元给出

$$
\left|\frac{\partial q_k}{\partial x_{k+1}}\right|
=\left|1-\frac\varepsilon2b'(\bar x_k)\right|,
\qquad
\prod_k\left|1-\frac\varepsilon2b'(\bar x_k)\right|
\simeq e^{-\frac12\sum_k\varepsilon b'(\bar x_k)}.
$$

因而常用的中点 Onsager–Machlup 表示写为

$$
\boxed{S_{\mathrm{OM}}[x]
=\int_0^Tdt\left[
\frac{(\dot x-b(x))^2}{4D}+\frac12b'(x)
\right],}
$$

并配合该规则的测度。Jacobian 补偿了变量换元时邻域体积的伸缩；它没有增加一种新的物理力。对这里的加性噪声，前点与中点方案可以收敛到同一扩散，但各自的漂移求值、测度和 Jacobian 必须一致使用。

对于光滑曲线附近的小管概率，Onsager–Machlup 泛函也有相应解释，仍需指定邻域与参考测度。这里不展开一般几何随机微积分，只保留关键区别：**平方噪声代价给出了核心权重结构，完整的有限噪声路径表示还包含换元和测度信息。** 关于这些构造的进一步讨论，可参见 [Chow 与 Buice 的综述](https://link.springer.com/article/10.1186/s13408-015-0018-5)。

**返回：**[路径积分正文](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/path_integrals/) · [专题首页](/notes/nonequilibrium_statistical_physics/stochastic_dynamics/)。
