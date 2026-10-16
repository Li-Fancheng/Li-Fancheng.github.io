# 基因调控网络如何编码空间位置
## 从分子计数噪声到位置解码与受约束的信息优化

> 西湖大学 2026 年秋季《定量生物学导论》课程学习笔记。

> **核心问题：** 在分子数有限、表达有噪声、蛋白质能够扩散的条件下，胚胎中的细胞核能从局部基因表达中多精确地识别自己的位置？

## 概述 {#main-lecture}

果蝇胚胎沿前后轴形成一系列基因表达区域，后续的分节过程再读取这些表达。一个位于胚胎中部的细胞核，不能直接测量自己距离前端多少微米，却能接触局部的调控蛋白。于是，空间位置通过分子的浓度和组合，转化为可由下游网络读取的信号。

先看一个具体的输入信号：**Bicoid（Bcd）是果蝇早期胚胎中的母源转录因子，其蛋白质浓度从前端向后端降低。** 母源 *bcd* mRNA 定位于前端，Bcd 蛋白参与调控 *hunchback*（*hb*）等下游基因。需要区分 mRNA 定位、蛋白质浓度梯度和下游转录图案三个层次；[Porcher 和 Dostatni（2010）的综述](https://doi.org/10.1016/j.cub.2010.01.026)可作为认识这个系统的背景读物。

<div id="bicoid-system"></div>

![Bicoid 前后轴浓度梯度与早期 hunchback 前部转录区域的文献示意图。](figures/lucas_2018_fig01.png)

*Fig. 1. 从连续梯度到空间表达边界。* 引自 [Lucas 等（2018），原文 Fig. 1](https://doi.org/10.1371/journal.pgen.1007676.g001)，完整原图，未改动，[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)。A：红色表示 Bicoid 浓度，A、P 分别指前端和后端。B：根据 RNA FISH 观测绘制的早期 *hb* 转录示意，绿色核正在转录，白色核未检测到转录；下方曲线表示活跃核的比例。两幅对应核分裂周期 nc8 与 nc11，图中的 30 min 是这两个阶段的间隔。

这张图帮助提出问题：平缓变化的输入，怎样形成较陡的表达边界？图中的细胞核排列是示意，真实胚胎的荧光测量见第 1 节的 [Fig. 2](#experimental-bicoid)。这一早期 *hb* 响应只是网络的一部分；后文多基因的位置编码还要考虑其他输入与基因间调控。此时胚胎仍是合胞体，多个细胞核共享细胞质。

这个过程首先是一个物理测量问题。同一位置的分子数会波动；同一种表达水平也可能出现在不同位置。要判断一个表达图案能否可靠指定位置，既要知道平均表达怎样随位置变化，也要知道噪声多大，以及不同位置的表达分布是否重叠。

[Petkova 等（2019）](https://pmc.ncbi.nlm.nih.gov/articles/PMC6526179/)用四种 gap gene 的联合表达构造位置后验分布，发现其定位精度可达约胚胎长度的 1%。更早的[位置信息测量](https://doi.org/10.1073/pnas.1315642110)已经表明，把这些表达压成开关变量，会丢失连续表达水平所携带的信息。这些结果提出了两个连续的问题：**如何从有噪声的表达计算位置精度？在给定分子资源下，什么样的网络能够产生这样的编码？**

正文依次建立梯度、局部估计、全局解码和信息优化四个层次。第 1 节说明梯度产生与可靠定位的区别；第 2、3 节把分子噪声换算成位置误差；第 4 节处理局部估计遗漏的空间歧义；第 5 节推导一个可手算的最优表达图案；第 6 节讨论扩散对信号和噪声的双重影响；第 7、8 节连接实际网络优化，并检验平均信息量能说明什么。

较长的推导放在[技术附录](#technical-appendices)。Fig. 1 用文献示意图说明系统，Fig. 2 给出真实胚胎的测量；Fig. 3–8 是明确假设下的模型计算，用于逐步解释定位机制。

### 记号与空间区间 {#notation}

用 $L$ 表示整个胚胎的前后轴长度，$x$ 表示绝对位置，$z=x/L$ 表示相对位置。统计分析可以只覆盖子区间 $\Omega=[x_a,x_b]$，其长度为 $W=x_b-x_a$；若讨论整个轴，才有 $W=L$。除特别说明外，位置在所分析的区间内均匀取样，即 $P_X(x)=1/W$。

| 记号 | 含义 |
|---|---|
| $c(x,t)$、$\mathbf c(x,t)$ | 单个或多个母源形态发生素输入 |
| $\mathbf G$、$\mathbf g$ | 随机表达向量及其一次观测值 |
| $\boldsymbol\mu(x)$ | 固定位置的平均表达 |
| $\Sigma(x)$ | 同一位置的表达噪声协方差 |
| $\sigma_x(x)$ | 局部位置不确定性的标准差，使用时需注明估计条件 |
| $\mathcal J(x)$ | 单次表达观测对位置的 Fisher information |
| $I(X;\mathbf G)$ | 在指定位置先验下，表达携带的互信息，单位为 bits |
| $\theta$、$\mathcal C$ | 可调网络参数及其物理可行集合 |

---

## 1. 梯度如何产生，又为什么不足以保证精确定位 {#gradient}

### 1.1 源、扩散和降解 {#source-diffusion-degradation}

先把梯度对应到可测量的量。[Gregor 等（2007）](https://pmc.ncbi.nlm.nih.gov/articles/PMC2253670/)通过 Bcd 与 GFP 的融合蛋白，观察活胚胎并标定逐核浓度。

<div id="experimental-bicoid"></div>

![活果蝇胚胎中 Bcd-GFP 的双光子显微图，以及逐个细胞核测得的前后轴 Bicoid 浓度。](figures/gregor_2007_fig02.jpg)

*Fig. 2. 从果蝇胚胎到 Bicoid 浓度梯度。* 引自 [Gregor 等（2007），原文 Fig. 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC2253670/#F2)，完整原图，未改动。A：活胚胎中的 Bcd-GFP 荧光，前端在左，亮点对应富集 Bcd-GFP 的细胞核，标尺 $50\,\mu\mathrm m$；胚胎外的 GFP 溶液用于浓度标定。B：逐核浓度随相对位置 $x/L$ 的变化，红、蓝点分别来自背侧和腹侧，黑点及误差棒是无 GFP 胚胎的背景对照；插图展示扣除平均背景后的中后段信号。

读图时把 A、B 对应起来：前端较亮，逐核浓度也较高；向后亮度减弱，浓度随之下降。单个胚胎中各核之间的散布同时受到测量误差、空间差异等影响，不能直接当成同一位置跨胚胎的全部噪声。后文将用平均曲线描述空间趋势，用明确取样条件下的分布描述表达变异。

**下面建立一个源–扩散–降解模型（SDD）。** 仅凭近似指数的浓度曲线，尚不能唯一确定梯度形成机制。[Durrieu 等（2018）](https://doi.org/10.15252/msb.20188355)比较了多种能产生相近浓度图案的机制，并用荧光计时器测量蛋白质年龄，以额外观测区分它们。这里先用 SDD 分离源强、扩散和寿命的作用。

先考虑一维区间 $x\in[0,L]$。分子从前端注入，在空间中扩散，并以寿命 $\tau$ 降解：

$$
\partial_t c=D\partial_x^2c-\frac{c}{\tau},
\qquad -D\partial_xc(0,t)=R,
\qquad \partial_xc(L,t)=0.
\tag{1.1}
$$

这里 $R$ 是边界注入通量；其单位由一维浓度的定义决定。后端采用零通量边界。稳态下，扩散输运与降解平衡，特征长度为 $\lambda=\sqrt{D\tau}$。满足两个边界条件的解是

$$
\boxed{
c(x)=\frac{R\tau}{\lambda}
\frac{\cosh[(L-x)/\lambda]}{\sinh(L/\lambda)}.
}
\tag{1.2}
$$

当 $L/\lambda\gg1$，且观察位置距后端有数个 $\lambda$，反射边界的影响很小，得到常用的指数近似：

$$
c(x)\simeq c_0e^{-x/\lambda},\qquad
c_0=\frac{R\tau}{\lambda}.
\tag{1.3}
$$

因此，局部浓度可以随到前端的距离单调变化。若下游在 $c=c_*$ 处切换表达，就得到一条空间边界。但这只说明边界如何出现，还没有说明它在噪声下能有多稳定。

### 1.2 边界位置怎样响应源强变化 {#threshold-variability}

在指数近似内，阈值位置满足

$$
c_0e^{-x_*/\lambda}=c_*
\quad\Longrightarrow\quad
x_*=\lambda\ln\frac{c_0}{c_*}.
\tag{1.4}
$$

若只有 $R$ 改变，$D$、$\tau$ 和 $c_*$ 固定，则小扰动产生

$$
\delta x_*\simeq\lambda\frac{\delta R}{R}.
\tag{1.5}
$$

例如，假设 $\lambda/L=0.2$，源强的相对变化为 10%，则边界移动约 $0.02L$。这是用于说明敏感性的假设算例；不能据此断言实验中的源强波动必然就是 10%。有限变化的精确结果为 $\Delta x_*=\lambda\ln(1.1)\simeq0.0191L$。

如果 $\lambda$、$c_0$ 和阈值也改变，一阶变化一般写成

$$
\delta x_*=
\ln\frac{c_0}{c_*}\,\delta\lambda
+\lambda\left(\delta\ln c_0-\delta\ln c_*\right).
\tag{1.6}
$$

这些变化之间可以相关。因而，“浓度梯度存在”与“边界对扰动稳定”是不同的结论，后者需要计算具体参数怎样一起变化。

![指数梯度、源强变化引起的阈值移动、后端边界效应，以及 Poisson 分子计数下随位置变化的误差。](figures/01_gradient_and_counting.svg)

*Fig. 3. 梯度的形成与定位限制。* A：$\lambda/L=0.2$，源强增加 10% 时，同一阈值对应的边界向后移动。B：式（1.2）与指数近似；后端的零通量条件使曲线变平，指数近似在该处失效。C：第 2 节的独立 Poisson 计数模型，$m(x)=m_0e^{-x/\lambda}$；图中使用半无限指数近似，局部位置误差随分子数降低而增大。

### 1.3 精度、鲁棒性和尺度缩放 {#precision-robustness-scaling}

还需要区分三种测量。**精度**关心给定条件下重复观测的离散程度；**鲁棒性**关心源强、降解率等参数改变时图案怎样移动；**尺度缩放**关心改变胚胎长度后，边界能否保留相同的相对位置 $z_*=x_*/L$。

在指数模型里，$z_*=(\lambda/L)\ln(c_0/c_*)$。要让 $z_*$ 不随 $L$ 改变，需要 $\lambda/L$ 与浓度比的变化共同满足相应条件；仅仅把扩散方程放到一个更长的区间，并不会自动得到缩放。另一方面，即使不同大小胚胎的平均边界都对齐，单个胚胎中的读出仍可能有很大噪声。

下面先固定输入的平均图案，考察分子计数的随机性怎样限制精度。

## 2. 从表达噪声到位置误差 {#local-estimation}

### 2.1 测量对象是条件分布 {#conditional-expression}

同一个位置的表达量并不唯一。转录、翻译、降解和分子到达都具有随机性，所以编码的基本对象是 $P(g\mid x)$。平均曲线 $\mu(x)=\mathbb E[G\mid x]$ 描述信号，方差 $\sigma_g^2(x)=\operatorname{Var}(G\mid x)$ 描述固定位置上的表达波动。

先在真实位置 $x_0$ 附近写成

$$
g=\mu(x_0)+\eta,\qquad
\mathbb E[\eta\mid x_0]=0,
\qquad \operatorname{Var}(\eta\mid x_0)=\sigma_g^2(x_0).
$$

若平均曲线在此可逆，且噪声足够小，反推位置得到 $\mu(\hat x)=g$。将左侧在 $x_0$ 展开：

$$
\mu(x_0)+\mu'(x_0)(\hat x-x_0)
\simeq\mu(x_0)+\eta.
$$

所以

$$
\boxed{
\sigma_x^2(x_0)\simeq
\frac{\sigma_g^2(x_0)}{[\mu'(x_0)]^2}.
}
\tag{2.1}
$$

分子噪声沿着平均表达曲线换算成了位置误差。同样大小的表达波动，遇到更大的空间斜率，对应较小的位置变化；平台区的斜率接近零，则难以依靠均值区分相邻位置。

式（2.1）是局部线性误差传播公式。若曲线有多个分支，必须先知道观测来自哪一支；若 $\mu'=0$，需要重新分析高阶变化，不能直接把发散的近似当作完整后验的标准差。

### 2.2 Poisson 计数给出一个可计算的极限 {#poisson-counting}

设细胞核从固定采样体积得到一次独立分子计数 $N$，其条件分布为

$$
P(n\mid x)=\frac{m(x)^ne^{-m(x)}}{n!},
\qquad m(x)=m_0e^{-x/\lambda}.
\tag{2.2}
$$

Poisson 分布满足 $\operatorname{Var}(N\mid x)=m(x)$，而 $m'=-m/\lambda$。代入局部误差公式：

$$
\boxed{\sigma_x(x)\simeq\frac{\lambda}{\sqrt{m(x)}}.}
\tag{2.3}
$$

这里有两个同时存在的限制：有限分子数使相对计数噪声为 $m^{-1/2}$；有限梯度长度又把相对浓度误差乘上 $\lambda$，转化成位置误差。

假设 $\lambda/L=0.2$，要求某个位置有 $\sigma_x/L\lesssim0.01$，则该处需要 $m(x)\gtrsim400$。这指局部有效独立计数的数量，不是整个胚胎的蛋白质总数，也不是所有调控机制共有的分子数下界。沿指数梯度向后，$m(x)$ 降低，误差按 $e^{x/(2\lambda)}$ 增长。

Poisson 模型还可以不经高斯近似直接计算 Fisher information，得到 $\mathcal J=(m')^2/m$；因此式（2.3）也对应局部估计界。其可达到性仍需要足够的有效计数，见[附录 B](#appendix-b)。

### 2.3 多观察几次，是否总能消除噪声 {#averaging-and-extrinsic-noise}

如果每次计数独立，平均 $K$ 次就使方差降低为原来的 $1/K$，位置误差相应降低为 $1/\sqrt K$。实际中，连续采样会反复遇到同一批分子，启动子也可能在一个状态停留很久；采样次数与独立样本数因而不同。

对于相关时间为 $\tau_c$ 的平稳过程，当观察时长 $T\gg\tau_c$，时间平均的方差通常按 $2\tau_c/T$ 缩减，而非按相机帧数缩减。这个因子由完整时间相关函数的积分决定，推导见[附录 D](#appendix-d)。有限的发育时间会限制能够利用的独立观测。

胚胎间的共同波动还有另一层问题。假设源强 $R$ 未知，$\ln m=\ln m_0(R)-x/\lambda$；一次局部计数只能约束这两个变量的组合。观测到较低浓度，既可能是位置更靠后，也可能是源强较弱。若没有关于 $R$ 的先验、额外输入或其他观测，单个读数不能同时确定位置与源强。增加分子数可以减少计数误差，却不能消除这种参数混淆。

## 3. 多个基因如何联合提高位置分辨能力 {#multigene-fisher}

### 3.1 表达空间中的信号与噪声 {#expression-space}

四个 gap genes 的表达构成向量，例如

$$
\mathbf g=(g_{\mathrm{Hb}},g_{\mathrm{Kr}},g_{\mathrm{Gt}},g_{\mathrm{Kni}}).
$$

当 $x$ 改变时，平均表达 $\boldsymbol\mu(x)$ 在四维表达空间里形成一条曲线。同一位置上的重复观测围绕曲线散开；协方差 $\Sigma(x)$ 描述散开的尺度和方向。空间上相隔 $\delta x$ 的两个位置，其平均表达之差为 $\boldsymbol\mu'\delta x$。

先假设局部噪声为多元高斯，且 $\Sigma$ 在所考察的位置范围内近似不变。高斯似然中，表达差异以 $\Sigma^{-1}$ 加权，因此可区分性的平方随位置间距变化为

$$
(\delta\boldsymbol\mu)^{\mathsf T}\Sigma^{-1}
(\delta\boldsymbol\mu)
\simeq
\underbrace{\boldsymbol\mu'^{\mathsf T}\Sigma^{-1}\boldsymbol\mu'}_{\mathcal J_\mu(x)}
\,\delta x^2.
\tag{3.1}
$$

噪声大的方向对表达差异的权重较低，噪声小的方向权重较高。所以，只列出每个基因的方差还不够，还需要知道噪声椭球相对于平均曲线切向的朝向。

### 3.2 Fisher information 与局部误差界 {#fisher-bound}

对一般条件分布，定义

$$
\mathcal J(x)=
\mathbb E_{\mathbf G\mid x}
\left[\left(\partial_x\ln P(\mathbf G\mid x)\right)^2\right].
\tag{3.2}
$$

在适当的正则条件下，对无偏位置估计量，Cramér–Rao 不等式给出

$$
\operatorname{Var}(\hat X\mid x)\geq\frac1{\mathcal J(x)}.
\tag{3.3}
$$

在小噪声、单一可识别分支及缓慢变化先验下，局部后验宽度也可近似为 $\sigma_x\simeq\mathcal J^{-1/2}$。但这是一个可达到的近似尺度，不表示任意实际解码器都达到此界；对有偏估计量，式（3.3）也不能原样套用。

若高斯噪声的协方差固定，则 $\mathcal J=\mathcal J_\mu$。若它依赖位置，完整结果为

$$
\boxed{
\mathcal J=
\boldsymbol\mu'^{\mathsf T}\Sigma^{-1}\boldsymbol\mu'
+\frac12\operatorname{Tr}
\left(\Sigma^{-1}\Sigma'\Sigma^{-1}\Sigma'\right).
}
\tag{3.4}
$$

第二项说明，噪声分布的形状本身也可能随位置变化并提供信息。实际使用“斜率除以噪声”时，通常只保留第一项；应明确这一近似，不能把它当作所有分布下 Fisher information 的定义。

### 3.3 不同基因可以在空间上互相补充 {#complementary-slopes}

如果各基因的噪声在给定位置条件下独立，且忽略方差变化项，则

$$
\mathcal J_\mu(x)=
\sum_{\alpha=1}^{K}
\frac{[\mu_\alpha'(x)]^2}{\sigma_\alpha^2(x)}.
\tag{3.5}
$$

一个基因在表达平台上几乎没有斜率，另一个基因可能正处于表达边界。它们的贡献可以互相补充，使沿整个轴都有表达随位置明显变化。新增测量不会降低已有测量包含的信息，因为解码器始终可以忽略它；但新增基因若分走固定的分子预算，原有通道的噪声也会改变，收益就不能只靠式（3.5）逐项相加来判断。

### 3.4 相关噪声究竟有利还是有害 {#correlation-geometry}

考虑两个基因，局部噪声为

$$
\Sigma=\sigma^2
\begin{pmatrix}1&\rho\\\rho&1\end{pmatrix},
\qquad |\rho|<1.
$$

如果两个均值同向变化，$\boldsymbol\mu'=(a,a)$，则

$$
\mathcal J_\mu=\frac{2a^2}{\sigma^2(1+\rho)}.
\tag{3.6}
$$

共同升降的噪声与位置信号方向一致，正相关会降低可分辨性。但若两个基因反向变化，$\boldsymbol\mu'=(a,-a)$，结果为

$$
\mathcal J_\mu=\frac{2a^2}{\sigma^2(1-\rho)}.
\tag{3.7}
$$

此时正相关噪声可以在差值中抵消，反而有助于读取位置。因而，相关系数的正负本身不足以判断编码好坏，必须同时看信号方向。$|\rho|\to1$ 的发散对应理想模型里某个方向噪声趋零；真实系统中的额外噪声会限制这一增益。

![两个基因同向和反向变化时，相关噪声椭球相对于信号方向的几何关系，以及对应 Fisher information。](figures/02_noise_geometry.svg)

*Fig. 4. 噪声相关性的效果取决于信号方向。* A、B：蓝色椭圆是 $\rho=0.8$ 时的一标准差等高线，橙色箭头表示位置改变引起的均值移动。C：固定单基因斜率大小与噪声方差，分别计算式（3.6）和（3.7）。同一相关结构可以使一种位置编码恶化，同时改善另一种。

## 4. 从局部精度走向全局位置解码 {#global-decoding}

### 4.1 陡峭的表达曲线也可能有位置歧义 {#branch-ambiguity}

考虑以胚胎中部为中心的对称表达峰。左右两个位置 $x_1$ 和 $x_2$ 可以具有相同均值；若它们的条件分布也相同，单独观测这个基因就无法区分二者。每个分支上的斜率可以很大，因此局部 Fisher information 也可以很大，但这种局部可分辨性没有消除远处的另一种解释。

Bayes 公式把所有候选位置放在一起比较：

$$
\boxed{
P(x\mid\mathbf g)=
\frac{P(\mathbf g\mid x)P_X(x)}
{\int_\Omega P(\mathbf g\mid u)P_X(u)\,du}.
}
\tag{4.1}
$$

在对称峰的例子里，后验可能有两个分离的峰。此时每个峰的局部宽度很小，仍不能说明全局定位准确。若增加一个能区分左右分支的基因，表达向量便可能只对应一个位置。

![对称表达峰导致双分支解码，增加一个单调表达通道后，平均后验集中在真实位置附近。](figures/03_global_ambiguity.svg)

*Fig. 5. 局部精度与全局歧义。* A：两个示意表达通道，$\mu_1=\exp[-(z-0.5)^2/(2\times0.19^2)]$，$\mu_2=z$。B、C：分别用第一个通道和两个通道解码，噪声独立且 $\sigma_g=0.06$。每一列是在真实位置 $x$ 对表达噪声平均后的后验密度，定义见附录 C。单基因图中出现两条分支；第二个通道帮助识别左右位置。这些曲线是教学模型，不是 gap gene 的拟合。

这也解释了实验中的解码图为何有用：[Petkova 等](https://pmc.ncbi.nlm.nih.gov/articles/PMC6526179/#F3)比较不同数量基因的联合读出，用完整的位置后验同时展示局部宽度和远距离歧义。只报告一个位置误差数字，会丢掉其中一部分结构。

### 4.2 互信息把所有表达观测纳入同一个量 {#mutual-information}

在指定位置先验下，表达的边缘分布是

$$
P(\mathbf g)=\int_\Omega P_X(x)P(\mathbf g\mid x)\,dx.
$$

互信息定义为

$$
\boxed{
I(X;\mathbf G)=
\int_\Omega dx\,P_X(x)
\int d\mathbf g\,P(\mathbf g\mid x)
\log_2\frac{P(\mathbf g\mid x)}{P(\mathbf g)}.
}
\tag{4.2}
$$

它比较知道表达前后，位置的不确定性平均减少了多少。若使用连续位置，可以写成 $I=h(X)-h(X\mid\mathbf G)$。其中两个微分熵分别依赖坐标单位，但它们的差不依赖单位；互信息因此是有明确意义的 bits 数。

若任意两个位置都给出相同的条件分布，$P(\mathbf g\mid x)=P(\mathbf g)$，互信息就是零。反之，不同位置的表达分布越容易区分，知道表达后的位置后验通常越集中。互信息还会计入多个分离分支造成的剩余不确定性，局部斜率公式则不会自动计入。

### 4.3 什么时候能把 bits 换算成位置误差 {#bits-and-error}

假设位置在长度为 $W$ 的区间上均匀分布；典型位置后验只有一个窄峰，并且可以近似为高斯；边界影响与宽度随观测变化的高阶修正均可忽略。此时

$$
\begin{aligned}
I\simeq{}&\log_2 W\\
&-\left\langle
\log_2\left[\sqrt{2\pi e}\,\sigma_x(x)\right]
\right\rangle_x.
\end{aligned}
\tag{4.3}
$$

定义几何平均误差 $\sigma_{x,\mathrm{geo}}=\exp\langle\ln\sigma_x\rangle_x$，则

$$
\boxed{
I\simeq\log_2\frac{W}{\sqrt{2\pi e}\,\sigma_{x,\mathrm{geo}}}.
}
\tag{4.4}
$$

如果误差在空间上近似恒定，几何平均才可以直接替换为共同的 $\sigma_x$。因此，bits 所对应的“误差”既不是无条件的算术平均，也不是均方根，更不是最坏位置的误差。

| $I$ | 以分析区间计的 $\sigma_{x,\mathrm{geo}}/W$ | 若 $W=0.8L$，相应的 $\sigma_{x,\mathrm{geo}}/L$ |
|---|---|---|
| 3 bits | 3.02% | 2.42% |
| 4 bits | 1.51% | 1.21% |
| 4.2 bits | 1.32% | 1.05% |
| 4.3 bits | 1.23% | 0.98% |
| 5 bits | 0.76% | 0.60% |

表中仅计算式（4.4），没有加入实验拟合。它说明，空间区间的选择足以影响“4.2 bits 对应多少百分比”的数值。[Tkačik 等（2015）](https://pmc.ncbi.nlm.nih.gov/articles/PMC4286692/)的分析通常取前后轴中间 80% 的区域；比较文献结果时，必须保留各自的位置范围、取样方式和发育时刻。

### 4.4 bits 不是无误差分类的区域数 {#discrete-and-continuous}

若位置变量是 $N$ 个等概率细胞核的离散标签，互信息有严格上限 $I\leq H(X)=\log_2N$。达到上限意味着观测表达后，核标签的不确定性为零。若位置作为连续变量处理，则没有一个仅由 $L$ 决定的有限熵上限；此时有限信息来自观测噪声和其他约束。

因此，$2^{4.2}\simeq18.4$ 可以帮助理解平均不确定性的缩减尺度，却不能解读成“只能区分 18 个固定区域”，更不能据此给出相邻细胞核的分类正确率。分类误差还取决于核间距、后验形状、损失函数和解码规则。

另一个严格结论来自输出取值数。若每个基因只有一个二值读数 $B_i\in\{0,1\}$，四个基因至多有 16 种组合，因而

$$
I(X;B_1,B_2,B_3,B_4)\leq H(\mathbf B)\leq4\ \mathrm{bits}.
\tag{4.5}
$$

只要可靠估计的连续表达信息超过 4 bits，就表明这种单时刻、四二值通道的描述丢失了信息；并非说明四个基因在物理上只能传递四个 bits。

## 5. 信息最大化如何预测表达曲线 {#optimal-encoding}

### 5.1 固定噪声和表达范围的单基因问题 {#linear-optimum}

先构造一个可以完全求解的问题。位置在 $[0,W]$ 上均匀分布，单基因读数满足

$$
G=\mu(x)+\eta,\qquad
\eta\sim\mathcal N(0,\sigma_g^2).
$$

要求 $\mu$ 严格单调递增，端点固定为 $g_{\min}$ 和 $g_{\max}$，表达范围 $\Delta g=g_{\max}-g_{\min}$ 固定。噪声小到可以使用局部单峰近似。代入 $\sigma_x\simeq\sigma_g/\mu'$，互信息为

$$
I\simeq\frac1W\int_0^W
\log_2\frac{W\mu'(x)}{\sqrt{2\pi e}\,\sigma_g}\,dx.
\tag{5.1}
$$

约束给出 $\int_0^W\mu'(x)\,dx=\Delta g$。由于 $\log$ 是严格凹函数，Jensen 不等式意味着

$$
\frac1W\int_0^W\log_2\mu'(x)\,dx
\leq\log_2\frac{\Delta g}{W}.
$$

等号要求所有位置的斜率相同，因此

$$
\boxed{
\mu^*(x)=g_{\min}+\frac{\Delta g}{W}x,
\qquad
I_{\mathrm{opt}}^{\mathrm{small\ noise}}
\simeq\log_2\frac{\Delta g}{\sqrt{2\pi e}\,\sigma_g}.
}
\tag{5.2}
$$

有限表达范围在所有位置之间均匀分配，给出线性梯度和均匀的局部位置精度。陡峭阶跃把很大一部分表达变化集中在很窄区域，虽然阈值附近定位更好，两侧平台却几乎无法区分位置；平均信息不一定增加。

这个解只是在上述单调、小噪声模型内的最优解。有限噪声下，允许更一般输入分布的幅度受限高斯通道，其精确容量问题并不等同于这里的变分问题。若式（5.2）算出负值，则小噪声近似已经失效；真实互信息始终非负。

![在相同表达动态范围和噪声下，线性梯度与陡峭 sigmoid 的编码表现，以及小噪声互信息近似的适用范围。](figures/04_optimal_gradient.svg)

*Fig. 6. 阶跃更陡，并不保证平均编码更好。* A：两个单调图案具有相同的 $[0,1]$ 动态范围与独立高斯噪声 $\sigma_g=0.025$；图例中的互信息由 $I=h(G)-h(\eta)$ 数值积分计算。B：局部误差传播公式显示 sigmoid 在中央更精确，在两侧明显变差；较大误差处局部线性近似不再可靠。C：线性图案的数值互信息与式（5.2）比较，噪声增大后可见有限区间边界修正。

### 5.2 噪声依赖表达时，应该均分什么 {#noise-weighted-range}

假设噪声标准差是表达水平的函数 $\sigma_g(\mu)>0$，而不是常数。仍在单调、小噪声条件下，引入新的表达坐标

$$
q(g)=\int_{g_{\min}}^g\frac{du}{\sigma_g(u)}.
\tag{5.3}
$$

这个坐标以当地的噪声为单位衡量表达变化。因为 $dq/dx=\mu'/\sigma_g(\mu)$，互信息中的可调项变为 $\langle\log(dq/dx)\rangle_x$。固定端点意味着总的 $q$ 范围固定，于是同样的凹性论证给出

$$
\boxed{
\frac{\mu'(x)}{\sigma_g(\mu(x))}
=\frac1W\int_{g_{\min}}^{g_{\max}}
\frac{dg}{\sigma_g(g)}.
}
\tag{5.4}
$$

最优图案均匀分配的是以噪声为单位的表达变化。若 $\sigma_g(g)\propto\sqrt g$，则 $q\propto\sqrt g$，最优的 $\sqrt{\mu(x)}$ 随位置线性变化，$\mu(x)$ 为相应的二次曲线。这个推导要求表达范围内的计数足够大，不能把连续高斯近似外推到零分子数。

### 5.3 多基因优化还需要动力学约束 {#why-networks}

多基因情形中，$\boldsymbol\mu(x)$ 是表达空间里的曲线，$\Sigma$ 定义了沿不同方向的噪声尺度。忽略协方差变化项，式（4.3）中的可调部分近似为

$$
\frac12\int_\Omega P_X(x)\log_2\mathcal J_\mu(x)\,dx.
\tag{5.5}
$$

这一表达鼓励网络在整个区间保留可分辨的变化。但它没有要求曲线在远处保持分离：不同位置可以在表达空间中重新靠近，甚至相交。所以，单独最大化局部 Fisher information 不能替代完整互信息。

此外，表达曲线不能任意画。每个位置遵循同一套调控规则，表达有合成速率和降解时间，蛋白质可以交换，调控强度也有上限。某条几何上很好的曲线，未必能在给定时间内由可行的基因网络产生。只有把这些限制放入可行集合，最优图案才具有机制上的含义。

## 6. 扩散为什么既能改善精度，也会抹平图案 {#diffusion}

### 6.1 平均场中的空间平滑 {#diffusion-signal}

在合胞体胚胎中，细胞核尚未由完整细胞膜逐一隔开，蛋白质能够在相邻核区域之间交换。考虑一个用于分析的线性方程：

$$
\partial_t g=s(x)-\gamma g+D\partial_x^2g+\xi(x,t).
\tag{6.1}
$$

先忽略噪声。若产生率具有 Fourier 模式 $s_k$，稳态平均表达为

$$
\bar g_k=\frac{s_k}{\gamma+Dk^2}.
\tag{6.2}
$$

较大的 $D$ 对高波数、短尺度的结构抑制更强。相应的平滑长度为 $\ell_D=\sqrt{D/\gamma}$。这可以使边界变宽、降低局部斜率，从而削弱位置信号。

### 6.2 噪声是否一起减小，取决于噪声来源 {#diffusion-noise}

若 $\xi$ 代表各处独立产生的额外表达噪声，扩散会把局部偏高、偏低的表达混合。第 $k$ 个线性涨落模态满足

$$
\partial_t\delta g_k=-(\gamma+Dk^2)\delta g_k+\xi_k.
$$

当 $\langle\xi_k(t)\xi_k(t')^*\rangle=Q_k\delta(t-t')$，其稳态方差为

$$
\left\langle|\delta g_k|^2\right\rangle
=\frac{Q_k}{2(\gamma+Dk^2)}.
\tag{6.3}
$$

若 $Q_k$ 固定，扩散能压低非均匀涨落。但是，分子的随机跳跃也会产生噪声，因而真实反应–扩散系统中的 $Q_k$ 往往随 $D$ 改变。对理想的独立产生、降解和随机迁移分子，均匀稳态可以仍保持 Poisson 计数分布；不能从“有扩散”直接推出局部 Poisson 方差变小。

扩散对额外噪声、共同噪声和迁移噪声的作用需要分别计算。均匀的共同波动对应 $k=0$，空间平均尤其不能将它消除。

![线性格点模型中，扩散压低局部额外噪声，也减弱平均图案；二者竞争产生有限的局部精度最优值。](figures/05_diffusion_tradeoff.svg)

*Fig. 7. 扩散的两种效果。* 周期性 64 格点模型，格距 $a$，降解率 $\gamma$，产生率包含第 4 个余弦模态。各点受到固定强度、相互独立的加性白噪声；扩散在此仅为确定性耦合，没有额外迁移噪声。A：平均模式振幅和单点噪声标准差均相对于 $D=0$ 归一化。B：最大斜率处的 $\mathcal J_\mu(D)/\mathcal J_\mu(0)$ 出现有限最优值。C：强扩散使表达差异消失。该例只展示额外局部噪声的平均与信号平滑之间的竞争，并未计算完整胚胎模型的 $I(D)$。

### 6.3 局部精度的最优值不等于全局编码最优值 {#diffusion-global}

图中的余弦模式在空间重复。即使某个 $D$ 使最大斜率处的定位尺度最好，单个余弦表达仍然有不同位置取相同值的歧义。因此，图 7 的纵轴必须理解为局部 Fisher information，而不能改标为完整位置信息。

在实际网络里改变 $D$，均值、协方差、远处位置的可区分性和达到图案所需的时间都会改变。如果每次改变 $D$ 后重新优化其余参数，得到的是“该扩散条件下能实现的最好网络”；如果固定网络只改变 $D$，测量的则是该网络对扰动的响应。两种曲线回答不同问题。

## 7. 从随机调控动力学到可优化的编码通道 {#network-optimization}

### 7.1 同一调控函数在所有位置运行 {#stochastic-network}

一个机制模型需要解释 $P(\mathbf g\mid x)$ 从哪里来。作为一般的两级表达示意，令 $m_{\alpha i}$、$p_{\alpha i}$ 分别为位置 $i$ 的基因 $\alpha$ 的 mRNA 和蛋白质数，可以写成

$$
\begin{aligned}
\dot m_{\alpha i}
&=r_\alpha(\mathbf p_i,\mathbf c_i;\theta)
-\gamma_{m,\alpha}m_{\alpha i}+\xi^m_{\alpha i},\\
\dot p_{\alpha i}
&=k_\alpha m_{\alpha i}-\gamma_{p,\alpha}p_{\alpha i}
+\frac{D_\alpha}{a^2}\sum_{j\sim i}(p_{\alpha j}-p_{\alpha i})
+\xi^p_{\alpha i}.
\end{aligned}
\tag{7.1}
$$

这里的随机项代表由反应计量与倾向函数产生的涨落，并非可以独立任取的附加误差。迁移造成相邻位置一增一减，因此还会在噪声协方差中留下非对角项。各位置使用同一个 $r_\alpha$，差异来自局部输入、当前状态与空间交换。

一种有界的调控参数化是

$$
\begin{aligned}
r_\alpha&=\frac{r_{\max,\alpha}}{1+e^{-F_\alpha}},\\
F_\alpha&=b_\alpha
+\sum_\beta H_{\alpha\beta}\ln\left(1+\frac{p_\beta}{K_{\alpha\beta}}\right)\\
&\quad+\sum_\kappa H^c_{\alpha\kappa}\ln\left(1+\frac{c_\kappa}{K^c_{\alpha\kappa}}\right).
\end{aligned}
\tag{7.2}
$$

正、负 $H$ 表示激活与抑制，$K$ 决定响应的浓度范围。这只是说明如何给一类网络建立有限维参数空间；具体模型还需规定启动子近似、噪声机制、初态、几何与边界条件。

### 7.2 均值和协方差怎样由动力学给出 {#linear-noise}

把所有基因、核位置及保留的分子自由度合并成向量 $\mathbf y$。在一个稳定的平均轨迹附近，线性噪声近似写成

$$
\dot{\boldsymbol\mu}=\mathbf F(\boldsymbol\mu,t;\theta),
\qquad
d\delta\mathbf y=A(t)\delta\mathbf y\,dt+B(t)\,d\mathbf W_t,
\tag{7.3}
$$

其中 $A=\partial\mathbf F/\partial\mathbf y$ 在平均轨迹上计算。协方差 $C=\mathbb E[\delta\mathbf y\delta\mathbf y^{\mathsf T}]$ 满足

$$
\boxed{\dot C=AC+CA^{\mathsf T}+Q,\qquad Q=BB^{\mathsf T}.}
\tag{7.4}
$$

该式把反馈对涨落的放大或抑制、扩散引起的空间耦合和反应噪声的持续注入放到一起。非线性系统中“确定性轨迹等于精确随机均值”并不普遍成立；这里采用的是线性噪声展开的近似。若存在显著多稳态切换或低分子数效应，需要更完整的随机模拟或分布演化。

在读取时刻 $T$，取出单个核的蛋白质均值与对应的协方差子块 $\Sigma_i$，近似构成 $P(\mathbf g\mid x_i,T)$。全空间协方差必须在动力学中计算，但如果解码器只观察本核表达，计算信息时使用本核的边缘分布；不能把邻近核的表达也悄悄交给解码器。

### 7.3 优化目标和约束必须一起给出 {#objective-and-constraints}

于是完整计算链条是

$$
\theta
\longrightarrow
\{\boldsymbol\mu_i(T),\Sigma_i(T)\}
\longrightarrow P(\mathbf g\mid x_i,T)
\longrightarrow I(X;\mathbf G(T)).
\tag{7.5}
$$

优化问题为

$$
\boxed{
\theta^*=\operatorname*{arg\,max}_{\theta\in\mathcal C}
I(X;\mathbf G(T;\theta)).
}
\tag{7.6}
$$

物理约束 $\mathcal C$ 至少要说明合成速率和分子数、实际平均资源使用、反应时间尺度、读取时间、扩散和几何，以及允许的调控函数与强度。没有这些条件，增加分子数就能不断降低相对噪声，“最大化信息”本身不会给出有意义的生物学预测。

这里优化的是网络参数，不是在一个给定通道上自由改变位置先验。后者的 $\max_{P_X}I(X;\mathbf G)$ 才是通常所说的通道容量；本文的主要问题固定细胞核取样方式，优化编码机制。

### 7.4 一个真实网络的研究实例 {#sokolowski-2025}

[Sokolowski 等（2025）](https://pmc.ncbi.nlm.nih.gov/articles/PMC11725783/)把上述思路用于果蝇 gap gene 网络，在分子资源与发育时序约束下优化 50 余个调控参数。其 WT-like 解约为 4.2 bits，与实验的 4.2–4.3 bits 接近；同时比较了表达图案、网络相互作用与解码表现。

该模型的基准尺度包括 $D=0.5\,\mu\mathrm m^2/\mathrm s$、蛋白质寿命 10 min、mRNA 寿命 20 min、核间距 $8.5\,\mu\mathrm m$，平均资源使用约为最大表达的 20%。用前两个数估算 $\sqrt{D\tau_P}\simeq17\,\mu\mathrm m$，约为两个核间距；这是无反馈的尺度估计。

论文还报告了扩散的宽最优区间，以及三基因网络可以有较高总信息、却保留局部定位缺陷的结果。搜索中可以使用弱图案偏置，也考察无偏置的优化。因此结果应理解为**给定模型和实测约束下的优化预测**；既不是没有参数，也不是已经证明唯一的全局最优解。以上结果和约束见原文 [Box 1 与 Figs. 2–4](https://pmc.ncbi.nlm.nih.gov/articles/PMC11725783/)。

## 8. 高互信息还不能替代哪些发育要求 {#limits-and-predictions}

### 8.1 平均性能可以掩盖局部缺陷 {#mean-versus-worst}

在单峰、小噪声近似内，最大化互信息等价于压低 $\langle\ln\sigma_x\rangle$。这会惩罚误差较大的区域，但一个占比较小的区域，仍可能被其他区域较好的精度补偿。于是下面两个目标一般不同：

$$
\max_\theta I(X;\mathbf G),
\qquad
\min_\theta\max_{x\in\Omega}\sigma_x(x;\theta).
\tag{8.1}
$$

也可以针对发育任务规定损失，例如把跨越细胞命运边界的错误赋予更大代价，或要求 $P(|\hat X-X|>\epsilon)<\delta$。这些目标都需要明确“哪些位置的哪些错误最不能接受”。

![两个假设的误差分布具有相同的几何平均，却有不同的均方根误差和最坏位置误差。](figures/06_average_and_defects.svg)

*Fig. 8. 相同的平均对数精度，不保证相同的局部可靠性。* 蓝线误差恒为分析区间的 1.5%；橙线在 10% 区域为 6%，其余区域约为 1.29%，因此两者具有相同几何平均误差。B 比较几何平均、均方根和最坏位置误差。这里直接构造误差函数，展示式（4.4）的平均方式；没有指定完整 $P(\mathbf g\mid x)$，所以不能据此声称两种实际通道具有完全相同的精确互信息。

### 8.2 读出时间的不确定性会变成额外表达噪声 {#readout-time}

[Lucas 等（2018）](https://doi.org/10.1371/journal.pgen.1007676)用 MS2 活体转录报告发现，所研究的早期 *hb* 报告系统在每次间期约 3 min 内建立陡峭且位置稳定的边界。这与 Fig. 1 跨两个发育阶段的 30 min 不同；它提醒我们，实际读出还受到单次间期的可用时间限制。

在固定时刻计算的高信息，并不保证下游在一段时间内都能读取同样的位置。如果实际读出时刻为 $T+\delta t$，局部表达变化为

$$
\delta\mathbf g_{\mathrm{time}}
\simeq\partial_t\boldsymbol\mu(x,T)\,\delta t.
$$

若时间抖动独立于其他噪声，均值为零、方差为 $\sigma_t^2$，则小抖动下新增的协方差为

$$
\Sigma_{\mathrm{time}}
\simeq\sigma_t^2
(\partial_t\boldsymbol\mu)
(\partial_t\boldsymbol\mu)^{\mathsf T}.
\tag{8.2}
$$

它对位置读取的影响仍取决于方向：沿着空间切向的时间漂移尤其容易被误认成位置变化。若读出时刻本身也可被细胞估计，应该研究联合解码；不能简单把时间噪声视为不可消除的常数。

### 8.3 编码中的信息还需要可实现的解码器 {#implementable-decoder}

互信息衡量表达包含多少关于位置的信息，并不指定下游怎样实现 Bayes 后验。若下游输出 $Y$ 仅由本核当前表达 $\mathbf G$ 和独立内部噪声产生，形成 Markov 链 $X\to\mathbf G\to Y$，数据处理不等式要求

$$
I(X;Y)\leq I(X;\mathbf G).
\tag{8.3}
$$

有限的分子数、时间和反应步骤可能使下游只取出其中一部分信息。因此，一个完整解释还应将下游阈值、组合调控或时间积分写成机制，比较它们实际保留多少信息。

这也不意味着下游某个时刻的信息必然小于上游任意一个瞬时快照。若下游整合了输入历史或邻域，它的输入变量就应写成相应的历史和空间观测；对不同观测集合的 bits 数直接套用式（8.3），会造成错误比较。

### 8.4 怎样使优化原则产生可检验的预测 {#testable-predictions}

最有力的检验需要区分物理扰动与重新优化。在现有胚胎里改变扩散或母源剂量，检验的是固定调控网络对扰动的响应；删除一个基因后，把其余参数全部重新优化，检验的是该模型类别在替代资源配置下能够达到什么性能，不能直接当作急性突变表型。

可以依次提出三个可计算的预测。第一，在固定网络里改变源强和输入噪声，同时预测平均边界与完整后验如何变化。第二，改变分子预算后，分别计算均值、协方差与信息，而非只比较表达曲线是否相似。第三，在相同预算下比较多个网络，要求它们同时满足平均信息、局部大误差概率和读取时间稳定性，再检查哪些相互作用仍然反复出现。

这些比较把优化假说变成一组可被数据否定的约束。若不同物理可行网络都达到类似信息量，说明功能目标允许多种机制；若加入独立测得的限制后只有一小类结构保留下来，网络结构的预测才会更有辨识力。

---

## 技术附录 {#technical-appendices}

## 附录 A 有限区间梯度与缩放条件 {#appendix-a}

### A.1 用后端边界直接选取解的形式

稳态方程为 $c''-c/\lambda^2=0$。写成

$$
c(x)=A\cosh\frac{L-x}{\lambda}
+B\sinh\frac{L-x}{\lambda}.
$$

在 $x=L$ 处，$c'(L)=-B/\lambda=0$，所以 $B=0$。前端条件给出 $DA\sinh(L/\lambda)/\lambda=R$，即

$$
A=\frac{R\lambda}{D\sinh(L/\lambda)}
=\frac{R\tau}{\lambda\sinh(L/\lambda)}.
$$

这就得到式（1.2）。积分稳态方程，还能检验总的分子收支：

$$
0=D[c'(L)-c'(0)]-\frac1\tau\int_0^Lc(x)\,dx
\quad\Longrightarrow\quad
\int_0^Lc(x)\,dx=R\tau.
\tag{A.1}
$$

总存量等于输入通量乘平均寿命，与空间如何分布一致。

### A.2 指数近似怎样接近后端时失效

有限区间解可以改写成

$$
c(x)=\frac{R\tau}{\lambda}
\frac{e^{-x/\lambda}+e^{-(2L-x)/\lambda}}
{1-e^{-2L/\lambda}}.
\tag{A.2}
$$

第二个指数相对于第一个的比例是 $e^{-2(L-x)/\lambda}$。因此，除了 $L\gg\lambda$，还需要观察点离后端足够远。恰在后端，两个指数相等，梯度为零；用指数斜率估计该处的位置精度会出错。

### A.3 无量纲形式中的缩放

令 $z=x/L$、$\kappa=L/\lambda$，则

$$
c(z;L)=\frac{R\tau}{\lambda}
\frac{\cosh[\kappa(1-z)]}{\sinh\kappa}.
\tag{A.3}
$$

固定浓度阈值所对应的 $z_*$ 由 $\kappa$ 与无量纲阈值 $c_*\lambda/(R\tau)$ 共同决定。保持这两个组合不变是一组充分的缩放条件；只保持其中一个，通常不足以使完整图案重合。

## 附录 B Fisher information 与有效位置估计 {#appendix-b}

### B.1 Poisson 计数的精确 Fisher information

对式（2.2）取自然对数并求导：

$$
\partial_x\ln P(n\mid x)
=\left(\frac n{m(x)}-1\right)m'(x).
$$

利用 $\mathbb E[N]=m$、$\operatorname{Var}(N)=m$：

$$
\mathcal J(x)
=\frac{[m'(x)]^2}{m(x)^2}
\mathbb E[(N-m)^2]
=\frac{[m'(x)]^2}{m(x)}.
\tag{B.1}
$$

这个 Fisher information 不需要把 Poisson 近似为高斯。指数梯度下为 $m/\lambda^2$。但单次计数可能为零，此时形式上的对数反演 $\hat x=\lambda\ln(m_0/N)$ 不存在；在有限区间内需要边界处理或先验。因此，精确 Fisher information 的存在不代表单次估计量一定无偏或达到下界。

### B.2 固定协方差时的最优局部线性读出

令 $\mathbf v=\boldsymbol\mu'(x_0)$，局部模型为

$$
\mathbf G-\boldsymbol\mu(x_0)
=\mathbf v\,\delta x+\boldsymbol\eta,
\qquad \boldsymbol\eta\sim\mathcal N(0,\Sigma).
$$

取线性估计 $\widehat{\delta x}=\mathbf w^{\mathsf T}[\mathbf G-\boldsymbol\mu(x_0)]$，局部无偏要求 $\mathbf w^{\mathsf T}\mathbf v=1$。最小化方差 $\mathbf w^{\mathsf T}\Sigma\mathbf w$，得到

$$
\mathbf w_*=
\frac{\Sigma^{-1}\mathbf v}{\mathbf v^{\mathsf T}\Sigma^{-1}\mathbf v},
\qquad
\operatorname{Var}(\widehat{\delta x})=
\frac1{\mathbf v^{\mathsf T}\Sigma^{-1}\mathbf v}.
\tag{B.2}
$$

这个结果具体说明了局部界怎样达到：不同基因按斜率和完整协方差加权。它要求已知展开位置附近的分支，不能作为全局歧义的解决办法。

### B.3 高斯协方差变化时的附加项

写 $\mathbf y=\mathbf g-\boldsymbol\mu(x)$，则

$$
\ln P=-\frac K2\ln(2\pi)-\frac12\ln\det\Sigma
-\frac12\mathbf y^{\mathsf T}\Sigma^{-1}\mathbf y.
$$

利用 $\partial_x\ln\det\Sigma=\operatorname{Tr}(\Sigma^{-1}\Sigma')$ 和 $\partial_x\Sigma^{-1}=-\Sigma^{-1}\Sigma'\Sigma^{-1}$，score 为

$$
\partial_x\ln P
=\boldsymbol\mu'^{\mathsf T}\Sigma^{-1}\mathbf y
+\frac12\left[
\mathbf y^{\mathsf T}\Sigma^{-1}\Sigma'\Sigma^{-1}\mathbf y
-\operatorname{Tr}(\Sigma^{-1}\Sigma')\right].
\tag{B.3}
$$

第一项是关于 $\mathbf y$ 的线性函数，第二项是中心化二次函数。零均值高斯分布的三阶矩为零，二者的交叉贡献消失；再用高斯四阶矩公式，得到式（3.4）。这一额外项描述分布宽度和方向的空间变化；当它不可忽略时，简单的线性均值解码不再利用全部局部信息。

### B.4 共同波动与位置混淆

考虑一个小的共同扰动 $u$，其表达方向为 $\mathbf b$：

$$
\mathbf G=\boldsymbol\mu(x)+\mathbf b\,u+\boldsymbol\epsilon,
\quad \operatorname{Var}(u)=s_u^2,
\quad \operatorname{Cov}(\boldsymbol\epsilon)=\Sigma_0.
$$

若 $u$ 与其他噪声独立，边缘化后 $\Sigma=\Sigma_0+s_u^2\mathbf b\mathbf b^{\mathsf T}$。Sherman–Morrison 公式给出均值贡献

$$
\mathcal J_\mu=
\mathbf v^{\mathsf T}\Sigma_0^{-1}\mathbf v
-\frac{s_u^2(\mathbf v^{\mathsf T}\Sigma_0^{-1}\mathbf b)^2}
{1+s_u^2\mathbf b^{\mathsf T}\Sigma_0^{-1}\mathbf b}.
\tag{B.4}
$$

共同噪声和位置信号在去除基底噪声尺度后的空间里越重合，损失越大。若 $\mathbf b=b\mathbf v$，则 $\mathcal J_\mu=J_0/(1+s_u^2b^2J_0)$；即使基底测量越来越精确，也留下位置等效扰动 $s_u|b|$ 的尺度。这里 $\Sigma_0$ 正定，方向与协方差都取局部常数。

## 附录 C 互信息、位置后验与离散化 {#appendix-c}

### C.1 小噪声公式中的几何平均

若位置先验均匀，$h(X)=\log_2W$。若给定观测的后验近似为单高斯，方差为 $s^2(\mathbf g)$，则

$$
h(X\mid\mathbf G)
\simeq\int d\mathbf g\,P(\mathbf g)
\frac12\log_2[2\pi e\,s^2(\mathbf g)].
\tag{C.1}
$$

小噪声下，每个典型观测集中在其真实位置附近，若局部宽度平滑变化，可以把对表达观测的平均近似为对位置的平均，得到式（4.3）。这个替换和高斯后验都是近似；边界截断、多峰后验和大误差都会破坏它。

两个等权、相隔很远的窄高斯峰就是一个反例。每个峰宽为 $s$，互相几乎不重叠，则后验熵约为单峰高斯熵再加 1 bit。只用局部峰宽会漏掉这个分支选择的不确定性。与此同时，后验总体方差还包含峰间距，不能用任一个峰的 $s^2$ 代替。

### C.2 一般先验与不同的优化问题

若先验不是均匀分布，局部近似为

$$
I\simeq h(P_X)+\frac12\int dx\,P_X(x)
\log_2\frac{\mathcal J(x)}{2\pi e}.
\tag{C.2}
$$

假设通道固定，且这个近似成立，变分优化 $P_X$ 会得到 $P_X^*(x)\propto\sqrt{\mathcal J(x)}$。这只是固定通道上选择输入分布的结果；不能转而断言胚胎应该改变细胞核密度来符合它，也不能把它和固定先验下优化网络参数混为一谈。

### C.3 核位置离散化不会改变推断问题的定义

若计算采用 $N$ 个位置，取概率 $p_i$，则

$$
P(\mathbf g)=\sum_{j=1}^Np_jP(\mathbf g\mid x_j),
$$

$$
I=\sum_{i=1}^Np_i\int d\mathbf g\,P(\mathbf g\mid x_i)
\log_2\frac{P(\mathbf g\mid x_i)}{\sum_jp_jP(\mathbf g\mid x_j)}.
\tag{C.3}
$$

这时 $I\leq-\sum_i p_i\log_2p_i$。若位置间距远小于噪声造成的位置宽度，离散计算可以逼近连续结果；若噪声远小于核间距，连续小噪声公式可能超过 $\log_2N$，应改用离散表达，而不是接受一个不可能的核标签信息值。

图 5 使用均匀格点与 Gauss–Hermite 求积。对于真实位置 $x_i$，计算

$$
\mathcal D(x_j\mid x_i)
=\int d\mathbf g\,P(\mathbf g\mid x_i)
P(x_j\mid\mathbf g).
\tag{C.4}
$$

每列概率和为 1，绘图时除以格点间距得到密度。$\mathcal D$ 描述在真实位置下，对表达噪声平均后得到的后验；它并不自动等于某个点估计器输出的频率分布。MAP 估计、后验均值以及从后验抽样，会给出不同的点估计误差。

### C.4 从实验样本估计信息时需要控制什么

表达数据需要先限定发育时刻、分析区间和取样方式。同一位置上的胚胎间波动既可能来自真实生物差异，也可能来自成像增益、背景和空间配准；它们对 $P(\mathbf g\mid x)$ 的估计有直接影响。若用每条表达曲线本身进行过强的对齐，还可能消除原本需要测量的位置变异。

有限样本下，直接将高维表达分箱会产生偏差；高斯拟合虽然降低估计维度，也引入分布假设。应检查样本量外推、重采样误差、协方差估计稳定性，并将训练编码分布与检验解码表现的数据分开。相关方法及其在 gap genes 上的比较见 [Tkačik 等（2015）](https://pmc.ncbi.nlm.nih.gov/articles/PMC4286692/)。一个信息数字只有连同这些统计定义才可比较。

## 附录 D 有限时间内的噪声平均 {#appendix-d}

设 $g(t)$ 为均值恒定的平稳过程，涨落相关函数为 $C_g(s)=\langle\delta g(t+s)\delta g(t)\rangle$。定义时长 $T$ 的平均 $\bar g_T=T^{-1}\int_0^Tg(t)\,dt$，则

$$
\begin{aligned}
\operatorname{Var}(\bar g_T)
&=\frac1{T^2}\int_0^T\!dt\int_0^T\!dt'\,C_g(t-t')\\
&=\frac2{T^2}\int_0^T(T-s)C_g(s)\,ds.
\end{aligned}
\tag{D.1}
$$

若 $C_g(s)=\sigma_g^2e^{-|s|/\tau_c}$，精确结果为

$$
\operatorname{Var}(\bar g_T)=
\frac{2\sigma_g^2}{T^2}
\left[T\tau_c-\tau_c^2(1-e^{-T/\tau_c})\right].
\tag{D.2}
$$

当 $T\ll\tau_c$，方差趋向 $\sigma_g^2$，几乎没有独立平均；当 $T\gg\tau_c$，方差近似为 $2\sigma_g^2\tau_c/T$。一般情形可用积分相关时间 $\tau_{\mathrm{int}}=\sigma_g^{-2}\int_0^\infty C_g(s)\,ds$ 替代 $\tau_c$，前提是该积分存在且窗口内平均信号近似不变。

若 $C_g$ 还包含不随 $s$ 衰减的共同分量 $\sigma_{\mathrm{ext}}^2$，式（D.1）会留下同样大小的方差，而不是随 $T$ 消失。若平均表达在窗口内明显漂移，则又会引入动态读出的偏差，不能直接使用平稳公式。

## 附录 E 从反应噪声到协方差演化 {#appendix-e}

对于反应通道 $r$，令计量变化向量为 $\boldsymbol\nu_r$，反应倾向为 $a_r(\mathbf y)$。在化学 Langevin 与线性噪声近似中，沿平均轨迹的漂移和噪声矩阵为

$$
\mathbf F=\sum_r\boldsymbol\nu_r a_r,
\qquad
Q=\sum_r\boldsymbol\nu_r\boldsymbol\nu_r^{\mathsf T}a_r.
\tag{E.1}
$$

产生与降解通道都贡献正的单点方差。对于从核 $i$ 到核 $j$ 的迁移，$\boldsymbol\nu=-\mathbf e_i+\mathbf e_j$，所以 $\boldsymbol\nu\boldsymbol\nu^{\mathsf T}$ 的 $ij$、$ji$ 分量为负；同一次迁移会使源位置减少一个分子，目标位置增加一个分子。这解释了扩散噪声为何不能当作各点独立的白噪声处理。

对线性化随机微分方程应用 Itô 乘法公式：

$$
d(\delta\mathbf y\delta\mathbf y^{\mathsf T})
=(d\delta\mathbf y)\delta\mathbf y^{\mathsf T}
+\delta\mathbf y(d\delta\mathbf y)^{\mathsf T}
+(d\delta\mathbf y)(d\delta\mathbf y)^{\mathsf T}.
$$

取平均，使用 $\mathbb E[d\mathbf W_t\,d\mathbf W_t^{\mathsf T}]=I\,dt$，就得到式（7.4）。若 $A$ 固定且所有特征值实部为负，稳态协方差满足 $AC+CA^{\mathsf T}+Q=0$。不稳定的平均态、明显的分支切换，以及有色外源噪声，都需要修改这一步描述。

对于图 7 的周期格点，离散 Laplacian 的正特征值为 $\Lambda_j=4\sin^2(\pi j/N)/a^2$。设噪声强度在所有格点相同且独立，$Q=qI$，稳态单点方差为

$$
\operatorname{Var}(g_i)
=\frac q{2N}\sum_{j=0}^{N-1}\frac1{\gamma+D\Lambda_j}.
\tag{E.2}
$$

第 $k$ 个平均模式的振幅则乘以 $1/(\gamma+D\Lambda_k)$。因此最大斜率处的 Fisher information 相对于 $D=0$ 的比值可直接由信号振幅平方与单点方差计算。$D\to\infty$ 时，零模态的噪声仍然保留，非零平均模式却消失，所以局部信息最终降低。

## 附录 F 一个可以继续计算的最小问题 {#appendix-f}

完整网络优化包含很多参数。要分离机制，可以先固定长度为 $L$ 的一维格点、两个母源输入、两到四个输出基因，以及相同的总平均分子预算。各位置采用相同的有界产生率，包含降解、蛋白交换及由反应得到的噪声；在固定时间 $T$ 读取本核的表达向量。

对每组参数，分别记录四种量：完整 $I(X;\mathbf G)$、局部 $\mathcal J_\mu(x)$、完整后验的大误差概率，以及读取窗口内的图案漂移。随后改变扩散、共同输入噪声或基因数，并区分固定参数扰动和重新优化。这样的比较能够检验：提高平均局部精度时，是否也消除了远处分支；压低独立噪声后，共同噪声是否成为瓶颈；总信息相近时，局部缺陷是否仍然存在。

这组问题已经保留了统计物理上的主要竞争：有限分子数限制读出，扩散改变空间信号与涨落，调控网络分配表达范围和分子资源，而全局后验决定这些局部变化能否组成一个可靠的位置编码。

## 参考文献与复现 {#references}

### 按问题选择背景阅读 {#biological-reading}

果蝇前后轴图案并不只对应一种统计问题。综述提供系统背景，而具体实验分别检验梯度的形成、启动子的响应，以及多个基因能否共同指定位置。

| 想理解的问题 | 阅读入口 | 与本笔记的联系 |
|---|---|---|
| Bicoid 系统包含哪些生物学过程？ | [Porcher & Dostatni, 2010](https://doi.org/10.1016/j.cub.2010.01.026)，系统综述 | 先区分母源输入、转录调控与图案的鲁棒性 |
| 梯度怎样形成，如何区分机制？ | [Durrieu et al., 2018](https://doi.org/10.15252/msb.20188355)，蛋白质年龄测量 | 第 1 节：浓度图案之外还需要动态观测 |
| 下游转录怎样在有限时间内形成边界？ | [Lucas et al., 2018](https://doi.org/10.1371/journal.pgen.1007676)，活体转录成像 | Fig. 1 与第 8.2 节：输入、响应及读出时间 |
| 有噪声的表达能指定多精确的位置？ | [Gregor et al., 2007](https://pmc.ncbi.nlm.nih.gov/articles/PMC2253670/)；[Dubuis et al., 2013](https://doi.org/10.1073/pnas.1315642110)；[Petkova et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6526179/) | 从单一梯度读取，走向多基因的信息测量与位置后验 |
| 给定资源和动力学，哪些网络能产生高信息编码？ | [Sokolowski et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11725783/) | 第 7 节：受约束的网络优化 |

### 引用文献

1. Dubuis, J. O., Tkačik, G., Wieschaus, E. F., Gregor, T. & Bialek, W. (2013). [Positional information, in bits](https://doi.org/10.1073/pnas.1315642110). *PNAS* **110**, 16301–16308. 位置信息的实验测量，以及连续表达与二值表达的比较。
2. Tkačik, G., Dubuis, J. O., Petkova, M. D. & Gregor, T. (2015). [Positional Information, Positional Error, and Readout Precision in Morphogenesis: A Mathematical Framework](https://pmc.ncbi.nlm.nih.gov/articles/PMC4286692/). *Genetics* **199**, 39–59. 信息、位置误差与估计方法的关系。
3. Petkova, M. D., Tkačik, G., Bialek, W., Wieschaus, E. F. & Gregor, T. (2019). [Optimal decoding of cellular identities in a genetic network](https://pmc.ncbi.nlm.nih.gov/articles/PMC6526179/). *Cell* **176**, 844–855.e15. Bayes 位置解码及其与下游图案的比较。
4. Sokolowski, T. R., Gregor, T., Bialek, W. & Tkačik, G. (2025). [Deriving a genetic regulatory network from an optimization principle](https://pmc.ncbi.nlm.nih.gov/articles/PMC11725783/). *PNAS* **122**, e2402925121. 在物理约束下优化随机空间调控网络。
5. Gregor, T., Tank, D. W., Wieschaus, E. F. & Bialek, W. (2007). [Probing the limits to positional information](https://pmc.ncbi.nlm.nih.gov/articles/PMC2253670/). *Cell* **130**, 153–164. Bicoid 浓度的活胚胎测量与位置读取精度；Fig. 2 的原图来源。

6. Porcher, A. & Dostatni, N. (2010). [The Bicoid Morphogen System](https://doi.org/10.1016/j.cub.2010.01.026). *Current Biology* **20**, R249–R254. 系统背景综述。
7. Lucas, T., Tran, H., Perez Romero, C. A., et al. (2018). [3 minutes to precisely measure morphogen concentration](https://doi.org/10.1371/journal.pgen.1007676). *PLoS Genetics* **14**, e1007676. 早期 *hb* 转录响应及其时间尺度；Fig. 1 的原图来源。
8. Durrieu, L., Kirrmaier, D., Schneidt, T., et al. (2018). [Bicoid gradient formation mechanism and dynamics revealed by protein lifetime analysis](https://doi.org/10.15252/msb.20188355). *Molecular Systems Biology* **14**, e8355. 以蛋白质年龄作为梯度形成机制的额外检验。

### 图源与计算复现

六幅教学图由本目录的 [generate_figures.py](figures/generate_figures.py) 生成，依赖 NumPy、SciPy 和 Matplotlib。图 5 对表达噪声作数值求积；图 6 计算单通道互信息；图 7 使用线性系统的解析稳态协方差。均未使用或生成胚胎实验数据。两幅文献原图的下载地址、许可说明与文件校验值记录在 [paper_figure_sources.json](figures/paper_figure_sources.json)。

在仓库根目录运行：

```bash
python3 notes/qbio_course/positional_information/figures/generate_figures.py
```
