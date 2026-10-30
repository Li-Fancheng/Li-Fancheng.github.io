> **Hierarchical Learning · 理论 03**　[← 理论 02 · 梯度可达性](/notes/hierarchical_data/gradient_accessibility/) · [总论](/notes/hierarchical_data/) · [理论 04 · 跨层反馈 →](/notes/hierarchical_data/coupled_feature_learning/)

**本篇问题：** 发现一个非线性中间表示，为什么能比直接拟合最终目标节省样本和计算？深度在这个过程中究竟做了什么？

考虑 $f_*(x)=g(p(x))$。如果 $p$ 是 $k$ 次多项式、$g$ 是 $q$ 次多项式，直接看 $x$，目标最高可达 $kq$ 次；若先得到标量 $p(x)$，后续只需解一个一维问题。困难在于训练器没有 $p$ 的标签，也不能预先知道这个分解。以下文献分别证明特定算法可以提取什么表示，以及这种表示如何改变后续复杂度。它们给出了比表达能力更进一步的结果，同时尚未构成一般深网训练的统一理论。

## 1. 四个不同层次的深度优势 {#four-advantages}

### 1.1 表达、统计、计算和轨迹 {#levels-of-claim}

表达优势问“是否存在小网络”；统计优势问“给定函数类需要多少样本”；计算优势问“是否有资源可控的算法找到预测器”；训练动力学则问“指定 GD/SGD 如何走到那里”。一个组合函数可以有短表达，却因低阶相关全为零而缺少梯度入口。[Staircase 理论](/notes/hierarchical_data/gradient_accessibility/)已经给出这种区别。

一种常见的核方法比较来自旋转对称输入：固定核不能根据任务旋转其表示，一个次数 $K$ 的方向可能分布在约 $d^K$ 维的基底空间中。学习特征可以让任务落入更小空间。但“核方法有某个下界”“某两层训练协议较差”“所有浅层网络都难以表达”是不同命题，后文只在文献确实证明的范围内使用。

### 1.2 降低维数不等于随意压缩 {#sufficiency}

若已经得到 $z=\phi(x)$，只有当 $f_*(x)\simeq\widetilde g(z)$ 时，丢弃其他坐标才不会损害预测。需要同时衡量：表示对任务是否充分、是否能从样本学到、误差经后续非线性放大多少、形成表示花了多少资源。一个对训练标签非常相关的标量也可能只是记忆，不能因此认定它恢复了生成潜变量。

本篇的模型通过具体目标分解和误差界回答这些问题。RHM 中的[有效状态与任务充分性](/notes/hierarchical_data/generative_learning/#effective-states)是离散版本：不同叶节点实现可以被归并，前提是后续任务只需要共同的父状态。

### 1.3 深度节省的不是次数本身，而是重复搜索 {#cost-accounting}

考虑一个简单但有代表性的分解：$p(x)$ 是难以从高维输入中发现的非线性标量，最终目标为 $f_*(x)=p(x)+\lambda p(x)^2$。即使 $p$ 是二次函数，直接把最终目标展开到原始坐标上也包含四次项。如果先恢复 $\widehat p$，后续的外层拟合只在一个标量上进行，不必为不同的 $p^j$ 分别在高维特征空间中搜索。

但必须计入建立表示的成本。若用独立样本集先学 $\widehat p$，再学 $\widehat g$，总样本预算是 $P_{\rm rep}+P_{\rm outer}$；若两阶段复用数据，不能简单相加，却须分析自适应选择引入的依赖。假设 $p(x)$ 和 $\widehat p(x)$ 几乎处处落在同一区间，且外层 $g$ 在该区间是 $L_g$-Lipschitz，则一个基本误差分解为

$$
\mathbb E[(g(p(x))-\widehat g(\widehat p(x)))^2]
\le 2L_g^2\mathbb E[(p-\widehat p)^2]
+2\mathbb E[(g(\widehat p)-\widehat g(\widehat p))^2].
$$

第一项是表示误差通过外层传播，第二项是给定表示后的学习误差。对于无界输入上的多项式 $g$，不能无条件使用全局 Lipschitz 常数，需要限制到高概率有界区域或直接控制高阶矩。这个技术细节也有物理意义：若外层对某些输入极其敏感，一个均方误差看似很小的中间表示，仍可能造成较大的最终预测误差。

因此真正有用的比较对象不是“已知 $p$ 的一维拟合”对“未知 $p$ 的高维拟合”，而是**发现 $p$ 的成本加上利用 $p$ 的成本**，再对照一个直接端到端方法的成本。已有结构若被架构硬编码、被额外标签提供，或者由数据随机特征免费保留，这三个问题的起点本就不同。

### 1.4 高阶特征空间为何容易昂贵 {#polynomial-dimension}

在各向同性输入上，$k$ 次齐次多项式的对称系数空间维数为

$$
D_k=\binom{d+k-1}{k}\asymp d^k
\qquad(k\ \text{固定},\ d\to\infty).
$$

这个计数解释了为何旋转不变的固定核往往要面对大量等价的多项式方向：不预先知道目标落在哪个方向时，同一个次数层中有约 $d^k$ 个候选系数。高维困难并不直接来自“把 $k$ 次函数写在纸上很长”，而是相关信号分散在一个很大的对称表示空间中。

不过 $D_k$ 是维数计数，**不是每种核方法都必须使用 $d^k$ 个样本的普适下界**。样本需求还取决于核特征值、目标系数的分布、精度、噪声与允许的算法。深层方法的实质优势是通过先学一个任务相关子空间，让后续高次操作只在少量坐标上发生，而不是纯粹减少网络参数个数。

## 2. Nichani 2023：一次梯度步提取一个非线性核投影 {#nichani}

### 2.1 固定第一层也能形成新的可训练表示 {#kernel-projection}

Nichani、Damian 与 Lee 研究三个权重层：

$$
f(x)=\frac1{m_1}a^\top\sigma_1(W\sigma_2(Vx)+b).
$$

最内层 $V$ 为固定随机特征，$W$ 从零开始；先在第一份数据上只更新 $W$ 一步，再在独立数据上训练输出权重 $a$，$V,b$ 保持固定。设 $h_0(x)=\sigma_2(Vx)$，在无限特征近似下定义

$$
K(x,x')=\mathbb E_v[\sigma_2(v^\top x)\sigma_2(v^\top x')],\quad
(\mathcal K f_*)(x)=\mathbb E_{x'}[K(x,x')f_*(x')].
$$

初始预测接近零时，第一步 $\Delta w_j\propto a_j\sigma_1'(b_j)\widehat{\mathbb E}[y h_0(x)]$，所以新的预激活沿着 $\mathcal Kf_*(x)$ 变化。可训练外层接着拟合 $q(\eta\mathcal Kf_*)$。学到的不是一个原始线性投影，而是任务相关的非线性函数。[正文 §2–3](https://papers.nips.cc/paper_files/paper/2023/hash/236b6a814a1d2c0ff504ca7bf380f7ff-Abstract-Conference.html) [[1]](#ref-1)。

### 2.2 低频部分为什么可能足以重建高频目标 {#low-to-high}

主定理的误差结构可以概括为

$$
\|f_{\rm learned}-f_*\|_2^2
\lesssim\|q(\eta\mathcal Kf_*)-f_*\|_2^2
+\widetilde O\!\left(\frac{\|q\|_{1,\infty}^2\|\mathcal Kf_*\|_2^{-2}}
{\min(P,m_1,m_2)}\right)
+\text{外层近似与估计误差}.
$$

这是保留机制的概写，常数、正则化和外层光滑性条件见原定理。第一项要求目标可由该投影再组合；第二项说明先提取它需要多强的信号。固定核回归要恢复目标的全部谱分量，包括小特征值上的分量；这里先抓住能被核突出显示的分量，再通过可训练非线性生成其余部分。两者使用同一个初始核，却解决不同任务。

在 $f_*(x)=g_*(x^\top Ax)$ 的适当高维、扩散二次型模型中，核投影近似正比于 $x^\top Ax$，于是可用约 $\widetilde O(d^4)$ 的样本和宽度提取它，维数指数不随外层函数的多项式次数逐次增加。这个量级是论文相应设定中的保证，不是二次特征恢复的最优极限。

论文还给出特定旋转二次型加 ReLU 外层的逼近下界：若两层网络要达到随维数缩小的误差，其宽度或最大参数必须超多项式增长。这个下界针对构造的函数族、误差区间和参数控制，不能推广成“所有非线性组合目标对两层网络都指数难”。

### 2.3 把第一次梯度步真正代入后，会出现什么 {#one-step-derivation}

令初始 $W=0$，第 $j$ 个中间神经元的预激活为 $b_j$。平方损失对它的输入行 $w_j$ 的梯度，在输出误差近似 $y$ 的初期满足

$$
\Delta w_j\simeq
\frac{\eta a_j\sigma_1'(b_j)}{m_1}\,
\frac1P\sum_{\mu=1}^P y_\mu h_0(x_\mu),
\qquad h_0(x)=\sigma_2(Vx).
$$

将这个更新与一个新输入的固定随机特征相乘：

$$
\Delta w_j^\top h_0(x)
\simeq\frac{\eta a_j\sigma_1'(b_j)}{m_1P}
\sum_{\mu=1}^P y_\mu
\underbrace{\langle h_0(x_\mu),h_0(x)\rangle}_{K_{m_2}(x_\mu,x)}.
$$

右侧正是经验标签加权的随机特征 Gram 平均。按固定特征宽度归一化，即令 $m_2^{-1}K_{m_2}(x,x')$ 逼近第 2.1 节的 $K(x,x')$，并将相应常数吸收到有效步长后，它趋向 $\mathcal K f_*$；再由 $a_j,b_j$ 的分布和外层非线性变成可调的标量函数族。这给“梯度提取特征”一个具体含义：新表示不是任意猜测的潜变量，而是**目标经过当前核算子过滤后的函数**。

若核的本征函数为 $\phi_n$，$\mathcal K\phi_n=\lambda_n\phi_n$，且 $f_*=\sum_n c_n\phi_n$，则

$$
\mathcal Kf_*=\sum_n\lambda_nc_n\phi_n.
$$

小本征值上的原始目标成分会被压低，甚至完全消失。这个做法能成功，必须因为被保留下来的低阶或大本征值部分足以通过后续非线性恢复目标；若目标的重要信息全在核的零空间里，外层再强也无法从一个常数表示中恢复它。所谓低频到高频，不是逆转核算子的谱损失，而是借助目标的组合关系**重新生成**高频函数。

### 2.4 何时低频种子真的足够 {#kernel-seed-example}

一个说明性构造是 $f_*(x)=p(x)+\lambda[p(x)^2-\mathbb E p(x)^2]$，其中 $p(x)$ 为归一化的非线性特征。假定当前核滤波对目标的主要保留分量正比于 $p$，即 $\mathcal Kf_*\simeq c\,p$ 且 $c\ne0$。那么后续仅需用标量函数

$$
q(z)=\frac zc+\lambda\left[\left(\frac zc\right)^2-\mathbb E p^2\right]
$$

就能重建目标的线性和二次组合。即使 $p$ 是原始输入上的二次特征，$p^2$ 含四次项，外层仍不必再次独立搜索全部四次张量方向。

这里的假设 $\mathcal Kf_*\simeq c\,p$ 是**需要验证的统计条件**，不是组合函数天然具备的事实。特别是 $p^2$ 也可能向与 $p$ 无关的低阶方向投影，或 $p$ 自身的核本征值很小。一个有判别力的实验是直接估计 $\mathcal Kf_*$ 与目标潜变量 $p$ 的重合，并同时检查给定这个投影后的最优标量读出误差。只测训练后预测准确度，无法判断成功来自预期的种子还是其他偶然的核方向。

## 3. Wang 2024：反复提取中间特征，改进一次步的代价 {#wang}

### 3.1 目标的非退化低阶分量 {#polynomial-assumptions}

Wang、Nichani 与 Lee 考虑 $x\sim\mathcal N(0,I_d)$，$f_*=g\circ p$。$p$ 是 $k$ 阶纯 Hermite 多项式，写成很多独立、正交低维块的均衡和；块数随 $d$ 增长，使 $p$ 具有近 Gaussian 行为。$g$ 次数 $q$ 固定、系数有界，并要求

$$
\mathbb E_{Z\sim\mathcal N(0,1)}g'(Z)=\Theta(1).
$$

这一条件让目标在 $p$ 方向保留一次有效信号。近似 Stein 关系将组合目标的低阶投影联系到 $p$ 本身；若该系数消失，正文给出的学习路线不能原样使用。所有 $k$ 次多项式都能学习，并不是这篇定理的意思。[原文 Assumptions 1–6](https://proceedings.iclr.cc/paper_files/paper/2024/hash/af536fa4312b5b2e5c541156cd31e93b-Abstract-Conference.html) [[2]](#ref-2)。

### 3.2 瓶颈表示的两阶段训练 {#polynomial-training}

学生具有标量瓶颈与残差连接，可写成

$$
\nu(x)=u^\top\sigma_1(Vx+s),\qquad
f_\theta(x)=\nu(x)+c^\top\operatorname{ReLU}(a\nu(x)+b).
$$

$V,s,a,b$ 固定；先对 $u$ 多步梯度下降，以得到 $\nu\simeq c_0p$，再在独立数据上训练 $c$ 拟合外层。训练器没有 $p$ 的标签，第一阶段仍用最终标签 $y$，依靠组合投影的非退化性提取 $p$。与 Nichani 的一次步相比，反复更新能更有效累计中间信号。

对任意小的固定 $\alpha>0$，原定理以 $P\ge d^{k+3\alpha}$、第一模块宽度至少 $d^{k+\alpha}$ 等资源达到 $\widetilde O(d^{-\alpha})$ 风险。常见的“约 $d^k$”是概述其维数指数，不能读成在任意消失误差下严格只有 $d^k$ 乘对数。对 $k=2$，合适的迹零、Frobenius 范数归一化且算子范数 $O(d^{-1/2})$ 的二次型，将先前约 $d^4$ 的保证改进到接近 $d^2$ 的指数。

相应旋转对称核方法面对次数 $kq$ 的目标会承担约 $d^{kq}$ 的困难。这里的改进来自先学 $p$ 再拟合 $g$，同时强烈依赖目标的扩散结构、Gaussian 输入、瓶颈架构和逐层训练。它不是普通端到端 SGD 自动发现任意组合分解的定理。

### 3.3 Stein 条件可以用两个外层函数直接比较 {#stein-worked}

若 $Z\sim\mathcal N(0,1)$，适当可积、可微的 $g$ 满足 Gaussian 分部积分恒等式

$$
\mathbb E[Zg(Z)]=\mathbb E[g'(Z)].
$$

这说明 $\mathbb E[g'(Z)]$ 实际是 $g$ 在一阶 Hermite 方向上的投影系数。对 $g_+(Z)=Z+\lambda(Z^2-1)$，无论 $\lambda$ 取何固定值，$\mathbb E[g_+'(Z)]=1$；即使最终目标含高次组合，仍留下与中间指标方向对齐的一阶种子。相反，$g_0(Z)=Z^2-1$ 满足 $\mathbb E[g_0'(Z)]=0$，一阶投影完全消失。

这并不意味着 $g_0(p(x))$ 绝对不能学习；只是依赖一阶有效投影提取 $p$ 的这条路线失去了主要信号。其他高阶统计或利用不同监督的算法可能仍有效。外层函数只做一个小改动，所需学习机制便可能改变，这比单纯说“两个目标都具有相同组合深度”更能揭示训练难度。

实际 $p(x)$ 在该论文中只是渐近接近 Gaussian，因此必须控制用真实 $p(x)$ 代替理想 $Z$ 时的近似误差。扩散、正交块和归一化条件正是为这一近似服务；不能把上述一维 Stein 恒等式直接当成任意高维非线性 $p$ 的精确定理。

## 4. Dandi 2025：潜在搜索空间逐次降低维数 {#dandi-depth}

### 4.1 从 $d$ 到 $d^\varepsilon$ 再到标量 {#sight-model}

Dandi、Pesce、Zdeborová 与 Krzakala 引入 SIGHT：

$$
z_*=W_*x\in\mathbb R^{d^\varepsilon},\quad
h_* =d^{-\varepsilon/2}a_*^\top p_k(z_*),\quad
y=g_*(h_*),\qquad 0<\varepsilon<1,
$$

其中 $W_*$ 行正交，$p_k$ 逐坐标作用。输入为各向同性 Gaussian，层级藏在目标中。与固定潜在维数的 multi-index 不同，中间空间维数也随 $d$ 增长。更深的 MIGHT 将这类结构递归组合。[NeurIPS 2025 正文 §1–3](https://papers.nips.cc/paper_files/paper/2025/hash/eb3878c1dcbfff9ee95d5d033e5f5942-Abstract-Conference.html) [[3]](#ref-3)。

其核心机制可表示为

$$
\underbrace{d}_{\text{输入搜索空间}}
\longrightarrow\underbrace{d^\varepsilon}_{\text{学到的潜在子空间}}
\longrightarrow\underbrace{1}_{\text{非线性指标}}.
$$

第二步的次数仍为 $k$，但统计对象已在 $d^\varepsilon$ 维而非 $d$ 维中。由此学习一个 $k$ 次结构的成本从以 $d^k$ 为尺度，变成以 $d^{k\varepsilon}$ 为尺度。表示没有消除非线性，而是改变了非线性搜索所处空间。

### 4.2 严格算法比“逐层 GD”更具体 {#sight-theorem}

浅 SIGHT 的证明采用：第一层在相关损失上进行逐神经元球面 SGD、用独立批次；接着重新初始化部分权重，对第二层做一次经过特征协方差预条件的梯度步；最后 ridge 拟合读出。对均匀目标权重、外层信息指数 1、内部多项式有非零二阶 Hermite 分量，以及激活的解析性、Hermite 对齐和额外消失矩条件，得到三个阶段的样本指数：

$$
P_1\gtrsim d^{1+\varepsilon+\delta},\qquad
P_2\gtrsim d^{k\varepsilon+\delta},\qquad
P_3\gtrsim d^\delta,
$$

$\delta>0$ 为任意小的固定余量；宽度与步骤还须满足原定理要求。总瓶颈指数是 $\max(1+\varepsilon,k\varepsilon)$。预条件在证明中有作用，不能省略以后声称同一界适用于普通未预条件 SGD。对更深目标，论文给出扩展分析；联合训练的可见维数降低还通过数值实验检验，证据等级与浅模型的严格结果不同。

文中把这一过程类比粗粒化很有启发，但它已经有完整的概率与谱描述。若要进一步声称 RG，仍需指定被消去的自由度、保持的观测量与可迭代的尺度变换，参见[综合讨论](/notes/hierarchical_data/mechanism_comparison/#effective-theory)。

### 4.3 一个具体指数例子说明逐层收益与瓶颈 {#sight-exponents}

令内部多项式次数 $k=4$、中间维数 $r=d^{1/3}$。第 4.2 节的受控算法给出三个阶段的样本指数，忽略任意小余量后分别为 $4/3,4/3,0$。总样本的主导指数因而接近 $4/3$，而不是直接在 $d$ 维空间中搜索一般四次方向所提示的 $4$。这里的关键不是四次函数突然变成一次，而是先要找的子空间有 $d^{1/3}$ 维，第二阶段的四次操作也只在这个空间里执行。

若把 $\varepsilon$ 改为 $1/2$，两个主要指数成为 $3/2$ 和 $2$，总瓶颈改由第二阶段控制。若把 $\varepsilon$ 再减小，第二阶段更便宜，但第一阶段始终有约 $d^{1+\varepsilon}$ 的成本。这给出一个直观的交叉尺度：两项指数在 $1+\varepsilon=k\varepsilon$，即 $\varepsilon=1/(k-1)$ 处相等。这个交叉是**规定模型与算法的资源配平**，不意味着真实数据会通过学习自动选择最优潜在维数。

这个比较也说明“深度越多越好”没有直接从定理得出。每引入一个中间空间，就多一次需要估计的表示和可能放大的误差；如果下一层并未显著降低搜索维数或信息指数，分层的开销可以超过收益。MIGHT 的价值在于展示某些目标族确实允许递归减少有效搜索空间，不能据此推断任意组合函数的每个隐层都会产生同样的指数节省。

### 4.4 中间维数是统计对象，不是神经元个数的同义词 {#dimension-distinction}

SIGHT 的 $d^\varepsilon$ 是教师目标中的真实潜在子空间维数；学生网络的宽度是算法资源，二者并不自动相等。一个很宽的网络可以只学到低秩的任务相关子空间；一个狭窄瓶颈也可能过早丢失对预测必要的信息。要说学习降低了“有效维数”，最好同时提供潜在子空间重合、预测充分性或条件风险的证据，而不是只统计隐藏层用了多少个单元。

在数学上，若 $W_*$ 的行空间已经被恢复到小主角度，之后可以在该子空间坐标里重写目标；这个方向的成功不要求逐行找回 $W_*$ 的原始坐标标签。事实上 $W_*$ 经过允许的内部正交变换后，后续函数可同步改写，而最终输入输出函数不变。恢复哪一种对象、允许什么规范自由度，是判断“中间表示是否正确”的先决条件。

## 5. Tabanelli 2026：在多项式特征空间里逐层做谱恢复 {#hierarchical-spectral}

### 5.1 谱估计面对的是哪一个大空间 {#spectral-model}

Tabanelli、Dandi、Pesce 与 Krzakala 的 ICML 论文考虑密集随机张量生成的非线性中间特征。记 $\Phi_k(x)$ 为展平并归一化的 $k$ 阶 Hermite 特征，维数 $D_k\asymp d^k$，有

$$
h^{(1)}=A^{(1)}\Phi_k(x)\in\mathbb R^{d^\varepsilon},\qquad
h^{(2)}=\langle A^{(2)},H_2(h^{(1)})\rangle,\qquad y=g_*(h^{(2)}).
$$

不是先寻找少量原始线性坐标，而是从巨大的多项式空间恢复一个较小的**非线性子空间**。第一步计算标签加权二阶 Hermite 矩阵

$$
\widehat C_k^{(1)}=\frac1P\sum_{\mu=1}^P y_\mu H_2(\Phi_k(x_\mu)),
$$

然后做次数、秩选择和谱清理，提取 $A^{(1)}$ 的行空间；在恢复的坐标中再估计二阶形式，最后拟合一维 $g_*$。[原文 Algorithm 1、§3–4](https://proceedings.mlr.press/v306/tabanelli26a.html) [[4]](#ref-4)。

### 5.2 信号—噪声配平与旋转歧义 {#spectral-recovery}

其 Gaussian 等价分析在受控随机张量类中，将信号联系到 $A^{(1)\top}A^{(2)}A^{(1)}$。恢复中间维数 $d^\varepsilon$ 的主导样本尺度是 $P\gg d^{k+\varepsilon}$，第二步在理想中间坐标中的尺度为 $P\gg d^{2\varepsilon}$。前者主导总成本，反映多项式特征维数乘中间秩，而非最终目标的全次数。

需要保留三个限制。第一，证明依赖随机密集张量、可控收缩项与中间近 Gaussian 性，不是任意确定张量。第二，正文严格论证集中在 $\varepsilon<1/2$ 等条件，并要求相应估计器算子范数受控；更宽参数区间的表现还包含数值证据。第三，原文 Theorem 3.3 首先控制信号子空间上的投影，完整的谱恢复还利用谱清理与后续集中分析，不能将投影误差式直接冒充全空间算子误差界。

恢复的是子空间而非每个潜变量的唯一命名。若 $\widehat h^{(1)}\simeq Oh^{(1)}$，后层可同时改写二次型来吸收正交旋转，所以预测成功不需要 $O=I$。这与[共享表示的对齐问题](/notes/hierarchical_data/shared_representations/)一致。该算法是层级谱估计，虽可与梯度低阶信号建立联系，定理没有直接追踪端到端 SGD。

### 5.3 为什么标签加权二阶矩能显露隐藏方向 {#wick-example}

先暂时离开复杂随机张量，考虑一个可完全算出的 Gaussian 例子。令 $z\sim\mathcal N(0,I_D)$，$A=A^\top$，目标是一个中心化二次型

$$
y=z^\top Az-\operatorname{tr}A.
$$

用 Wick 恒等式 $\mathbb E[z_i z_j z_k z_l]=\delta_{ij}\delta_{kl}+\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}$，直接得到

$$
\mathbb E\big[y(zz^\top-I)\big]=2A.
$$

也就是说，**对标签加权的二阶 Hermite 矩求谱，恰好恢复目标二次型的方向**。例如 $D=2$，取 $A_{12}=A_{21}=a$、对角为零，则 $y=2az_1z_2$，估计矩的总体值为

$$
2A=\begin{pmatrix}0&2a\\2a&0\end{pmatrix}.
$$

两条正交特征方向 $(1,1)/\sqrt2$、$(1,-1)/\sqrt2$ 的特征值分别为 $2a,-2a$。如果只保留最大**正**特征值，可能丢掉同样有信息的负方向；根据目标和算法，必须明确是按正谱、负谱还是奇异值选择。具体论文的归一化 Hermite 定义会改变常数因子，但不改变这一机制。

Tabanelli 的层级算法在更大的 $\Phi_k(x)$ 空间使用类似的标签加权二阶对象。困难随之从“为什么二阶矩含信息”转为“在维数约 $d^k$、信号秩约 $d^\varepsilon$ 且存在随机张量收缩时，如何把总体谱与经验噪声分开”。后者才需要论文的 Gaussian 等价、谱清理和精细集中分析，不能被上述两行 Wick 计算替代。

### 5.4 哪些谱恢复误差会继续传到下一层 {#spectral-error-propagation}

设第一阶段得到的子空间投影为 $\widehat\Pi$，真实投影为 $\Pi_* $。若下一阶段只使用 $\widehat\Pi\Phi_k(x)$，与真实中间坐标的偏差可以由 $\|(\widehat\Pi-\Pi_*)\Phi_k(x)\|$ 衡量。它不只取决于投影误差的算子范数，还取决于输入在相应方向上的方差和高阶矩。即使两个投影在平均意义下接近，后续二次函数仍可能放大某些稀有但幅度较大的误差。

另一方面，预测器通常不需要恢复潜变量的固定命名。若 $\widehat h\simeq Oh$ 且 $O$ 正交，二次型 $\langle A,H_2(h)\rangle$ 可以通过同时将后层系数变换到相应坐标系而保持不变。因此合理的数值指标是子空间主角度、下游最优风险及新坐标下的规则可恢复性，而不是逐神经元均方距离。比较时还应保留一个独立验证集，避免先选谱方向、再用同一噪声矩报告过于乐观的恢复误差。

这类误差传递也解释了为何逐层算法喜欢样本分割：下一阶段看到的表示估计误差在证明中较易与新估计噪声分开。但数据分割有资源成本，严格证明并不自动说明这种资源安排对实际联合训练是最优的。若要取消分割，需要处理表示与下游估计在同一数据上的相关，而不只是把样本数替换成总预算。

## 6. 两种容易被误归类的补充结果 {#other-representations}

### 6.1 Chen 2020：表示有用，不代表表示由任务学得 {#chen-random-representation}

Chen 等人构造固定随机非线性表示 $g(x)=\sigma(Vx+b)$，用无标签样本估计协方差并白化为 $h(x)=\widehat\Sigma^{-1/2}g(x)$。然后在 $h$ 上训练宽两层网络的二阶 Taylor 模型，带适当正则化并寻找二阶稳定点。白化改善表示的各向同性，随机非线性表示提供可供二次模型组合的基函数。[原文 Algorithm 1、Theorems 1–3](https://papers.nips.cc/paper/2020/hash/fb647ca6672b0930e9d00dc384d8b16f-Abstract.html) [[5]](#ref-5)。

对指定低秩 $p$ 次多项式类，仅看维数依赖，文中 Quad-Neural 上界约为 $\widetilde O(d^{\lceil p/2\rceil})$，所比较的 Quad-Raw 已知上界约为 $\widetilde O(d^{p-1})$。这是特定模型与界的比较，不是对一切 raw-input 学习器的下界。论文还证明无限宽 NTK 的相应比较不能获得同样改进，说明表示必须和其后的学习模型一起判断。

它在本综述中说明“表示能改变可用组合与统计几何”，并未证明监督训练从数据中发现这个随机中间表示。它与 Nichani 的任务相关核投影、Wang 的中间 $p$ 学习，是相连但不同的贡献。

### 6.2 Daniely 2026：可观察的丰富标签提供中间学习入口 {#daniely-label-hierarchy}

Daniely 的 ICML 工作观察到，若只看一个最终标签，任意深组合可以非常难学；若同时观测许多粒度不同的标签，问题会改变。其训练样本是 $(x,f_1^*(x),\ldots,f_n^*(x))$，存在算法未知的集合

$$
L_1\subseteq L_2\subseteq\cdots\subseteq L_r=[n],
$$

$L_1$ 中标签是输入的有界复杂度多项式阈值函数，$L_i$ 中标签又是 $L_{i-1}$ 标签的简单函数，并有 margin/鲁棒性条件。逐层残差网络利用已预测的标签，为下一层创建容易拟合的输入，按层优化带正则化的凸子问题，获得固定局部复杂度下的多项式学习保证。[Definitions 3.1–3.3、Algorithm 1、Theorem 4.2](https://proceedings.mlr.press/v306/daniely26a.html) [[6]](#ref-6)。

这允许非常深的计算结构，而算法不需要知道哪一个标签属于哪一级。但中间相关标签本身是可观察监督；未知的是标签之间的层级。它不能被改写成“仅有一个最终标签的普通深网络可以发现任意深电路”。与 staircase 相比，前者把多个有效梯级放在一个标量目标的 Fourier 分量中，Daniely 则把它们放在同时观测的输出任务中。提供什么监督，本身就是学习机制的一部分。

### 6.3 同一个中间变量，三种获取方式的含义不同 {#access-to-intermediates}

设真实预测可写为 $y=g(z)$，$z=p(x)$。现在比较三种训练条件：

1. 只给 $(x,y)$，并要求算法通过最终标签找到 $z$。Nichani 与 Wang 的相关机制属于这一类，因而必须解释为什么最终标签对 $z$ 留下可提取的低阶投影。
2. 给 $(x,z,y)$ 或多个与 $z$ 对应的中间任务标签。此时训练可以直接监督中间表示，困难由“发现中间目标”转为“拟合多个已标注目标并正确利用它们”；Daniely 的框架属于更一般的可观察标签层级。
3. 预先提供固定随机非线性映射 $\phi(x)$，再证明一个后续算法可利用它。Chen 的结果说明有用基函数可以不经过同一任务的监督学习就存在，但这与恢复教师潜变量 $z$ 的结论不同。

这三种情形可能在最终测试误差上不可区分，却对应不同的学习机制。若只展示深网络能拟合 $y$，不能推断它采用第一种途径；它也可能依赖架构中预置的坐标变换、冗余随机特征或训练数据附带的额外任务。比较模型之前先核对可用信息，往往比争论“哪个架构更深”更重要。

## 7. 这些工作合在一起说明了什么 {#synthesis}

新的中间表示有三种不同价值：使目标的高频部分可由低频种子组合得到；把下一阶段搜索放进低维潜在空间；把原本高阶的输入统计改写成当前坐标中的低阶统计。它们可以同时发生，但维数下降不是自动可学性的充分条件。若外层响应对所需中间特征完全不敏感、低阶矩消失，或者误差经非线性放大，链条仍会中断。

这些证明大多通过逐层训练、样本分割、瓶颈或谱估计隔离机制，优点是能指出每一步到底提取什么，代价是减少了联合训练中的反馈与竞争。[下一篇](/notes/hierarchical_data/coupled_feature_learning/)将问：如果第一阶段表示有偏差，或者某个底层方向仅在上层响应存在时才显露，是否必须让各层共同修正？[RHM 专题](/notes/hierarchical_data/generative_learning/)则检验离散同义字典、局部结构和任务监督如何改变这些问题。

### 7.1 用两个独立检验避免把压缩当成机制 {#representation-tests}

对一个可控教师目标，至少可以把“学到中间表示”拆成两个实验。第一个是**任务充分性检验**：冻结网络中间层，只重新训练一个资源受限的读出，测量它在新样本上的风险是否接近使用真实 $p(x)$ 的读出。第二个是**几何或潜变量恢复检验**：对允许的旋转、缩放或非线性重参数化做对齐后，测量模型中间变量与真实子空间的关系。前者回答“够不够用”，后者回答“是否接近预期机制”，两者都不是另一个的自动推论。

再加一个因果干预：在保留输入分布和模型容量的前提下，打乱已学中间特征与下游样本的对应关系，观察预测是否显著退化。如果某个看似恢复良好的中间变量被打乱后，模型预测几乎不变，说明它可能只是一种可读出但未被实际使用的旁路信息。相反，若打乱后性能下降，也仍需排除简单的分布外破坏；应使用匹配边际分布、改变特定条件依赖的干预作为对照。

这些实验对 RHM 同样适用：子树父状态可被线性读出，不代表网络必须用它组合更高一级的规则；特定坐标的恢复也不必要求与生成器的内部命名完全一致。真正值得比较的是中间表示的任务充分性、形成成本，以及其在后续计算中是否承担了必要的因果作用。

## 8. 参考文献 {#references}

1. <span id="ref-1"></span>Nichani E, Damian A, Lee JD（2023）。[Provable Guarantees for Nonlinear Feature Learning in Three-Layer Neural Networks](https://papers.nips.cc/paper_files/paper/2023/hash/236b6a814a1d2c0ff504ca7bf380f7ff-Abstract-Conference.html)。*NeurIPS 36*。核投影作为非线性特征、一次步算法及特定深浅分离。
2. <span id="ref-2"></span>Wang Z, Nichani E, Lee JD（2024）。[Learning Hierarchical Polynomials with Three-Layer Neural Networks](https://proceedings.iclr.cc/paper_files/paper/2024/hash/af536fa4312b5b2e5c541156cd31e93b-Abstract-Conference.html)。*ICLR 2024*。近似 Stein 关系、瓶颈残差结构与两阶段训练。
3. <span id="ref-3"></span>Dandi Y, Pesce L, Zdeborová L, Krzakala F（2025）。[The Computational Advantage of Depth in Learning High-Dimensional Hierarchical Targets](https://papers.nips.cc/paper_files/paper/2025/hash/eb3878c1dcbfff9ee95d5d033e5f5942-Abstract-Conference.html)。*NeurIPS 38*。SIGHT/MIGHT、潜在维数降低及受控逐层训练。
4. <span id="ref-4"></span>Tabanelli H, Dandi Y, Pesce L, Krzakala F（2026）。[Deep Learning of Compositional Targets with Hierarchical Spectral Methods](https://proceedings.mlr.press/v306/tabanelli26a.html)。*ICML，PMLR 306:117916–117938*。Hermite 矩估计、随机张量与层级 Gaussian 等价。
5. <span id="ref-5"></span>Chen M, Bai Y, Lee JD, Zhao T, Wang H, Xiong C, Socher R（2020）。[Towards Understanding Hierarchical Learning: Benefits of Neural Representations](https://papers.nips.cc/paper/2020/hash/fb647ca6672b0930e9d00dc384d8b16f-Abstract.html)。*NeurIPS 33*。固定随机表示、白化与二阶 Taylor 学习模型。
6. <span id="ref-6"></span>Daniely A（2026）。[Deep Networks Learn Deep Hierarchical Models](https://proceedings.mlr.press/v306/daniely26a.html)。*ICML，PMLR 306:22909–22936*。多标签层级、残差架构与逐层学习；中间标签可观察。

---

[← 理论 02 · 梯度可达性](/notes/hierarchical_data/gradient_accessibility/) · [专题目录](/notes/hierarchical_data/#contents) · [下一篇：理论 04 · 跨层反馈与共同演化 →](/notes/hierarchical_data/coupled_feature_learning/)
