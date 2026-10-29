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
u(x)=u^\top\sigma_1(Vx+s),\qquad
f_\theta(x)=\nu(x)+c^\top\operatorname{ReLU}(a\nu(x)+b).
$$

$V,s,a,b$ 固定；先对 $u$ 多步梯度下降，以得到 $\nu\simeq c_0p$，再在独立数据上训练 $c$ 拟合外层。训练器没有 $p$ 的标签，第一阶段仍用最终标签 $y$，依靠组合投影的非退化性提取 $p$。与 Nichani 的一次步相比，反复更新能更有效累计中间信号。

对任意小的固定 $\alpha>0$，原定理以 $P\ge d^{k+3\alpha}$、第一模块宽度至少 $d^{k+\alpha}$ 等资源达到 $\widetilde O(d^{-\alpha})$ 风险。常见的“约 $d^k$”是概述其维数指数，不能读成在任意消失误差下严格只有 $d^k$ 乘对数。对 $k=2$，合适的迹零、Frobenius 范数归一化且算子范数 $O(d^{-1/2})$ 的二次型，将先前约 $d^4$ 的保证改进到接近 $d^2$ 的指数。

相应旋转对称核方法面对次数 $kq$ 的目标会承担约 $d^{kq}$ 的困难。这里的改进来自先学 $p$ 再拟合 $g$，同时强烈依赖目标的扩散结构、Gaussian 输入、瓶颈架构和逐层训练。它不是普通端到端 SGD 自动发现任意组合分解的定理。

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

## 7. 这些工作合在一起说明了什么 {#synthesis}

新的中间表示有三种不同价值：使目标的高频部分可由低频种子组合得到；把下一阶段搜索放进低维潜在空间；把原本高阶的输入统计改写成当前坐标中的低阶统计。它们可以同时发生，但维数下降不是自动可学性的充分条件。若外层响应对所需中间特征完全不敏感、低阶矩消失，或者误差经非线性放大，链条仍会中断。

这些证明大多通过逐层训练、样本分割、瓶颈或谱估计隔离机制，优点是能指出每一步到底提取什么，代价是减少了联合训练中的反馈与竞争。[下一篇](/notes/hierarchical_data/coupled_feature_learning/)将问：如果第一阶段表示有偏差，或者某个底层方向仅在上层响应存在时才显露，是否必须让各层共同修正？[RHM 专题](/notes/hierarchical_data/generative_learning/)则检验离散同义字典、局部结构和任务监督如何改变这些问题。

## 8. 参考文献 {#references}

1. <span id="ref-1"></span>Nichani E, Damian A, Lee JD（2023）。[Provable Guarantees for Nonlinear Feature Learning in Three-Layer Neural Networks](https://papers.nips.cc/paper_files/paper/2023/hash/236b6a814a1d2c0ff504ca7bf380f7ff-Abstract-Conference.html)。*NeurIPS 36*。核投影作为非线性特征、一次步算法及特定深浅分离。
2. <span id="ref-2"></span>Wang Z, Nichani E, Lee JD（2024）。[Learning Hierarchical Polynomials with Three-Layer Neural Networks](https://proceedings.iclr.cc/paper_files/paper/2024/hash/af536fa4312b5b2e5c541156cd31e93b-Abstract-Conference.html)。*ICLR 2024*。近似 Stein 关系、瓶颈残差结构与两阶段训练。
3. <span id="ref-3"></span>Dandi Y, Pesce L, Zdeborová L, Krzakala F（2025）。[The Computational Advantage of Depth in Learning High-Dimensional Hierarchical Targets](https://papers.nips.cc/paper_files/paper/2025/hash/eb3878c1dcbfff9ee95d5d033e5f5942-Abstract-Conference.html)。*NeurIPS 38*。SIGHT/MIGHT、潜在维数降低及受控逐层训练。
4. <span id="ref-4"></span>Tabanelli H, Dandi Y, Pesce L, Krzakala F（2026）。[Deep Learning of Compositional Targets with Hierarchical Spectral Methods](https://proceedings.mlr.press/v306/tabanelli26a.html)。*ICML，PMLR 306:117916–117938*。Hermite 矩估计、随机张量与层级 Gaussian 等价。
5. <span id="ref-5"></span>Chen M, Bai Y, Lee JD, Zhao T, Wang H, Xiong C, Socher R（2020）。[Towards Understanding Hierarchical Learning: Benefits of Neural Representations](https://papers.nips.cc/paper/2020/hash/fb647ca6672b0930e9d00dc384d8b16f-Abstract.html)。*NeurIPS 33*。固定随机表示、白化与二阶 Taylor 学习模型。
6. <span id="ref-6"></span>Daniely A（2026）。[Deep Networks Learn Deep Hierarchical Models](https://proceedings.mlr.press/v306/daniely26a.html)。*ICML，PMLR 306:22909–22936*。多标签层级、残差架构与逐层学习；中间标签可观察。

---

[← 理论 02 · 梯度可达性](/notes/hierarchical_data/gradient_accessibility/) · [专题目录](/notes/hierarchical_data/#contents) · [下一篇：理论 04 · 跨层反馈与共同演化 →](/notes/hierarchical_data/coupled_feature_learning/)
