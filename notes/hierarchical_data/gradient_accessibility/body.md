> **Hierarchical Learning · 理论 02**　[← 理论 01 · 谱动力学](/notes/hierarchical_data/spectral_dynamics/) · [总论](/notes/hierarchical_data/) · [理论 03 · 非线性中间表示 →](/notes/hierarchical_data/nonlinear_depth/)

**本篇问题：** 一个复杂目标明明可以由小网络表达，为什么梯度仍可能找不到它？什么结构能把学习变成一条可走通的路径？

Abbe 等人的 staircase 路线将“层级有助于学习”变成一个有关统计投影与支持集发现的问题。关键不是总多项式次数小，而是在发现部分坐标以后，剩余目标是否对新的坐标留下足够低阶的信号。2021 年提供专门深层算法的充分条件，2022 年研究两层 mean-field SGD 的近乎刻画，2023 年进一步分析有限网络的长时间在线 SGD。三篇的模型、算法与证明强度都不相同。

## 1. 用 Fourier 支持集表示学习的入口 {#fourier-signals}

### 1.1 Boolean 展开与未知坐标 {#boolean-expansion}

令 $x\sim\operatorname{Unif}\{\pm1\}^d$，$\chi_S(x)=\prod_{i\in S}x_i$，则

$$
f_*(x)=\sum_{S\subseteq[d]}\widehat f_*(S)\chi_S(x),\qquad
\widehat f_*(S)=\mathbb E[f_*(x)\chi_S(x)].
$$

$\chi_S$ 构成正交基。我们主要考虑只依赖少量未知坐标的目标：$f_*(x)=h_*(x_I)$，$|I|=r\ll d$。潜在维数小不等于目标容易发现，因为训练器还不知道 $I$ 在哪里。以下 $r$ 是潜在维数，样本量始终写 $P$；2022/2023 原文使用 $P$ 表示潜在维数时，已在这里换记号。

纯奇偶函数 $x_1\cdots x_k$ 对任意更低次数的 Fourier 特征都正交。一个低阶统计学习器不能从均值或二阶相关中获知其支持。这个障碍依赖输入分布、可用统计与计算模型；例如有专门算法可以用无噪 Boolean 方程求解某些奇偶问题，所以“纯奇偶对一切算法都困难”不是这里的结论。

### 1.2 加入低阶项为什么改变学习路线 {#conditional-signal}

比较

$$
f_{\rm parity}=x_1x_2x_3,\qquad
f_{\rm stair}=x_1+x_1x_2+x_1x_2x_3.
$$

两者最高次数相同，最初的低阶投影却不同。对第二个目标，$\mathbb E[yx_1]=1$。若已经学到 $z=x_1$，则

$$
\mathbb E[yzx_2]=1.
$$

发现 $x_2$ 以后，$\mathbb E[y(x_1x_2)x_3]=1$。这组等式直接来自 Fourier 正交性，是说明统计机制的算例，不是声称每个 SGD 神经元都逐字实现这些乘法。网络的非线性导数和已有权重可以产生类似的乘积投影，严格论文要证明这种投影确实在指定动力学中形成。

这里不一定产生新的标签信息：这些相关性已经包含在联合分布中。发生改变的是学习器当前能以低复杂度提取的统计。较简单表示把原本高阶的相关变成当前坐标中的低阶相关，提供新的梯度方向。

## 2. Staircase：每个高阶项都有可接续的前项 {#staircase}

### 2.1 结构条件到底要求什么 {#staircase-definition}

2021 年论文对稀疏 Fourier 多项式提出 $[1/M,M]$ staircase 条件：非零系数的绝对值在 $[1/M,M]$ 内；每个非零 $S$ 若 $|S|\ge2$，都存在非零 $S'\subset S$ 满足 $|S\setminus S'|=1$。因此每个高阶项沿着非零支持可以一直追溯到一次项。它允许多条链，既不是只许一个总链，也不是说所有低阶子集必须出现。[Abbe et al., 2021，Definition 1.1](https://papers.nips.cc/paper_files/paper/2021/hash/e2db7186375992e729165726762cb4c1-Abstract.html) [[1]](#ref-1)。

系数下界有机制意义：作为梯级的项如果小到淹没于采样噪声，形式上存在的链并不提供高效入口。系数上界与稀疏度则控制函数规模和算法资源。结构条件只是保证的一部分。

### 2.2 证明的是哪一种网络训练 {#staircase-algorithm}

论文构造稀疏随机连接的深网络，使用二次激活、输入到后续层的连接及特定正则化。训练按层进行，层内逐神经元运行随机块坐标下降，并包含阈值、扰动等算法细节。对 $s$ 个非零项的 staircase 多项式，该构造在 $d,s,M,1/\epsilon,1/\delta$ 的多项式时间和样本内，以至少 $1-\delta$ 的概率达到误差 $\epsilon$。

二次神经元可以将已知项与一个新坐标组合，随机连边则以足够概率把所需前项送到同一模块，训练与正则化避免大量不相关连接主导。这个证明把“可表达的组合”变成“能由指定局部优化步骤得到的组合”。它不是普通稠密 ReLU 网络端到端 SGD 的一般收敛定理。原文同时展示后者学习 staircase、难以学习纯奇偶的实验；实验与构造性证明应该并列阅读。

## 3. Merged staircase：已有坐标可以来自不同分支 {#merged-staircase}

### 3.1 支持集的合并比嵌套链更一般 {#msp-definition}

2022 年令非零 Fourier 支持组成 $\mathcal S=\{S_1,\ldots,S_m\}$。Merged-staircase property（MSP）要求存在排序，使

$$
\left|S_i\setminus\bigcup_{j<i}S_j\right|\le1
\quad\text{对每个 }i.
$$

注意并集取的是所有已经出现的支持，而不只是前一个支持。例如 $x_1+x_2+x_1x_2x_3$ 满足 MSP：前两项各打开一个坐标，第三项只带来 $x_3$。它不满足前一节严格的 staircase 条件，因为三次项缺少非零二次前项。$x_1+x_1x_2x_3$ 则不满足 MSP，最后一步同时引入两个未知坐标。

“合并”的意义是多个支路已学到的坐标可以共同作为后续信号的条件。它并不要求网络先显式构造所有前项的乘积，也没有规定一个非零单项式对应一层网络。

### 3.2 两层网络中为什么也会有层级学习 {#msp-dynamics}

论文研究两层网络

$$
f_\rho(x)=\int a\,\sigma(w^\top x+b)\,\rho(da,dw,db),
$$

及有限宽批量、单遍 SGD 在 mean-field 标度中的极限，潜在维数 $r$ 固定、环境维数 $d\to\infty$。它构造一个与 $d$ 无关的潜在参数分布动力学，并将 $O(d)$ 样本标度中的强 SGD 可学习性联系到该动力学能否从潜在方向权重为零的初态达到零风险。[Abbe et al., 2022，§1–3](https://proceedings.mlr.press/v178/abbe22a.html) [[2]](#ref-2)。

直观上，许多无关坐标形成各向同性背景；新潜在坐标能否从该背景中长出，取决于现有支持是否产生非零漂移。非 MSP 目标在这一极限中有无法打开的方向，风险停在正值。MSP 则给出逐次打开方向的候选路径。层级在这里发生于训练轨迹中的支持发现，一个隐层足以执行多次发现。

### 3.3 “必要且近乎充分”中的限定 {#msp-boundary}

必要性针对论文定义的强 $O(d)$ 标度与训练极限，不是说非 MSP 函数永远不能用更大样本、其他初始化或其他算法学习。充分性针对适当激活和一般位置的 Fourier 系数：固定支持结构后，排除测度为零的退化系数组合。某些对称性会使漂移精确抵消，即便支持满足 MSP，也不能直接保证学习。

因此支持图提供可能路径，系数与对称性决定路径上是否有力。判断一组数据是否“满足层级学习条件”，不能只画支持集的箭头，还要检查相应统计矩是否非零、激活的相关导数是否存在，以及采用的极限是否保留该信号。

## 4. Leap complexity：最难的一步需要跨过多少未知方向 {#leap-complexity}

### 4.1 Boolean 定义与 Gaussian 的次数重数 {#leap-definition}

对 Boolean 支持集定义

$$
\operatorname{Leap}(h_*)=
\min_\pi\max_i\left|S_{\pi(i)}\setminus\bigcup_{j<i}S_{\pi(j)}\right|.
$$

MSP 恰为 leap 至多 1 的结构条件。纯 $k$ 次奇偶的 leap 是 $k$；完整 staircase 的 leap 是 1；$x_1+x_1x_2x_3+x_2x_3x_4$ 的 leap 是 2。高次已经熟悉的变量不计作新方向。

对 Gaussian 输入，用正交 Hermite 多项式的多重指标 $\nu\in\mathbb N^r$ 替代 Boolean 支持。每次对**从未出现的坐标**计入该项的次数重数：

$$
\operatorname{Leap}(h_*)=
\min_\pi\max_i\sum_{k=1}^r\nu_{\pi(i),k}
\mathbf1\{\nu_{\pi(j),k}=0\text{ 对所有 }j<i\}.
$$

所以单变量 $\operatorname{He}_k(z_1)$ 的 leap 也是 $k$，并非 1；一旦 $z_1$ 已经出现，其后更高次幂不会再次按新方向计费。Gaussian 的坐标可以旋转，因而 leap 与选择的潜在基底有关；原文还讨论旋转不变版本。不要把任意一个展开的支持图当成数据分布的绝对不变量。[Abbe et al., 2023，Definition 1](https://proceedings.mlr.press/v195/abbe23a.html) [[3]](#ref-3)。

### 4.2 有限维初始化如何留下弱梯度 {#weak-gradient}

随机单位尺度权重在一个未知方向上的初始重合通常约为 $d^{-1/2}$。若需要同时接触 $k$ 个新方向，梯度的有用投影会包含多个小重合的乘积。论文给出的 Boolean 单项式例子中，支持内梯度相关量级为 $d^{-(k-1)/2}$，支持外约为 $d^{-(k+1)/2}$。支持内信号并非严格为零，而是随维数减弱。

这解释了两个不同表述为什么可以共存：在 2022 年特定大维极限的有效动力学中某方向完全关闭；在有限 $d$、更长时间的 SGD 中，它仍可能通过弱初始相关长出来。改变时间与维数极限的顺序，会改变可见的过程。严格为零的对称屏障、渐近弱信号、有限样本中的噪声淹没应分别检验。

## 5. 从结构可达性到 saddle-to-saddle 轨迹 {#saddle-dynamics}

### 5.1 三种复杂度不能共用一个“时间” {#online-complexity}

2023 年提出的一般猜想是：对固定潜在维数、适当一般位置的目标，规则初始化网络的在线 SGD 步数满足

$$
T_{\rm steps}\asymp\widetilde\Theta\!\left(d^{\max(\operatorname{Leap}(h_*)-1,1)}\right),
$$

其中略去了精度因子。每步使用一个新样本，因此这套在线协议中 $P=T_{\rm steps}$；若每步计算至少要处理 $d$ 个输入，总运算成本相应为 $\widetilde\Theta(d^{\max(\operatorname{Leap},2)})$。这三者的联系来自在线协议，不适用于固定训练集的多轮重复学习。

**这个面向一般目标和普通训练的刻画是猜想。** 论文的严格正面结果针对 Gaussian 输入上的代表性函数类、两层光滑网络、逐层训练和投影 SGD 等技术设定，有限宽、离散时间分析超过了连续 mean-field 流的时间范围。对一般 Boolean 目标与任意网络不能直接引用上述等式为定理。CSQ/noisy population-GD 的下界也不是针对全部在线 SGD 的无条件下界。

### 5.2 鞍点之间移动的具体含义 {#saddle-escape}

已经学到一个支持子集时，网络先拟合与它有关的目标分量，损失快速下降；尚未发现的方向重合仍小，使参数停留在近似部分拟合的鞍区。下一次弱漂移累积到足够大以后，新坐标进入可学习表示，损失再次下降。最大 leap 常常决定最慢的逃逸。

原文对其受控类分析了这种支持逐次发现，并特别说明：严格的 saddle-to-saddle 行为证明要求发现过程中的 leaps 按所需条件增长；leaps 可以降低的一般情形还没有同样完整的结论。不能把任意目标的损失平台都解释为经过一串已经识别的鞍点。要支持该解释，至少要测量潜在方向重合、相应 Fourier/Hermite 分量和局部梯度，而不是只显示训练误差。

## 6. 如何与其他理论和 RHM 比较 {#comparison}

### 6.1 与谱速率的区别 {#spectrum-vs-accessibility}

[Saxe 模态](/notes/hierarchical_data/spectral_dynamics/#feature-times)从初始时刻就由固定 $s_\alpha$ 驱动。Staircase 的下一步则依赖已发现支持改变有效漂移。两者都能产生快慢曲线，但干预预测不同：若冻结前一步特征，后一个 Saxe 独立模态仍可增长；若消去 staircase 所需入口，相关路径可能关闭或变得极慢。

“促进”也未必意味着严格门控。一个已有方向可以把 $d$ 的困难指数降低，却未必把原来的非零梯度变成零与非零的二值开关。[跨层反馈](/notes/hierarchical_data/coupled_feature_learning/)还会让高层敏感度参与选择，超出单向支持合并。

### 6.2 RHM 有局部字典，却没有天然的 Fourier 梯级 {#rhm-link}

RHM 的随机语法、有限词表、同义替换与局部卷积结构，和这里的独立 Boolean 输入、未知稀疏坐标并不相同。[RHM 的标签条件统计](/notes/hierarchical_data/generative_learning/#sample-complexity)会随生成层数衰减，学习器也需要把局部元组归为父状态；这在机制上同样涉及“前一级表示让后一级任务可用”。但将一个语法树节点直接认作一个 Fourier 支持，或将 RHM 的样本阈值写成 leap 时间公式，都缺少依据。

已有[奇偶组合模型笔记](/notes/neural_network_theory/hierarchical_learning/)计算预先指定乘法架构中的链式梯度：它适合看弱幅度与长等待，却没有解决未知支持或配对的发现。本篇补上的正是这一步。架构已经知道正确配对、数据提供有效 Fourier 梯级、算法证明能发现支持，是三种不同的知识来源。

## 7. 参考文献 {#references}

1. <span id="ref-1"></span>Abbe E, Boix-Adserà E, Brennan MS, Bresler G, Nagaraj D（2021）。[The Staircase Property: How Hierarchical Structure Can Guide Deep Learning](https://papers.nips.cc/paper_files/paper/2021/hash/e2db7186375992e729165726762cb4c1-Abstract.html)。*NeurIPS 34*。正文 §1–2：支持链与逐层随机坐标训练；常规 ResNet 的观察为另列数值证据。
2. <span id="ref-2"></span>Abbe E, Boix-Adserà E, Misiakiewicz T（2022）。[The merged-staircase property: a necessary and nearly sufficient condition for SGD learning of sparse functions on two-layer neural networks](https://proceedings.mlr.press/v178/abbe22a.html)。*COLT，PMLR 178:4782–4887*。MSP、无维度潜在动力学及一般位置限制。
3. <span id="ref-3"></span>Abbe E, Boix-Adserà E, Misiakiewicz T（2023）。[SGD learning on neural networks: leap complexity and saddle-to-saddle dynamics](https://proceedings.mlr.press/v195/abbe23a.html)。*COLT，PMLR 195:2552–2623*。正文 Definition 1、Conjecture 2、§3 与证明：在线步数、受控函数类及鞍区逃逸；一般结论保留为猜想。

---

[← 理论 01 · 谱动力学](/notes/hierarchical_data/spectral_dynamics/) · [专题目录](/notes/hierarchical_data/#contents) · [下一篇：理论 03 · 非线性中间表示与深度优势 →](/notes/hierarchical_data/nonlinear_depth/)
