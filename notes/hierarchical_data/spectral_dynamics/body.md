> **Hierarchical Learning · 理论 01**　[总论](/notes/hierarchical_data/) · [理论 02 · 梯度可达性 →](/notes/hierarchical_data/gradient_accessibility/)

**本篇问题：** 数据中的相关结构，怎样变成内部表示的集体模态，又为什么在不同训练时刻显现？

Saxe、McClelland 与 Ganguli 的工作给出一个重要起点：即使网络表示的输入输出函数完全线性，因子化参数的学习仍然是非线性动力学。它能产生长平台、快速增长、语义分化以及模态之间的合作与竞争。但这里的“出现”首先是已有数据模态被放大，不等于网络发现了一个原先无法接触的非线性组合。后者要在[梯度可达性](/notes/hierarchical_data/gradient_accessibility/)和[非线性中间表示](/notes/hierarchical_data/nonlinear_depth/)中另行研究。

## 1. 一个线性函数为什么具有非线性学习过程 {#linear-map-nonlinear-flow}

### 1.1 数据、目标与时间的约定 {#setup}

考虑输入 $x\in\mathbb R^d$、多维目标 $y\in\mathbb R^q$ 和一个隐层的线性网络

$$
\widehat y=W_2W_1x,\qquad
\mathcal L=\frac12\mathbb E\|y-W_2W_1x\|^2.
$$

以下把学习率与样本平均的归一化吸收到连续时间 $t$；不同论文中的 $\tau$ 可以通过重标度恢复。假定输入已经白化，$\mathbb E[xx^\top]=I$，记 $C=\mathbb E[yx^\top]$。期望既可以指总体，也可以指满足白化条件的固定训练集平均；这两种选择的泛化含义不同。梯度流为

$$
\dot W_1=W_2^\top(C-W_2W_1),\qquad
\dot W_2=(C-W_2W_1)W_1^\top.
$$

方程在权重中包含三次项。直接训练一个矩阵 $M$ 时，$\dot M=C-M$；训练两个因子的乘积时，每个因子的更新都乘着另一个因子。小初始化于是同时提供低初始预测和小初始速度。这是平台期的来源之一。[Saxe et al., 2014，§1](https://arxiv.org/abs/1312.6120) [[1]](#ref-1)。

### 1.2 数据谱把高维矩阵变成少量协作变量 {#mode-reduction}

写 $C=U\operatorname{diag}(s_\alpha)V^\top$，$s_\alpha\geq0$。在与奇异向量对齐、不同模态占据正交隐方向的子流形上，可写

$$
W_1=R\operatorname{diag}(a_\alpha)V^\top,\qquad
W_2=U\operatorname{diag}(b_\alpha)R^\top,
\quad R^\top R=I.
$$

每个模态满足

$$
\dot a_\alpha=b_\alpha(s_\alpha-a_\alpha b_\alpha),\qquad
\dot b_\alpha=a_\alpha(s_\alpha-a_\alpha b_\alpha).
$$

数据提供的 $s_\alpha$ 是驱动，$a_\alpha b_\alpha$ 是已经解释的信号，二者的差是残差。两侧权重协作：输入侧把信息送入隐层，输出侧把同一信息读出。任一侧过小都会限制另一侧增长。$a_\alpha^2-b_\alpha^2$ 沿轨道守恒，所以“平衡”也是一个明确的初始化条件，而不是随时都可施加的近似。

这种约化在对齐初始化下精确；对一般随机初始化，初期还包含模态旋转与相互作用。小随机权重经常趋向接近该约化的轨道，但不能因此把标量方程称为任意矩阵初始化下的精确解。输入未白化时，输入协方差也进入动力学，单独知道 $C$ 的奇异值不再充分。

## 2. 平台、快速转变与特征时间 {#feature-times}

### 2.1 两个权重矩阵的精确平衡解 {#balanced-solution}

取 $a_\alpha=b_\alpha>0$，定义端到端幅度 $m_\alpha=a_\alpha b_\alpha$。于是

$$
\dot m_\alpha=2m_\alpha(s_\alpha-m_\alpha),\qquad
m_\alpha(t)=\frac{s_\alpha}
{1+\left(s_\alpha/m_\alpha(0)-1\right)e^{-2s_\alpha t}}.
$$

早期 $m_\alpha\ll s_\alpha$ 时，$m_\alpha\simeq m_\alpha(0)e^{2s_\alpha t}$；靠近目标后由残差项饱和。达到 $m_\alpha=\theta s_\alpha$ 的时间是

$$
t_{\alpha,\theta}=\frac{1}{2s_\alpha}
\log\!\left[\frac{\theta}{1-\theta}
\left(\frac{s_\alpha}{m_\alpha(0)}-1\right)\right],\quad 0<\theta<1.
$$

只要初始化足够小，主导行为是 $s_\alpha^{-1}\log[s_\alpha/m_\alpha(0)]$。较大奇异值通常更早形成；从 10% 到 90% 的转变宽度是 $\log9/s_\alpha$。初始化趋于零时，等待相对于转变宽度变长，曲线看起来越来越突然。有限初始化、有限系统中的这条平滑曲线，并不构成热力学奇点。

这里每个 $s_\alpha>0$ 的模态从一开始就有增长率。小模态尚未可见，不表示它必须等待大模态完成以后才得到梯度。判别这一点的方法是检查早期增长率，而不是只看损失曲线中的台阶。

### 2.2 深度改变早期幂次，而不只有一个学习率因子 {#depth-initialization}

对 $D$ 个权重矩阵，在奇异方向对齐、每层正幅度相等时，$m_\alpha=\prod_{j=1}^D a_{j\alpha}$，有

$$
\dot a_{j\alpha}=(s_\alpha-m_\alpha)\prod_{k\ne j}a_{k\alpha},\qquad
\dot m_\alpha=Dm_\alpha^{\,2-2/D}(s_\alpha-m_\alpha).
$$

$D=1$ 是直接线性拟合，$D=2$ 是指数起步，$D>2$ 的早期增长还被额外的小幅度幂次抑制。令每层初始幅度为 $\varepsilon$，则 $m_0=\varepsilon^D$；早期积分给出 $D>2$ 时的等待量级 $\varepsilon^{-(D-2)}/[D(D-2)s_\alpha]$。这是固定归一化、固定平衡分支上的估计，不能与按端到端幅度固定的深度比较混用。

因此“深度总会减慢训练”也不成立。Saxe 2014 同时研究预训练和正交初始化，说明保留信号传播与梯度传播的初始化可以显著改变深度依赖。比较时应先说明固定的是每层方差、端到端奇异值，还是 Jacobian 的条件数。已有[深度优化笔记](/notes/neural_network_theory/deep_optimization/)讨论了这些量与梯度传播的关系。

## 3. 合作与竞争从哪里来 {#cooperation-competition}

### 3.1 对齐前的完整模态相互作用 {#full-modes}

把输入模态 $\alpha$ 到隐层的向量记为 $\boldsymbol a^\alpha$，隐层到输出模态 $\alpha$ 的向量记为 $\boldsymbol b^\alpha$。白化后的精确方程包含

$$
\dot{\boldsymbol a}^{\alpha}
=(s_\alpha-\boldsymbol a^\alpha\!\cdot\!\boldsymbol b^\alpha)\boldsymbol b^\alpha
-\sum_{\beta\ne\alpha}(\boldsymbol a^\alpha\!\cdot\!\boldsymbol b^\beta)\boldsymbol b^\beta,
$$

以及交换 $a,b$ 的另一式。相应的损失除去与参数无关的常数，是

$$
\mathcal L=\frac12\sum_\alpha(s_\alpha-\boldsymbol a^\alpha\!\cdot\!\boldsymbol b^\alpha)^2
+\frac12\sum_{\alpha\ne\beta}(\boldsymbol a^\alpha\!\cdot\!\boldsymbol b^\beta)^2.
$$

第一项使同一模态两侧配对；第二项惩罚把一个输入模态送到错误输出的串扰，推动不同模态正交化。所谓竞争在这里有具体对象：隐层方向与输出解释权。隐层维数不足时，稳定解保留最强的可容纳奇异模态。这不等价于两个无关技能在任意网络中必然互相抑制。

### 3.2 解耦本身消去了什么 {#decoupling-limit}

标量 logistic 方程已经把串扰降到零。它适合计算各独立模态的特征时间，却不能再用来证明非平凡的模态竞争。若希望研究模态相互促进或排斥，必须保留交叉重合、容量限制，或者引入共享路径。把所有耦合都约掉以后，再根据几条错开的 S 曲线推断“相互作用”，会把建模假设误当作理论发现。

## 4. 从层级语义到奇异模态 {#semantic-hierarchy}

### 4.1 树结构怎样进入数据矩阵 {#tree-spectrum}

Saxe 2019 研究项目到属性的映射：输入可以是项目的 one-hot 编码，输出是该项目的多个语义属性。树上的分支生成相关属性，近亲项目共享较多属性。任务是从项目—属性经验中学到预测和隐层表示，不是从连续图像中自动识别部件。层级结构先进入项目之间的相关矩阵，再通过奇异值分解进入学习方程。[Saxe et al., 2019](https://www.pnas.org/doi/10.1073/pnas.1820226116) [[2]](#ref-2)。

一个可直接计算的说明性例子是四个项目，前两个属于一类，后两个属于另一类。设归一化后的属性相关矩阵为

$$
Q=\begin{pmatrix}
c_0&c_1&c_2&c_2\\c_1&c_0&c_2&c_2\\
c_2&c_2&c_0&c_1\\c_2&c_2&c_1&c_0
\end{pmatrix}.
$$

其均匀模态 $(1,1,1,1)/2$ 的特征值为 $c_0+c_1+2c_2$；类别对比 $(1,1,-1,-1)/2$ 为 $c_0+c_1-2c_2$；两个类内对比分别为 $c_0-c_1$。当 $c_0\ge c_1\ge c_2\ge0$ 时，粗类别对比不弱于细对比。$Q$ 与数据矩阵的 Gram 矩阵成比例，所以这些特征值对应奇异值的平方，比例由样本归一化决定。

这个例子解释的是“为什么粗模式可能更早”，并没有证明“任何分类树都从粗到细学”。若一个细分类别具有极强、稳定的属性信号，而某个上层分类统计很弱，谱排序可以改变。真正控制速率的是任务相关谱，树的层号只是通过生成参数影响该谱。

### 4.2 隐层几何、语义错觉与典型性 {#semantic-geometry}

在平衡解上 $W_1=R\operatorname{diag}(\sqrt{m_\alpha})V^\top$，因此项目的隐层 Gram 矩阵是

$$
H^\top H=V\operatorname{diag}(m_\alpha)V^\top.
$$

早期只保留强模态，项目主要按粗类别分组；细模态增长后，同一类别内部再分化。输出预测却是 $U\operatorname{diag}(m_\alpha)V^\top$，一个具体项目的属性由多个模态相加。不同模态在最终答案中可能相互抵消；只学到一部分时，预测会暂时偏向错误属性，形成论文讨论的语义错觉。它可以出现在总体损失持续下降的同时，说明单项预测不必单调改善。

“典型性”对应项目在类别相关奇异向量上的投影，而非事先附加的心理量。模型也分析了类别一致性、归纳投射与不同网络隐表示的相似关系。隐藏坐标可以相差正交变换，所以相似几何不要求逐神经元对应。这与[跨语言共享几何](/notes/hierarchical_data/shared_representations/#允许一次对齐后应该测量什么)中的坐标问题直接相连。

## 5. Neural Race：共享路径为什么可能更快 {#neural-race}

### 5.1 门控网络把架构变成路径统计 {#gated-networks}

Saxe、Sodhani 与 Lewallen 2022 将深线性约化扩展到门控网络。每个样本有规定的路径开关；门控后活动沿所选子网络线性传播，因此能用路径条件下的输入与输入输出协方差描述梯度。可解约化要求这些路径统计具有兼容的模态基底和相应对称性。门控使整体输入输出行为比单个线性映射丰富，但论文并未一般性证明网络会从零发现正确门控。[原文 §2–4](https://proceedings.mlr.press/v162/saxe22a.html) [[3]](#ref-3)。

在多输入域、多输出域的路由任务中，一个共享隐层会收到许多训练路径的相容更新。独占路径每次只能利用部分任务，强共享路径则积累更多一致证据。这里的“快”来自架构与数据共同决定的梯度复用，不能只由一个原始数据奇异值解释。

### 5.2 路径之间的竞争与抽象的边界 {#pathway-race}

在原文全路由、对称化的三矩阵约化中，$M$ 为域数，$q$ 为流经某共享模块的路径数，$B_1$ 表示相等的输入和输出模态幅度，$B_2$ 为共享中间幅度，$s,d_0$ 为相应对角统计。对单个模态，原文式 (16)–(17) 可写成

$$
\tau\dot B_1=\frac{\sqrt q}{M^2}B_2B_1(s-d_0B_2B_1^2),\qquad
\tau\dot B_2=\frac{q}{M^2}B_1^2(s-d_0B_2B_1^2).
$$

这里将原文的路径计数 $P$ 改为 $q$，避免与本专题样本数混淆。更多共享路径提升了两个更新中的系数，尤其提升共享层的增长。多种解法共存并解释同一残差时，较快的共享解可能先消耗可用误差，使独占解难以继续增长。论文的约化给出速率机制；随机初始化的完整网络与约化的吻合、未训练路径上的迁移还包含数值证据。不能把其模拟中的比例阈值当成通用常数。

抽象在这个任务中有可操作定义：不同输入域的同一对象被映射到相容表示，以便输出到未见过的域组合。若任务彼此冲突、门控不相容或初始化过大，共享不一定最快，更不保证迁移。它与 RHM 的权重共享可以相互启发，但前者共享的是路由表示，后者还包含空间位置复用与同义元组归并，不能直接套同一个计数公式。

## 6. 这一理论究竟解释了哪种“分阶段” {#mechanism-boundary}

本篇的基本机制是：数据已经定义了多个可增长模态，学习按不同速率显露它们。深度与初始化决定放大规律，路径复用进一步改变速率和残差分配。它没有一般解决非线性隐特征的发现，也没有把数据层级与网络层数一一对应。Saxe 语义模型里，一个隐层已经可以逐次形成许多语义分辨尺度。

与 RHM 比较时尤其要固定横轴。[语言统计](/notes/hierarchical_data/language_statistics/)中的不同样本阈值，决定经验统计何时可辨；本篇的 $t_{\alpha,\theta}$ 则决定固定算子下模态何时增长。有限数据也会污染经验奇异值与奇异向量，进而改变训练顺序；但要从这种联系得到联合的 $P,t$ 预测，必须另外分析谱估计误差。不能将一个随 $P$ 的平台换名为随 $t$ 的平台。

下一篇将保留“分阶段”现象，却更换其来源：有些新坐标起初在有效梯度中不可见，只有已经发现的坐标改变统计投影后，学习路径才打开。

## 7. 参考文献 {#references}

1. <span id="ref-1"></span>Saxe AM, McClelland JL, Ganguli S（2014）。[Exact solutions to the nonlinear dynamics of learning in deep linear neural networks](https://arxiv.org/abs/1312.6120)。*ICLR 2014；arXiv:1312.6120*。本文使用其白化、模态约化、合作竞争与初始化分析。
2. <span id="ref-2"></span>Saxe AM, McClelland JL, Ganguli S（2019）。[A mathematical theory of semantic development in deep neural networks](https://doi.org/10.1073/pnas.1820226116)。*PNAS 116(23):11537–11546*；[作者预印本](https://arxiv.org/abs/1810.10531)。树生成语义、谱分化、典型性及表示几何。
3. <span id="ref-3"></span>Saxe AM, Sodhani S, Lewallen SJ（2022）。[The Neural Race Reduction: Dynamics of Abstraction in Gated Networks](https://proceedings.mlr.press/v162/saxe22a.html)。*ICML，PMLR 162:19287–19309*。路径条件统计、共享表示及路径竞争。

---

[← 总论](/notes/hierarchical_data/) · [专题目录](/notes/hierarchical_data/#contents) · [下一篇：理论 02 · 梯度可达性与学习瓶颈 →](/notes/hierarchical_data/gradient_accessibility/)
