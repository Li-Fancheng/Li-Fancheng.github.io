> **Hierarchical Learning · 理论 04**　[← 理论 03 · 非线性中间表示](/notes/hierarchical_data/nonlinear_depth/) · [总论](/notes/hierarchical_data/) · [RHM 01 · 层次生成 →](/notes/hierarchical_data/generative_learning/)

**本篇问题：** 上层学到的响应，怎样反过来改变底层应该学习的特征？能否把这种共同演化写成可分析的方程？

逐层构造把“先得到一个表示，再解下一阶段问题”隔离得很清楚，却也冻结了已经形成的表示。联合训练还包含另一条通路：上层决定最终预测对某个神经元、某个输入有多敏感，这个响应参与底层梯度；底层改变后又改变上层输入与响应。本篇先用 Neural LoFi 说明前向谱选择，再研究 backward feature correction 与 BackLoFi 对反馈的两种描述。2026 年两篇 LoFi 工作均按截至 2026-10-29 可核查的 arXiv 预印本版本讨论，不能将其特殊替代模型结论提升为一般训练定律。

## 1. Neural LoFi：低阶是相对于当前表示而言 {#forward-lofi}

### 1.1 从层内梯度展开得到相关算子 {#local-expansion}

令当前表示为 $z_{\ell-1}(x)$，一层的预激活为 $h_\ell=W_\ell z_{\ell-1}$。在逐层训练、各层初始化尺度分离、待学权重很小的条件下，Dandi 等人展开初期梯度，得到单个输入权重向量的形式

$$
w_i(t+1)-w_i(t)
\simeq\eta c_{0i}m_\ell+\eta c_{1i}C_\ell w_i(t),
$$

$$
m_\ell=\widehat{\mathbb E}[y z_{\ell-1}],\qquad
C_\ell=\widehat{\mathbb E}[y z_{\ell-1}z_{\ell-1}^\top].
$$

$c_{0i},c_{1i}$ 与激活导数、上层初始权重有关，正负号也重要；省去项受初始化阶数和局部时间范围控制。第一项沿标签相关均值漂移，第二项在标签加权二阶算子的特征方向上放大或缩小。在线性化范围内，若 $C_\ell v_j=\lambda_jv_j$，对应幅度含 $e^{c_{1i}\lambda_jt}$。[Neural LoFi v2，Proposition 1](https://arxiv.org/abs/2605.13612v2) [[1]](#ref-1)。

$C_\ell$ 不是普通无监督协方差，可以不半正定。算法会按采用的正负分支保留相关谱方向；讨论按绝对值的谱选择时，也要容许翻转标签或对应的增长分支。不能将它替换成“寻找输入方差最大的主成分”。

### 1.2 从局部方程到算法，中间还有一步设计 {#filter-lift}

Neural LoFi 将这种放大离散化为可计算原语：选取 $m_\ell$ 的方向及 $C_\ell$ 的若干相关特征向量，组成 $V_\ell$；得到滤波坐标 $q_\ell=V_\ell^\top z_{\ell-1}$；再用独立 Gaussian 随机矩阵做非线性提升：

$$
z_\ell(x)=\frac1{\sqrt{p_\ell}}\sigma(R_\ell q_\ell(x)),
$$

最后训练读出。硬选前几个谱方向、选择保留秩以及重新随机提升，是谱算法的设计，不能从局部短时 GD 展开直接推出任意时刻网络都等价于这一算法。

滤波和提升各有作用。滤波减少无关方向和估计噪声；非线性提升使已有坐标的乘积、曲率等进入下一层的简单函数空间。下一层仍然使用一次或二次相关，但这些相关作为原始输入的函数可以是高阶的。若只投影而不非线性提升，反复线性投影无法产生一般非线性组合表示。反之，若先丢掉了之后仍必需的方向，随机提升也不能恢复已经丢失的信息。

### 1.3 核几何是由标签逐层改写的 {#adaptive-kernel}

记 $K_{\ell-1}(x,x')=\langle z_{\ell-1}(x),z_{\ell-1}(x')\rangle$。在其诱导 RKHS 中，一个单位复杂度特征 $\phi=\langle v,z_{\ell-1}\rangle$ 满足

$$
v^\top C_\ell v=\widehat{\mathbb E}[y\phi(x)^2].
$$

所以谱选择是在当前几何中寻找“平方激活与标签最相关”的简单特征，而不要求 $\mathbb E[y\phi]\ne0$。例如对称偶函数的有用线性坐标可能没有一次相关，却有非零二次相关。正负分支可用相应 Rayleigh 商或其绝对值表述；后续方向还须正交化。

在无限提升宽度和谱间隙条件下，算法可写为核递归

$$
K_\ell(x,x')=\mathbb E_r[\sigma(r^\top q_\ell(x))\sigma(r^\top q_\ell(x'))].
$$

它区别于固定核方法之处在于 $q_\ell$ 已由标签选择，因而函数空间的简单方向随层数变化。它仍然是逐层、单向的构造；后层读出拟合得再好，也不会自动回到第一层修改已经选定的子空间。

## 2. 表示涌现的统计阈值与算法盲点 {#lofi-emergence}

### 2.1 应比较的是算子噪声与谱间隙 {#spectral-noise}

设总体相关算子为 $C_\ell^*$，经验误差 $\delta C_\ell=\widehat C_\ell-C_\ell^*$。对当前候选函数类，相关估计的均匀波动是

$$
\tau_\ell(P)=\sup_{\|\phi\|_{\mathcal H_{\ell-1}}=1}
\left|\widehat{\mathbb E}[y\phi^2]-\mathbb E[y\phi^2]\right|.
$$

有限维时它对应算子范数误差。若有用信号强度淹没于该波动，经验首模态容易由噪声决定；要稳定恢复一个指定子空间，还需其**谱间隙**超过误差。在论文的局部复杂度分析中，波动可受当前 RKHS 的有效维数控制，结构为“局部尺度乘 $\sqrt{D_{\rm eff}/P}$”，并带技术条件与对数因子。前一级表示改变候选空间、信号大小和有效维数，于是改变后一级的统计可见性。

这里随 $P$ 的谱分离首先是样本阈值。局部梯度展开又给出短时增长率，但这两件事并不提供一般端到端网络的完整训练时间预测。要得到固定数据上的特征出现时间，还须控制初始化投影、谱随轨迹变化、残差、离散步长和噪声。

### 2.2 当前低阶统计看不到什么 {#forward-blindness}

若有用特征对 $m_\ell,C_\ell$ 都不可见，原始 LoFi 可能停止；若错误选择导致不可逆信息丢失，更多前向层也不能保证修复。高阶统计可扩大范围，但会面对张量估计与计算困难。模型因此既提供机制，也明确展示低阶谱原语的能力边界。

这种边界与 [staircase](/notes/hierarchical_data/gradient_accessibility/)有关：低阶性依赖当前坐标，已有表示可以打开下一步。区别是 LoFi 明确规定用哪个算子、怎么选方向与提升；staircase/leap 描述的是特定神经训练协议中支持如何被发现，不能直接由相同的“逐步”措辞判定两者定理等价。

## 3. Backward Feature Correction：上层拟合怎样消除下层偏差 {#backward-feature-correction}

### 3.1 前向近似为什么可能留下有结构的误差 {#bfc-mechanism}

Allen-Zhu 与 Li 研究二次激活的层级目标与带跳连的深学生。目标包含不同层贡献，较深项的系数较小，形成 information gap：$\alpha_{\ell+1}\ll\alpha_\ell$。先学浅项可以给高层提供初步表示，但较小的深项也会污染低层拟合。若把低层永久冻结，该误差会成为后续问题的一部分。[COLT 2023 扩展摘要及完整 v6，§1.1–1.2](https://proceedings.mlr.press/v195/allen-zhu23a.html) [[2]](#ref-2)。

原文的说明性目标含 $x_1^4+x_2^4$，再加幅度 $\alpha$ 的更深组合。只拟合前一层时，可得到带 $O(\alpha)$ 混入的特征；上层先用这些近似特征拟合深项，带到最终预测中的误差便可能缩到 $O(\alpha^2)$。联合更新以后，较好解释的深项从残差中被扣除，低层不必再用错误方向代偿它，从而也能改进。如此往复形成论文描述的 $\alpha\to\alpha^2\to\alpha^3$ 修正链。

这里不是“梯度必须等下一层成熟才第一次非零”，而是**下层应拟合的残差与误差容忍度随上层改变**。前向表示提供高层入口，高层再改善低层估计，构成合作循环。

### 3.2 证明的范围与不能省略的训练条件 {#bfc-guarantee}

完整论文对规定的层级二次目标、宽度和条件数、信息间隙、受控输入分布，以及特定正则化和阶段安排的 SGD 变体，证明多项式时间和样本学习保证。深度可以随维数缓慢增长，主结果要求 $L=o(\log\log d)$ 等条件。论文还给某些固定特征/核及特定浅层模型的超多项式下界；对所有逐层算法的无条件下界并没有由“未找到高效替代算法”自动得到。

因此 BFC 提供了联合训练能反复修正表示的严格实例。它没有证明 information gap 在所有现实数据中存在，也没有证明移除 gap 后相同训练仍有保证。COLT 页面只收录一页扩展摘要；这里的机制与假设核对使用[完整 arXiv v6](https://arxiv.org/abs/2001.04413v6)，而非仅据摘要转述。

## 4. BackLoFi：在真实梯度中隔离相邻层的配对模态 {#backlofi-local}

### 4.1 精确梯度先给出什么 {#exact-gradient}

取归一化网络 $z_\ell=\sigma(h_\ell)/\sqrt{p_\ell}$，$h_\ell=W_\ell z_{\ell-1}$，标量预测 $f$。定义 $g_\ell(x)=\partial f/\partial h_\ell(x)$，$\Delta=y-f$。平方损失梯度流精确满足

$$
\dot W_\ell=\widehat{\mathbb E}[\Delta\,g_\ell z_{\ell-1}^\top],\quad
 g_\ell=\frac{\sigma'(h_\ell)}{\sqrt{p_\ell}}
\odot W_{\ell+1}^\top g_{\ell+1}.
$$

输入相关的导数与上层响应共同决定一个底层特征如何被使用。$g$ 本身是训练得到的函数，不是固定标签变换。[BackLoFi v1，§2](https://arxiv.org/abs/2610.32849v1) [[3]](#ref-3)。

### 4.2 小相邻权重、固定周围网络的受控结果 {#paired-growth}

对第 $\ell$ 层一个神经元，记输入行为 $w$，下一层读出它的列为 $\omega$。冻结周围的 $z=z_{\ell-1}$ 与 $g=g_{\ell+1}$，令这对权重很小，并假设 $\sigma(0)=0,\sigma'(0)>0$。一阶近似给出

$$
\dot w\simeq\kappa\widetilde A_\ell^\top\omega,\qquad
\dot\omega\simeq\kappa\widetilde A_\ell w,\qquad
\widetilde A_\ell=\widehat{\mathbb E}[\Delta\,g z^\top],\quad
\kappa=\frac{\sigma'(0)}{\sqrt{p_\ell}}.
$$

足够接近零预测初始化时才能进一步将 $\Delta$ 近似为 $y$，得到 $A_\ell=\widehat{\mathbb E}[y g z^\top]$。若 $A_\ell v_i=s_i u_i$，对应两侧幅度满足

$$
\dot a_i=\kappa s_i b_i,\quad\dot b_i=\kappa s_i a_i,
\qquad
\frac d{dt}(a_i\pm b_i)=\pm\kappa s_i(a_i\pm b_i).
$$

增长的是配对的输入特征与输出响应，指数为 $\kappa s_i$；其差对应衰减分支。在多个神经元上，同一组初始系数 $q_i$ 同时出现在相邻权重的低秩项 $q_i v_i^\top$ 和 $u_iq_i^\top$ 中。不能独立地选一组底层特征、再随机选择一组无关的上层列，而声称保留了同一个梯度机制。

这对应原文 Theorem 1 的局部控制：小初始化、短时间、周围网络固定。上层与下层都明显变化以后，$A_\ell$ 随轨迹改变，指数解不能继续当作完整网络解。它与 Saxe 的协作放大有相似代数形式，但这里的算子还包含输入相关的非线性响应。

### 4.3 特征与响应为什么必须成对比较 {#paired-variational}

以 $K(x,x')=\langle z(x),z(x')\rangle$ 和 $G(x,x')=\langle g(x),g(x')\rangle$ 分别定义特征与响应空间。领先奇异对等价于

$$
s_1(A_\ell)=\max_{\|\chi\|_{\mathcal H_K}\le1,\ \|\psi\|_{\mathcal H_G}\le1}
\widehat{\mathbb E}[y\psi(x)\chi(x)].
$$

上层响应 $\psi=\langle u,g\rangle$ 把标签变成有效标签 $y\psi$，底层特征 $\chi=\langle v,z\rangle$ 由其一次相关选出。$\chi$ 或 $\psi$ 单独对 $y$ 都可以不相关，乘积却相关。这是前向 $\mathbb E[y\phi^2]$ 之外的选择原理；但它是否增加可见特征，还取决于当前两个函数空间，不能只凭多了一个因子就保证成功。

## 5. 谱算法与三层自洽动力学 {#backlofi-surrogate}

### 5.1 BackLoFi 算法还有哪些附加选择 {#backlofi-algorithm}

算法先拟合前向 LoFi，再从上到下计算 $A_\ell$ 的主要奇异对，加入配对的行与列修正、重算激活并重新拟合读出。它在已经拟合的网络上仍使用 $y$ 而非真实残差 $y-f$；这在近初始化的局部推导中有依据，在远离初始化处属于算法选择。新的随机提升、行列幅度分别归一化、扰动大小，以及采用精确链式敏感度还是近似敏感度，同样影响算法。

原文近似敏感度还把部分原有传输替换为训练均值，使用修正前的线性负载，并保留当前门控。这些操作有利于可计算和分析，却不是 backprop 的恒等改写。评估“捕捉了反馈机制”与“等价于 SGD”需要完全不同强度的证据。

### 5.2 自洽闭合包含有限维特征与函数态 {#self-consistent}

Theorem 2 研究三层网络的一种**连续谱更新替代过程**：前向滤波已经得到固定维数 $k_1,k_2$；相同选定扰动加入所有第一层滤波坐标；上层前向 LoFi 随下层核重新计算；行、列与读出使用相同微小步长。依次取 $dt\to0$、$p_1\to\infty$、$p_2\to\infty$，用固定 Gaussian 神经元索引 $r,\eta$ 写

$$
h_{1,t}(r,x)=r^\top h^*_{1,t}(x),\qquad
h_{2,t}(\eta,x)=\eta^\top h^*_{2,t}(x)+\mathbb E_r[C_t(\eta,r)z_{1,t}(r,x)].
$$

$C_t$ 保存历次列修正积累的跨层耦合，不能在每步重新独立抽样时被抹去。定义响应函数 $S_t(\eta,s)=w_{3,t}(\eta)\sigma'(s)$，于是 $g_{2,t}(\eta,x)=S_t(\eta,h_{2,t}(\eta,x))$。特征核 $K_{1,t},K_{2,t}$ 和响应核 $G_t$ 由当前激活与敏感度的同索引乘积平均得到。

若当前奇异变分的最优对为 $(\psi_t,\chi_t)$、奇异值为 $\rho_t>0$，则原文式 (17)–(18) 给出

$$
\dot h^*_{1,t,a}(x)=\chi_t(x),\quad a=1,\ldots,k_1,\qquad
h^*_{2,t}=\operatorname{FwdLoFi}_2(K_{1,t},y),
$$

$$
\partial_t C_t(\eta,r)=\frac1{\rho_t}
\widehat{\mathbb E}[y\chi_t g_{2,t}(\eta,X)]\,
\widehat{\mathbb E}[y\psi_t z_{1,t}(r,X)],
$$

$$
\partial_t S_t(\eta,s)=\sigma'(s)\widehat{\mathbb E}[(y-f_t(X))z_{2,t}(\eta,X)].
$$

读出输出还包含核梯度拟合、上层滤波坐标移动、列修正、以及已有耦合对下层激活变化的传输四种贡献，见原文式 (19)。前两式决定表示与耦合，最后一式来自读出拟合；下一时刻的最优对再由更新后的 $K,G$ 选择，形成闭环。

这是一个非平凡的自洽理论：不需把全部有限宽权重保留为状态，也不需逐次记录完整更新历史。但状态仍含 $C_t(\eta,r)$、$S_t(\eta,s)$ 等函数，**不是已经化成几个标量的普适集体动力学**；更不是一般端到端 SGD 的无限宽极限。极限次序、滤波秩固定、谱选择的正则性和原文 Appendix C 的假设都是结论的一部分。

### 5.3 深 staircase 数值例子能检验什么 {#deep-staircase}

论文用 Gaussian 输入的密集二次特征构造可见块 $h_U$ 与前向低阶统计中隐藏的块 $h_V,h_W$，目标包含可见二次组合 $Q(h_U)$、另一组合 $Z(h_U)$ 与隐藏交互 $S(h_V,h_W)$，形如

$$
y_{\rm raw}=\tanh\!\left(\frac{Q+ZS}{\sqrt2}+b_0\right).
$$

通过保留 $Q$、翻转 $Z$ 的对称性安排，可见块先提供入口，隐藏方向却难由普通前向算子取得。学到的上层响应改变有效标签后，隐藏块可以在跨层相关中显露。BackLoFi 与联合训练的比较用于检验这一反馈解释；它是精心构造的机制实验，而不是所有层级任务都需要后向修正的证明。

样本指数 $\log P/\log d$ 上的性能改善说明算法统计可达范围改变；固定 $P$ 下观察新重合随训练时间形成，才涉及动力学出现。两类图不能互换。沿完整轨迹控制经验跨层算子的误差、谱隙与模式切换，仍是把局部机制变成预测性训练时间理论的关键步骤。

## 6. 对 RHM 与集体动力学的含义 {#rhm-and-dynamics}

[RHM 局部学习规则](/notes/hierarchical_data/generative_learning/#local-learning-rules)表明输入相关的非线性 masking 会影响反馈方法能否形成隐层表示；BackLoFi 的 $g_\ell(x)$ 则明确把这类响应放进选择算子。它们都提醒我们只比较标签相关的原始输入统计不够。但前者比较局部目标与反馈近似的实验，后者研究谱替代过程，两者没有共同的通用定理。

BFC 已经证明受控模型里上层会修正下层，BackLoFi 已经提出响应参与选择与自洽核演化。因此“一个特征帮助另一个特征”本身不能作为尚无人研究的新机制。更具体的问题是：耦合何时只改变连续增长率，何时产生显著等待和突发；多个候选模态如何共享残差、争夺容量；有限样本噪声是否使可达方向重新不稳定；以及何时可以从函数态进一步约化到少量重合与响应变量。这些将在[综合比较](/notes/hierarchical_data/mechanism_comparison/)中区分已知、可检验假设与未解决问题。

## 7. 参考文献 {#references}

1. <span id="ref-1"></span>Dandi Y, Vilucchio M, Arnaboldi L, Tabanelli H, Krzakala F（2026）。[Deep Learning as Neural Low-Degree Filtering: A Spectral Theory of Hierarchical Feature Learning](https://arxiv.org/abs/2605.13612v2)。*arXiv:2605.13612，v2，2026-09-04，预印本*。局部梯度展开、前向谱算法、核递归与有效维数。
2. <span id="ref-2"></span>Allen-Zhu Z, Li Y（2023；初稿 2020）。[Backward Feature Correction: How Deep Learning Performs Deep (Hierarchical) Learning](https://proceedings.mlr.press/v195/allen-zhu23a.html)。*COLT，PMLR 195:4598，扩展摘要*；[完整论文 arXiv:2001.04413v6，2023-07-07](https://arxiv.org/abs/2001.04413v6)。正文 §1–4、Theorems 1–2：信息间隙、误差修正与规定 SGD 变体。
3. <span id="ref-3"></span>Dandi Y, Vilucchio M, Tabanelli H, Krzakala F（2026）。[Neural LoFi with Backward Coupling: A Spectral Theory of Cross-Layer Feature Learning](https://arxiv.org/abs/2610.32849v1)。*arXiv:2610.32849，v1，2026-10-23，预印本*。Theorem 1：局部相邻权重对；Algorithm 1：谱修正；Theorem 2 与 Appendix C：三层连续替代动力学。版本状态核查于 2026-10-29。

---

[← 理论 03 · 非线性中间表示](/notes/hierarchical_data/nonlinear_depth/) · [专题目录](/notes/hierarchical_data/#contents) · [下一篇：RHM 01 · 层次生成与表征学习 →](/notes/hierarchical_data/generative_learning/) · [直接进入综合比较 →](/notes/hierarchical_data/mechanism_comparison/)
