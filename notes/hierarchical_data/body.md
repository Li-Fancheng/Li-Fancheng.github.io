一张图像或一句话可以有指数多个表面实现，却由有限、可复用的局部规律生成。学习器怎样发现这些规律，把它们组织成表示，并在不同语言或模态之间保留共同关系？这个专题从层次结构数据出发，连接统计学习、生成推断和表示几何。

五篇笔记按问题展开，并采用篇内连续的章节编号。每篇从直观例子建立对象，再展开公式、中间推导与可检验的讨论，篇末给出独立参考文献和前后篇导航。精确概率计算、均场近似、简化算例与原论文的网络实验分别说明，便于判断一个结论依赖哪些假设。

## 1. 五篇笔记 {#contents}

| Note | 主要问题 | 关键计算 |
| --- | --- | --- |
| [01 · 层次生成与表征学习](/notes/hierarchical_data/generative_learning/) | 不同表面实现怎样被认成同一个部件？ | 小语法、标签统计、样本阈值与可证明逐层训练 |
| [02 · 语言统计与学习曲线](/notes/hierarchical_data/language_statistics/) | 更远的上下文何时真正有用？ | 关联平台、条件熵分解、有限样本分辨率与长尾积分 |
| [03 · 层级推断与扩散生成](/notes/hierarchical_data/inference_generation/) | 加噪以后哪些层次仍能恢复？ | 树上后验、联合采样、均场固定点与动力学关联 |
| [04 · 语法解析与隐变量学习](/notes/hierarchical_data/parsing_latents/) | 组合边界未知时怎样学习？ | 多重解析、区间动态规划、隐变量预测与可辨识性 |
| [05 · 跨语言表示与共享几何](/notes/hierarchical_data/shared_representations/) | 翻译句为何在中间层更相似？ | 共享语法、几何与坐标、邻域指标、仿射残差与概念谱 |

## 2. 先分清生成规律、推断与学习 {#starting-point}

可以把整个专题的任务放在同一个情境里：生成者用一套固定语法不断产生句子，学习者只看到最终字符串，有时还看到类别标签。研究者知道内部树，所以可以计算理想答案，也可以精确控制哪些因素改变。这个便利让许多在真实语言中混在一起的问题，能够被单独检验。

| 层面 | 已知与未知 | 需要回答的问题 |
| --- | --- | --- |
| 生成模型 | 给定字典、规则概率与树的组织 | 有多少合法样本，哪些变化共享同一原因？ |
| 已知规则的推断 | 规律已知，这次样本的隐状态未知 | 不完整观测支持哪些解释，它们各有多大概率？ |
| 有限数据学习 | 只能看到样本，字典或表示未知 | 哪些统计能暴露规律，多少样本足够？ |
| 网络实现 | 还指定了架构、参数与训练方法 | 所需统计能否被取得，表示能否逐层复用？ |
| 表示评价 | 只能通过读出、距离或干预观察内部状态 | 共享的是局部关系、坐标方向，还是实际功能？ |

生成空间很大，不代表必须记住每条样本；字典很小，也不代表一定容易发现。后验能够精确计算，不代表网络已经学会字典；表示邻域很相似，也不代表探针能够直接跨语言迁移。这些联系构成五篇的共同问题，而不是五个彼此独立的应用。

### 2.1 记号与默认条件 {#notation}

| 记号 | 含义 | 比较公式时需要留意 |
| --- | --- | --- |
| $L,\ell$ | 数据树深与节点高度；叶层为 $0$ | 原论文可能从根向下编号 |
| $s,d=s^L$ | 固定分支数与叶序列长度 | 第四篇的可变分支使句长成为随机量 |
| $v,m,n_c$ | 层内词表、每父规则数与根类别数 | 部分分析取 $n_c=v$，并要求局部无歧义 |
| $\rho=m/v^{s-1}$ | 规则对完整子元组空间的占有率 | 语言统计主要考虑 $\rho<1$ |
| $P$ | 独立生成的完整训练序列数 | 不直接等于相关 token 的数量 |
| $k,t$ | 网络计算层与噪声时间或上下文长度 | $t$ 的具体含义在各篇重新定义 |

默认先抽取一套随机字典，再在训练与测试中固定；同层规则跨位置共享，不同完整样本独立生成。前三篇先研究已知固定分块，第四篇才去掉边界提示，第五篇再把单语言改成共享高层、各自实现低层的多语言模型。改变这些条件时，样本标度、推断算法与可辨识性都需要重新检查。

## 3. 几条贯穿的联系 {#connections}

**规则与统计信号。** 生成树提供可复用的对象，但不保证它们能从有限数据辨认。第一篇建立同义关系，第二篇追踪信号怎样随距离衰减，第四篇再讨论边界和目标也未知时的困难。

**统计阈值与训练保证。** [可证明的逐层串联](/notes/hierarchical_data/generative_learning/#shallow-to-deep-chaining)将局部标签信号、可复用表示和误差控制接成深层学习保证，重点补充 [Ren、Dandi、Krzakala、Lee 的 COLT 2026 工作](https://proceedings.mlr.press/v336/ren26a.html)。这条监督学习路线与第四篇的隐变量预测使用不同目标和恢复条件，需要分别比较。

**推断与学习。** 第三篇先假设字典已知，计算真正的后验；其他各篇讨论如何从数据接近这个参照。最优后验仍有歧义，与模型没有学会正确后验，是两种不同限制。

**共享内容与表面实现。** 第五篇将前面的单语言模型扩展到共享高层规则的语言。重点纳入 [A theory of platonic representations in language models](https://arxiv.org/abs/2610.07168)，并结合概念层次、统计对称性与相似性校准，区分共享几何、坐标对齐和实际计算。

所有笔记使用叶层为 $0$、根层为 $L$ 的记号；原论文的层编号不同，在比较标度公式时重新对齐。数据树的层次、网络计算层和训练时间分别定义。

## 4. 怎样按问题阅读 {#reading-paths}

第一次阅读建议按 01–05 的顺序：先建立完整语法与同义关系，再看有限样本统计，随后计算已知规律下的最优后验，最后进入未知边界和共享表示。

如果关心学习曲线，可以先读 01 的样本阈值与 02 的上下文和长尾分析；如果关心扩散，先读 01 的生成模型，再进入 03 的 BP、噪声通道与生成评价；如果关心跨语言表示，读过 03 的上下行消息后进入 05，会更容易理解编码、解码与中间层共享的含义。04 则解释固定分块为什么是一个强假设，以及隐变量目标怎样改变学习信号。

每篇的结尾都把结论接向下一步。读公式时可沿“对象是什么—哪些条件已知—统计怎样估计—误差如何汇总”检查；读几何图时再补上“允许哪些坐标变换—使用哪个空模型—是否验证功能”三项问题。

## 5. 与其他笔记的联系 {#perspective}

[神经网络理论与可解模型](/notes/neural_network_theory/)提供核方法、特征学习和深度动力学背景；其中的[层级组合学习](/notes/neural_network_theory/hierarchical_learning/)用奇偶机制模型分析弱信号和链式梯度，与这里典型随机语法的统计信号形成对照。[Scaling Laws 机制与学习曲线](/notes/scaling_laws/)进一步比较统计估计、谱结构、离散技能和动态特征的标度机制，[扩散模型的动力学](/notes/diffusion_dynamics/)则从生成过程的时间演化进入。

## 6. 文献与研究脉络 {#references}

各篇在相关推导与结论旁链接原文，并在篇末列出作者、年份和出处。这里保留总表，便于沿同一问题比较不同模型；预印本与正式发表版本按同一项工作处理。

| 问题 | 原始文献 | 关注的联系 |
| --- | --- | --- |
| 层级生成与可学习性的早期形式化 | [Deep Learning and Hierarchal Generative Models](https://arxiv.org/abs/1612.09057)，2016 | 区分层级信息和受限低阶统计所能实现的分类 |
| 逐层学习算法 | [A Provably Correct Algorithm for Deep Learning that Actually Works](https://arxiv.org/abs/1803.09522)，2018 | 局部梯度更新、聚类与层级恢复的算法路线 |
| 随机语法与监督分类 | [How Deep Neural Networks Learn Compositional Data: The Random Hierarchy Model](https://arxiv.org/abs/2307.02129)，PRX 2024 | 随机层级模型、同义不变性和样本阈值 |
| RHM 的可证明深层学习 | [Provable Learning of Random Hierarchy Models and Hierarchical Shallow-to-Deep Chaining](https://proceedings.mlr.press/v336/ren26a.html)，COLT 2026 | 根标签监督、RBF 随机特征、逐层梯度训练与误差串联 |
| 稀疏性与空间稳定性 | [How Deep Networks Learn Sparse and Hierarchical Data](https://arxiv.org/abs/2404.10727)，ICML 2024 | 空白位置、离散形变与权重共享 |
| 语言关联与有效上下文 | [Towards a Theory of How the Structure of Language Is Acquired by Deep Neural Networks](https://arxiv.org/abs/2406.00048)，NeurIPS 2024 | 共同祖先、关联平台与有限样本分辨率 |
| 推断与卷积架构 | [U-Nets as Belief Propagation](https://arxiv.org/abs/2404.18444)，2024 年初稿 | 层级生成模型中的推断近似、分类与去噪 |
| Transformer 的层级计算 | [How Transformers Learn Structured Data](https://proceedings.mlr.press/v267/garnier-brun25a.html)，ICML 2025 | 用层级过滤检验不同关联尺度的恢复 |
| 最优去噪中的重建转变 | [A Phase Transition in Diffusion Models Reveals the Hierarchical Nature of Data](https://arxiv.org/abs/2402.16991)，PNAS 2025 | 前向—后向实验与不同层次的重合率 |
| 用协同变化探测隐结构 | [Probing the Latent Hierarchical Structure of Data via Diffusion Models](https://arxiv.org/abs/2410.13770)，ICLR 2025 | 动力学关联、易感性和临界长度 |
| 长尾规则与学习曲线 | [Learning Curves Theory for Hierarchically Compositional Data with Power-Law Distributed Features](https://proceedings.mlr.press/v267/cagnetta25a.html)，ICML 2025 | 规则稀有性、层级前因子与任务依赖 |
| 架构与语言标度 | [Scaling Laws and Representation Learning in Simple Hierarchical Languages](https://arxiv.org/abs/2505.07070)，PRE 2025 | 局部性、共享与高层规则的统计效率 |
| 有限数据下的生成学习 | [How Compositional Generalization and Creativity Improve as Diffusion Models Are Trained](https://arxiv.org/abs/2502.12089)，ICML 2025 | 从局部合法到全局连贯的逐尺度生成 |
| 可变树与解析歧义 | [Deep Networks Learn to Parse Uniform-Depth Context-Free Languages from Local Statistics](https://arxiv.org/abs/2602.06065)，ICML 2026 | 不固定组合边界时的统计学习 |
| 隐变量预测的效率 | [Learn from Your Own Latents and Not from Tokens](https://arxiv.org/abs/2605.27734)，2026 年预印本 | 将预测目标提升到隐层次及其适用条件 |
| 跨语言共享表示 | [A theory of platonic representations in language models](https://arxiv.org/abs/2610.07168)，2026 年 10 月预印本 | 共享高层语法、相似性与对齐、逐层仿射残差 |
| 表示收敛假说 | [Position: The Platonic Representation Hypothesis](https://proceedings.mlr.press/v235/huh24a.html)，ICML 2024 | 跨模型与跨模态表示收敛的经验假说 |
| 表示相似性的统计校准 | [Revisiting the Platonic Representation Hypothesis: An Aristotelian View](https://arxiv.org/abs/2602.14486)，ICML 2026 | 置换空模型、维度与选层效应、局部关系 |
| 概念层次与共现谱 | [Hierarchical Concept Geometry in Language Models Emerges from Word Co-occurrence](https://arxiv.org/abs/2605.23821)，2026 | 概念分类树中的由粗到细谱结构 |
| 连续统计对称性与几何 | [Symmetry in language statistics shapes the geometry of model representations](https://arxiv.org/abs/2602.15029)，2026 | 共现对称性、Fourier 模式与连续几何 |
| 条件统计的集中界 | [Hoeffding：Probability Inequalities for Sums of Bounded Random Variables](https://doi.org/10.1080/01621459.1963.10500830)，1963 | 规则出现次数、信号分离度与估计精度 |
| 概率图模型与消息传递 | [Wainwright 与 Jordan：Graphical Models, Exponential Families, and Variational Inference](https://doi.org/10.1561/2200000001)，2008 | 树上的边缘后验与 sum-product |
| 语法概率估计 | [Lari 与 Young：The estimation of stochastic context-free grammars using the Inside-Outside algorithm](https://doi.org/10.1016/0885-2308(90)90022-X)，1990 | 软解析、外部权重与规则期望计数 |
| 正交坐标对齐 | [Schönemann：A generalized solution of the orthogonal Procrustes problem](https://doi.org/10.1007/BF02289451)，1966 | 允许正交变换后的最小二乘匹配 |

---

从第一篇开始：[01 · 层次生成与表征学习 →](/notes/hierarchical_data/generative_learning/)
