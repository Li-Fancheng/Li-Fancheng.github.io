---
title: "神经网络统计物理：文献地图"
subtitle: "Scaling laws、生成模型、特征学习与物理方法"
date: 2026-10-30
lastmod: 2026-11-02
summary: "按科学问题组织 27 篇工作，比较观测量、自由度、理论转化与机制差异。"
layout: cell-cycle
article_class: neural-theory-note
math: true
toc: true
---

这页按研究问题组织文献。每组同时列出观测量、保留的自由度和理论对象，便于从一个解释缺口找到可比较的工作。核心案例与计算见[专题首页](/notes/statistical_physics_of_learning/)。文献中的最优推断、训练轨迹和现象论分别回答不同问题，相同的幂律或相图名称不代表相同机制。

## 资源增加为何改善泛化？

先区分训练后统计估计与训练时间。固定特征使有限样本和有限表示可单独计算，但不自动解释表征适应。

| 工作 | 观测与自由度 | 可计算问题与机制差异 |
| --- | --- | --- |
| Bahri, Dyer, Kaplan, Lee & Sharma，[Explaining neural scaling laws](https://doi.org/10.1073/pnas.2311878121)；*PNAS* **121**, e2311878121 (2024)。 | 数据量、特征数、总体与经验协方差。 | 随机特征回归；区分方差受限与分辨率受限，并联系两种资源轴。 |
| Bordelon, Atanasov & Pehlevan，[A Dynamical Model of Neural Scaling Laws](https://proceedings.mlr.press/v235/bordelon24a.html)；*ICML* (2024)。 | 随机特征模型的时间、宽度与样本复用。 | 学习动力学与模式松弛；说明时间标度和最终资源标度可以不同。 |
| Bordelon, Canatar & Pehlevan，[Spectrum Dependent Learning Curves in Kernel Regression and Wide Neural Networks](https://proceedings.mlr.press/v119/bordelon20a.html)；*ICML* (2020)。 | 核谱、目标投影及随机训练集。 | 连续学习曲线近似与 replica；样本数上的模式学习不同于时间上的指数松弛。 |
| Kramp, Lindner & Helias，[Dynamics of neural scaling laws in random feature regression with powerlaw-distributed kernel eigenvalues](https://proceedings.mlr.press/v306/kramp26a.html)；*ICML* (2026)。 | 固定谱模式、有效噪声和集体响应。 | DMFT 将梯度流、正则化和 Langevin 动力学放在共同描述中。 |

比较时应同时改变目标投影与数据谱，并记录核是否演化。仅匹配误差斜率，难以分辨固定表示的统计限制与训练适应。详见 [scaling law 案例](/notes/statistical_physics_of_learning/spectral_scaling/)。

## 深度、宽度与表示学习各留下什么？

这一组的取舍不同：深度线性保留乘法优化，NTK 保留固定切方向，平均场与比例极限尝试保留表示变化。

| 工作 | 观测与自由度 | 可计算问题与机制差异 |
| --- | --- | --- |
| Saxe, McClelland & Ganguli，[Exact solutions to the nonlinear dynamics of learning in deep linear neural networks](https://arxiv.org/abs/1312.6120)；*ICLR* (2014)。 | 分层矩阵与数据奇异模态。 | 特殊对齐轨道中的非线性 ODE；分离深度参数化与激活非线性。 |
| Bordelon & Pehlevan，[Deep Linear Network Training Dynamics from Random Initialization: Data, Width, Depth, and Hyperparameter Transfer](https://proceedings.mlr.press/v267/bordelon25a.html)；*ICML* (2025)。 | 随机初始化、宽度、数据与层间相关。 | 大宽度训练动力学；把典型随机初始化与特殊精确轨道对照。 |
| Jacot, Gabriel & Hongler，[Neural Tangent Kernel: Convergence and Generalization in Neural Networks](https://papers.nips.cc/paper_files/paper/2018/hash/5a4be1fa34e62bb8a6ec6b91d2462f5a-Abstract.html)；*NeurIPS* (2018)。 | 预测函数及参数梯度的内积。 | 精确函数动力学在特定极限中冻结核，转成谱松弛与核预测。 |
| Chizat, Oyallon & Bach，[On Lazy Training in Differentiable Programming](https://papers.nips.cc/paper_files/paper/2019/hash/ae614c557843b1df326cb29c57225459-Abstract.html)；*NeurIPS* (2019)。 | 参数尺度、曲率和线性化路径。 | 控制局部线性化偏差；缩放可以选择学习机制。 |
| Mei, Misiakiewicz & Montanari，[Mean-field theory of two-layers neural networks: dimension-free bounds and kernel limit](https://proceedings.mlr.press/v99/mei19a.html)；*COLT* (2019)。 | 神经元经验分布及其演化。 | 分布梯度流与近似界；核极限是进一步的尺度选择。 |
| Li & Sompolinsky，[Statistical Mechanics of Deep Linear Neural Networks: The Backpropagating Kernel Renormalization](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.11.031059)；*Physical Review X* **11**, 031059 (2021)。 | 训练后权重的 Gibbs 统计和层核。 | 逐层权重积分求有效核及泛化；平衡统计不同于逐时优化轨迹。 |
| Lauditi, Bordelon & Pehlevan，[Adaptive kernel predictors from feature-learning infinite limits of neural networks](https://proceedings.mlr.press/v267/lauditi25a.html)；*ICML* (2025)。 | 任务适应的层核、前激活分布或 DMFT 固定点。 | 用自洽核表达预测；核形式不必表示固定特征。 |
| Baglioni, Pacelli 等，[Predictive Power of a Bayesian Effective Action for Fully Connected One Hidden Layer Neural Networks in the Proportional Limit](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.133.027301)；*Physical Review Letters* **133**, 027301 (2024)。 | 温度、先验、样本与宽度的比例。 | 以 Langevin 实验检验贝叶斯有效作用量；整体核重标度是所考察范围内的重要效应。 |
| Ingrosso, Pacelli, Rotondo & Gerace，[Statistical Mechanics of Transfer Learning in Fully Connected Networks in the Proportional Limit](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.134.177301)；*Physical Review Letters* **134**, 177301 (2025)。 | 源—目标相关及比例极限。 | 单实例 Franz–Parisi 形式得到跨任务有效核，分析迁移收益。 |
| Barbier 等，[Statistical Physics of Deep Learning: Optimal Learning of a Multilayer Perceptron near Interpolation](https://journals.aps.org/prx/abstract/10.1103/56sb-pdh6)；*Physical Review X* **16**, 031014 (2026)。 | 教师结构、矩阵序参量与特化。 | replica 与随机矩阵方法计算最优推断；不能直接当作有限时间 SGD 的轨迹。 |

一个具体缺口是，相同测试误差能否对应完全不同的表示几何。应同时比较对齐、核方向及算法可达性，而非仅按是否存在闭式解评价模型。参见[深度线性](/notes/statistical_physics_of_learning/deep_linear_model/)与 [NTK/平均场](/notes/statistical_physics_of_learning/changing_variables/)。

## 数据结构怎样进入学习阈值？

固定一个网络仍不足以定义层次或语义；生成规则与任务结构先决定有哪些信号可供学习。

| 工作 | 观测与自由度 | 可计算问题与机制差异 |
| --- | --- | --- |
| Cagnetta, Petrini, Tomasini, Favero & Wyart，[How Deep Neural Networks Learn Compositional Data: The Random Hierarchy Model](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.14.031001)；*Physical Review X* **14**, 031001 (2024)。 | 同义展开、层级规则、局部相关与不变性。 | 组合概率和信号可检测性，与网络表示实验对照。 |
| Cui, Behrens, Krzakala & Zdeborová，[A Phase Transition between Positional and Semantic Learning in a Solvable Model of Dot-Product Attention](https://papers.nips.cc/paper_files/paper/2024/hash/3fefebc2d4e3c1c6ee9b892bd293117d-Abstract-Conference.html)；*NeurIPS* (2024)。 | 低秩共享 query/key 与受控 token 数据。 | 高维全局最优解在位置和内容策略间转变；模型中的语义具有指定含义。 |
| Erba, Troiani, Biggio, Maillard & Zdeborová，[Bilinear Sequence Regression: A Model for Learning from Long Sequences of High-Dimensional Tokens](https://journals.aps.org/prx/abstract/10.1103/l4p2-vrxt)；*Physical Review X* **15**, 021092 (2025)。 | 长序列、高维 token 与双线性结构。 | 最优泛化、消息传递及梯度算法比较，区分信息阈值与具体训练。 |

若算法学习早于选定统计量的检测阈值，可能是观测量遗漏了其他信号；若晚于阈值，则需定位算法障碍。这两种情况不能都称为数据不足。参见 [RHM](/notes/statistical_physics_of_learning/designing_hierarchy/)与[随机规则计算](/notes/statistical_physics_of_learning/hierarchy_exercise/)。

## 随机优化需要多少历史信息？

局部稳态、在线新样本和有限数据复用分别保留不同相关，不能由一个噪声强度统一替代。

| 工作 | 观测与自由度 | 可计算问题与机制差异 |
| --- | --- | --- |
| Mandt, Hoffman & Blei，[Stochastic Gradient Descent as Approximate Bayesian Inference](https://jmlr.org/papers/v18/17-214.html)；*Journal of Machine Learning Research* **18**, 134 (2017)。 | 局部 Hessian 与梯度协方差。 | OU 稳态、协方差与后验近似，研究稳定区域内的涨落。 |
| Mignacco, Krzakala, Urbani & Zdeborová，[Dynamical mean-field theory for stochastic gradient descent in Gaussian mixture classification](https://proceedings.neurips.cc/paper_files/paper/2020/hash/6c81c83c4bd0b58850495f603ab45a93-Abstract.html)；*NeurIPS* (2020)。 | 固定随机数据、采样过程、两时相关与响应。 | 无序平均导向带记忆的自洽场过程，计算整条训练轨迹。 |
| Mignacco & Urbani，[The effective noise of stochastic gradient descent](https://doi.org/10.1088/1742-5468/ac841d)；*Journal of Statistical Mechanics: Theory and Experiment* **2022**, 083405。 | 涨落—响应及不同噪声实现的轨迹差。 | 在指定稳态和模型中量化有效噪声；不是预先赋予 SGD 普适温度。 |
| Biehl 等，[Dynamics and Generalization Ability of LVQ Algorithms](https://www.jmlr.org/beta/papers/v8/biehl07a.html)；*Journal of Machine Learning Research* **8**, 323–360 (2007)。 | 原型与类别方向的重叠、范数。 | 高维在线独立输入使序参量闭合，连接宏观轨迹与泛化。 |

保持单时噪声协方差近似相同而改变采样相关时间，可以检验历史信息是否必要。参见 [SGD](/notes/statistical_physics_of_learning/stochastic_dynamics/)与[在线重叠闭合](/notes/statistical_physics_of_learning/online_order_parameters/)。

## 生成结构和样本记忆为何分开？

给定经验 score 后，问题转向加噪密度的能量统计与生成轨迹。

| 工作 | 观测与自由度 | 可计算问题与机制差异 |
| --- | --- | --- |
| Biroli, Bonnaire, De Bortoli & Mézard，[Dynamical regimes of diffusion models](https://www.nature.com/articles/s41467-024-54281-3)；*Nature Communications* **15**, 9957 (2024)。 | 协方差集体方向、样本能量和归一化权重。 | 谱分化与无序凝聚分别组织类别选择和样本坍缩；精确经验 score 是核心理想化。 |

可比较学习 score 偏离经验 score 的方式，观察类别恢复和样本凝聚是否分别移动。具体配分函数与典型能量计算见[扩散案例](/notes/statistical_physics_of_learning/diffusion_score/)。

## 相变究竟改变了哪个观测？

传播稳定性与约束可满足性都可出现阈值，但控制参量、序参量和与泛化的关系不同。

| 工作 | 观测与自由度 | 可计算问题与机制差异 |
| --- | --- | --- |
| Schoenholz, Gilmer, Ganguli & Sohl-Dickstein，[Deep Information Propagation](https://arxiv.org/abs/1611.01232)；*ICLR* (2017)。 | 随机初始化下的方差、输入相关与梯度传播。 | 平均场递推的固定点稳定性给出深度尺度，临界性主要描述初始化传播。 |
| Geiger, Spigler, d’Ascoli, Sagun, Baity-Jesi, Biroli & Wyart，[Jamming transition as a paradigm to understand the loss landscape of deep neural networks](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.012115)；*Physical Review E* **100**, 012115 (2019)。 | 样本约束、可调方向及 Hessian 谱。 | 以堵塞比较插值阈值和损失景观；训练约束满足不等于泛化转变。 |

区分机制需要测传播相关长度、约束边界或曲率，而不能只把误差下降称为相变。临界初始化是否也预测训练后的表示，需要另外检验。

## 哪些自由度可以被积分掉？

场论和 RG 必须规定场、尺度与保留的观测，逐层权重积分也不自动等同于 Wilson 模态积分。

| 工作 | 观测与自由度 | 可计算问题与机制差异 |
| --- | --- | --- |
| Cohen, Malka & Ringel，[Learning curves for overparametrized deep neural networks: A field theory perspective](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.3.023034)；*Physical Review Research* **3**, 023034 (2021)。 | 网络诱导的函数分布与训练数据约束。 | 场论、数据平均和图展开组织学习曲线，与多项式回归实验比较。 |
| Howard, Jefferson, Maiti & Ringel，[Wilsonian renormalization of neural network Gaussian processes](https://doi.org/10.1088/2632-2153/adc8fc)；*Machine Learning: Science and Technology* **6**, 025038 (2025)。 | 核模态截止与可学习方向。 | 积分当前难以学习的模态，追踪有效 ridge 等参数的流。 |

若训练显著改变特征，固定核模态的相关与非高斯修正是否仍可吸收到少数有效参数，是比笼统宣称网络等于 RG 更具体的问题。

## 微观推导不可得时，怎样约束宏观模型？

先定义网络可测量，再提出相互作用；拟合相图只是证据的一部分。

| 工作 | 观测与自由度 | 可计算问题与机制差异 |
| --- | --- | --- |
| Shi, Pan & Dokmanić，[Spring-Block Theory of Feature Learning in Deep Neural Networks](https://journals.aps.org/prl/abstract/10.1103/ys4n-2tj3)；*Physical Review Letters* **134**, 257301 (2025)。 | 层间数据分离、噪声与非线性。 | 随机弹簧—滑块链组织宏观相图，属于候选现象论。 |

需要用未拟合的负载、噪声或层间扰动检验关联预测，避免每条曲线独立调参。参见[层间现象论](/notes/statistical_physics_of_learning/phenomenological_model/)。
