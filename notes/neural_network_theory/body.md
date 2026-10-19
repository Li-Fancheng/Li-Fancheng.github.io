这里整理我学习神经网络理论时的文献综述、模型推导和数值记录。整套笔记按五卷、20 篇组织，从线性谱动力学进入高维泛化、无限宽理论、特征学习、集体动力学与标度规律。

我关心的是：神经网络中的哪些现象能够在简单、可解析的模型中得到理解，这些解释依赖什么假设，又在什么地方失效。各篇围绕具体问题连接相关文献，保留关键推导、物理解释、数值结果和仍未厘清的问题。

## 五卷目录 {#contents}

### I · 泛化理论与高维可解模型

| Note | 核心问题 | 理论层次 |
| --- | --- | --- |
| [01 · 线性回归的完整动力学](/notes/neural_network_theory/linear_regression/) | 模态怎样学习？为什么早停可能改善泛化？ | 有限维精确解 |
| [02 · Double Descent](/notes/neural_network_theory/double_descent/) | 插值附近为什么放大噪声？ | Gaussian 精确期望与高维谱 |
| [03 · Random Features](/notes/neural_network_theory/random_features/) | 固定非线性特征能解释多少现象？ | 有限特征求解与核极限 |
| [04 · Teacher–Student](/notes/neural_network_theory/teacher_student/) | 参数恢复怎样联系泛化与样本数？ | Hebbian 渐近解与 RS 容量分析 |

### II · 无限宽网络、核理论与特征学习

| Note | 核心问题 | 理论层次 |
| --- | --- | --- |
| [05 · Gaussian Process](/notes/neural_network_theory/gaussian_process/) | 随机网络怎样定义函数先验？ | 固定深度与有限输入的宽度极限 |
| [06 · Neural Tangent Kernel](/notes/neural_network_theory/ntk/) | 非线性函数为什么出现线性训练？ | 精确核恒等式与条件性冻结极限 |
| [07 · Kernel Learning Curves](/notes/neural_network_theory/kernel_learning_curves/) | 核谱和目标对齐怎样决定误差？ | 可解序列代理与随机设计自洽描述 |
| [08 · Lazy Training 与 Mean Field](/notes/neural_network_theory/lazy_and_mean_field/) | 哪种无限宽极限保留特征移动？ | 粒子到分布的极限理论 |

### III · 深度、优化动力学与表征形成

| Note | 核心问题 | 理论层次 |
| --- | --- | --- |
| [09 · 深度线性网络](/notes/neural_network_theory/deep_linear/) | 相同函数类为何有平台与快慢学习？ | 对齐、平衡分支精确解 |
| [10 · Implicit Bias](/notes/neural_network_theory/implicit_bias/) | 插值解和分类方向由什么选择？ | 线性模型定理与参数化机制 |
| [11 · 非线性特征涌现](/notes/neural_network_theory/nonlinear_teacher_student/) | 隐藏单元怎样向不同教师专门化？ | Gaussian 新样本在线闭合 |
| [12 · 层级组合学习](/notes/neural_network_theory/hierarchical_learning/) | 最终监督怎样传到潜在层级？ | 受限模型精确推导与研究型假设 |

### IV · 记忆、Grokking 与集体动力学

| Note | 核心问题 | 理论层次 |
| --- | --- | --- |
| [13 · Grokking 与延迟泛化](/notes/neural_network_theory/grokking/) | 为什么训练拟合远早于测试改善？ | 线性精确机制与实际模加法训练 |
| [14 · SGD 的统计物理](/notes/neural_network_theory/sgd/) | 随机梯度什么时候能用温度描述？ | 局部 OU 精确解与扩散近似 |
| [15 · Hopfield 记忆](/notes/neural_network_theory/hopfield/) | 模式如何成为吸引子？容量从哪来？ | 能量下降、低负载与自洽近似 |
| [16 · 随机循环网络](/notes/neural_network_theory/random_rnn/) | 非对称随机反馈怎样产生混沌？ | 动力学平均场与有限网络检验 |

### V · Scaling Laws 与损失曲线

| Note | 核心问题 | 理论层次 |
| --- | --- | --- |
| [17 · Time Scaling](/notes/neural_network_theory/time_scaling/) | 指数模态怎样叠加成损失幂律？ | 固定谱精确解与连续谱渐近 |
| [18 · Data Scaling](/notes/neural_network_theory/data_scaling/) | 数据指数怎样来自目标谱与噪声？ | 高斯序列代理的偏差–方差推导 |
| [19 · Model Scaling](/notes/neural_network_theory/model_scaling/) | 模型增大怎样扩展可表示模态？ | 嵌套正交特征的可解模型 |
| [20 · Compute-Optimal Scaling](/notes/neural_network_theory/compute_scaling/) | 固定计算量怎样配置模型与数据？ | 给定损失和成本后的精确优化 |

## 几条贯穿的理论主线

**谱与标度主线：** 01 → 02 → 03 → 06 → 07 → 17 → 18 → 19 → 20。先区分训练与测试，再把优化时间、数据规模、表示容量和计算预算分别定义。

**动力学与特征主线：** 01 → 09 → 10 → 04 → 11 → 08 → 12 → 13。从深度参数化进入重叠动力学，最后比较慢松弛与结构化表示的形成。

**统计物理主线：** 04 → 14 → 15 → 16，再回到 08 与 11。关注序参量、概率流、自洽相关函数与极限条件。

其中 **01、09、06、07、18、13** 连接了谱动力学、深度参数化、核学习曲线与延迟泛化，也是这套笔记中相互联系最紧密的几个主题。

## 统一记号与容易混淆的量 {#notation}

| 记号 | 约定 |
| --- | --- |
| $D$ | 训练样本数；Note 20 中的数据量须另按 token 或独立样本定义 |
| $d$ | 输入维数；不等于模型参数数 |
| $N$ | 隐藏宽度、随机特征数或表示秩，具体含义在每篇说明 |
| $X,H=X^\top X/D$ | 设计矩阵与经验二阶矩 |
| $\Sigma=\mathbb E[xx^\top]$ | 独立测试输入的二阶矩 |
| $\sigma_\varepsilon^2$ 或 $\sigma^2$ | 标签或序列观测噪声方差；Note 09 的 $\sigma_i$ 是信号奇异值 |
| $\eta_i,\psi_i,a_i$ | 总体核本征值、本征函数、目标 $L^2$ 系数 |
| $K_t=J_tJ_t^\top/D$ | 已吸收样本归一化的训练 NTK |
| $\lambda$ | 采用平均平方损失时的 ridge 参数；未平均损失使用 $\Lambda=D\lambda$ |
| $Q,R,T$ | 学生–学生、学生–教师、教师–教师重叠 |
| $\alpha,\beta$ | 按上下文表示高维比例、逆温度或标度指数；每篇重新定义 |

训练损失常带 $1/2$，测试 MSE 通常不带；正文在改变约定时明确写出。时间可以是梯度流时间、在线样本数除以维数、或离散优化步数，三者不能不换算就比较。初始函数、零空间与不可约噪声也始终保留。

## 可解性分四层

**有限模型精确可解**意味着在所列条件下可直接验证闭式解，例如线性梯度流与平衡两层标量网络。**特定极限渐近可解**需要先声明哪些规模趋于无穷、哪些比值固定，例如 MP 谱和 NNGP。**附加假设下闭合**要说明采样与分布条件，例如新样本在线委员会机。**机制模型或有效理论**解释一个可检验的部分现象，例如层级链式瓶颈、Grokking 的某种时间分离与经验损失拟合。

同一篇可能同时出现这四层中的几层，正文会逐步标明边界。一个漂亮的 toy model 不被用来宣称已经解决完整神经网络。

## 数值图与复现实验 {#reproduce}

全系列配有 **25 张数值图，涵盖模型计算、模加法训练和函数结构分析**，记录了本文模型的计算结果，提供 SVG 与 PNG。模加法实验还记录了四个未达到完整泛化的探索对照，避免只展示一个成功配置。

下载：[完整正文、图和复现包](reproducibility.zip) · [20 篇数值脚本](experiments.py) · [补充推导的数值核对](extended_checks.py) · [补充计算结果](extended_diagnostics.json) · [模加法训练脚本](modular_grokking.py) · [依赖列表](requirements.txt) · [篇目配置](series.json) · [文献来源表](sources.json)。

解压后在包根目录运行：

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python experiments.py --all
python experiments.py --note 9
python extended_checks.py
python modular_grokking.py
```

每篇实验写入对应目录的 `figures/`，保存 `experiment.svg`、`experiment.png` 和 `diagnostics.json`。模加法输出训练历史、超参数、权重与固定划分；全部随机实验使用显式种子。浮点库、硬件与依赖版本可能造成细小差异，尤其不能把同一随机种子当成跨平台逐位一致的保证。

图表脚本使用 NumPy、SciPy 与 Matplotlib，模加法无需 PyTorch 或 GPU。这里的轻量实验用来核对理论预测与观察机制；若估计统计规律，应增加种子、样本规模与误差条，并把超参数选择从测试集评价中分离。

## 从这些模型到真实网络

贯穿全系列的联系是：核谱同时影响泛化和时间尺度；深度参数化会改变有效迁移率；特征学习让核与目标对齐一起变化；而记忆、延迟泛化、集体相变和幂律必须分别指定可观测量。看见相似曲线，不足以认定相同微观机制。

研究入口也由此清楚：经验谱何时能代表总体谱？何种任务需要超出固定核？有限样本回用怎样破坏在线闭合？特征形成能否改变长期损失指数？计算最优拟合在哪些数据与架构区间有效？这些问题也构成后续阅读与修订的线索。
