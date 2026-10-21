# 神经网络理论与可解模型

个人学习笔记与文献综述，四卷 10 篇，网站入口 `/notes/neural_network_theory/`。

正文与附件位于本目录，Hugo 页面元信息在 `content/zh/notes/neural_network_theory/`，沿用 `shared-note` 和 `cell-cycle` / `series` 布局。`series.json` 记录 10 篇网页文章与顺序，`experiments.json` 独立记录 20 个基础数值实验；总目录 `body.md` 说明理论主线、统一记号及可解性层次。正文是完整推导笔记，包含模型假设、数值说明、适用边界和原始文献。

## 复现

```bash
python3 -m venv /tmp/neural-theory-venv
/tmp/neural-theory-venv/bin/pip install -r notes/neural_network_theory/requirements.txt
/tmp/neural-theory-venv/bin/python notes/neural_network_theory/experiments.py --all
/tmp/neural-theory-venv/bin/python notes/neural_network_theory/modular_grokking.py
/tmp/neural-theory-venv/bin/python notes/neural_network_theory/extended_checks.py
```

数值实验沿用原编号：`--note 1` 到 `--note 20` 可单独重画。17–20 现归属独立的相邻 `scaling_laws/` 专题，并不代表四篇网页。原 25 张模型、训练与函数分析图均提供 SVG、220 dpi PNG，源参数见脚本和图注，数值诊断见各篇 `figures/diagnostics.json`。不含论文原图，不需要图像生成工具、GPU 或下载数据集。

模加法默认平方激活、模数 19、隐藏宽度 128、训练占比 0.7、AdamW 学习率 0.003、解耦衰减 1、40000 步、种子 7。训练/测试划分和权重保存在 NPZ；完整 CSV 每 100 步记录交叉熵、准确率、范数和初步 Fourier 诊断。`modular_controls.json` 保留另外四次探索的结果；成功示例经过设置选择，不是无调参基准，也不声称所有随机种子都能 Grok。

## 数学与证据边界

- 线性回归采用平均平方损失；NTK 训练矩阵也除以 D。
- Double descent 的 Gaussian 公式明确固定信号范数，未伪造欠参数侧的经典下降段。
- 核学习曲线与数据标度中的高斯序列代理，与真实随机输入采样分开；保留 ridge 饱和条件。
- 深度线性闭式解限定平衡与对齐；矩阵因子化不写成无条件最小核范数定理。
- 委员会机保留在线学习率平方项；层级模型明确哪些特征已经由架构提供。
- Grokking 的线性玩具模型使用指定经验/测试度量，实际模加法另行训练。
- OU、Hopfield 低负载和 RNN DMFT 分别说明扩散近似、负载与耦合假设。
- Compute scaling 的优化是损失与成本模型成立后的条件性结论。

`sources.json` 提供检索日期、原始论文链接与各篇用途。自制推导及实验不宣称逐图复刻原论文。

## 网页维护

`data/neural_theory_figures.json` 为每张图声明宽度；现有图片 hook 负责原尺寸链接。所有正文在顶部、底部链接前后篇与总目录。主页和 notes 总目录提供入口。复现 ZIP 包含本目录除 ZIP 自身和缓存外的全部正文、脚本、图和数值记录，修改后需重新打包。

比较符号在数学源码中使用 `\lt` / `\gt`，避免 Markdown 把紧接字母的裸 `<` 误判为 HTML。页面的 `reading_time` 按中文长度及展示公式数量给出粗略阅读估计，供长文导航使用；逐步核对推导所需时间会更久。通用布局仅在存在该字段时覆盖 Hugo 对共享正文的默认估计。

## 展开后的推导与补充计算

原 20 篇均增加了中间推导和解释，现按 17 篇组织：条件风险与噪声传播、逆 Wishart 与 MP resolvent、球面感知机 RS 容量、有限采样核风险、平均场有限矩闭合、非平衡 OU 概率流、Hopfield 零温分支，以及两副本 RNN 稳定性。核心篇进一步连接参数化几何、模加法函数对称性、临界谱的训练与测试指数和标度领先常数。

`extended_checks.py` 独立核对早停最优值、MP 积分、有符号对角模型、两层不平衡闭式解、OU 概率流、Hopfield 容量、时间与数据标度常数，并分析已保存模加法权重；生成四张补充图和 `extended_diagnostics.json`。补充计算不更新网络权重，保留原训练与测试划分。

## 完整综述扩写

本轮以用户提出“每篇继续扩到 2–3 倍”前的版本为基线，逐篇补入背景、模型选择的来由、跨模型推导、解释与结语。`expansion_report.json` 保留合并前的历史扩写统计；当前 10 篇结构与长度见 `consolidation_report.json`。这里的长度统计不是单词数或学习时间；页面阅读时长仍按中文文字和展示公式估计。

新增的独立计算见 `review_checks.py` 与 `review_diagnostics.json`：固定惩罚的留一 ridge、GP 边缘似然导数、NTK 初始化协方差、Bayesian 秩一更新、深度线性 Gram 守恒、OU 概率流示例、模加法群投影与凸交叉熵、局部场噪声教师，以及额外训练轮数的计算最优条件。它们不更新任何网络权重。

```bash
python notes/neural_network_theory/review_checks.py
```

25 张独立图全部保留；深度与隐式偏置合并后，初始化偏置图只在相应推导处展示一次，网页共有 25 个图位置。所有正文与新核对结果一起进入下载 ZIP。

## 网页正文的展示约定

网页仅呈现综述、推导、公式、图表、实验参数与结果解释。运行命令、代码块、脚本文件名、依赖安装和输出目录说明不进入系列目录或文章正文；技术说明集中在本 README 中。

合并前清理了 20 篇中的运行命令和目录页代码块，并将图注、训练记录说明改写为正文表述。长度记录按清理后的正文更新，未为满足旧长度数字补入文字。

## Scaling Laws 合并

原 17–20 篇合为 `scaling_laws/body.md`，这一轮将网页系列从 20 篇收拢为 17 篇。新文以联合谱风险开篇，依次讨论时间、数据、模型、计算、经验外推与指数之间的联系；删去重复前言、尾部推导和跨篇复述，保留关键推导与五张图。旧四个网页入口使用不进入目录的重定向，分别落到新文的对应章节。

旧 `time_scaling`、`data_scaling`、`model_scaling`、`compute_scaling` 文件夹只保存原实验输出。基础实验目前将结果同步到相邻 `../scaling_laws/figures/`，补充计算同样同步临界谱图。实验编号与网页篇目分别维护，合并不改变数值设置。

统一谱风险的新增计算以系数 ODE 和独立带噪观测分别核对，结果记录在 `review_diagnostics.json` 的 `joint_time_data_model_risk` 中。

## 前半部分进一步合并

系列从 17 篇整理为 12 篇：线性回归与双下降、随机特征与核学习曲线、GP 与 NTK、深度线性与隐式偏置，以及两篇 Teacher–Student 分别合并。合并稿重写前言、过渡与结语，保留核心推导和全部独立数值图；原网页地址作为不进入目录的跳转入口。

网页编号按 01–12 重新排序，文献关联、前后篇导航、跨篇引用和首页同步更新。实验仍使用 1–20 的历史编号，`experiments.json` 的 `note_number` 与 `figure` 明确当前文章和发布图片路径。旧专题目录保留原始实验输出，不再保存独立正文。

历史扩写长度保留在 `expansion_report.json`，每轮合并和当前正文统计见 `consolidation_report.json`。不会为维持旧篇目或长度重复插入内容。

所有 25 张网页图声明展示宽度与固有宽高，SVG 读取前即预留空间，减少长文目录和旧地址跳转后的布局位移。

## 2026-10-20 Scaling Laws 独立专题

当前网页为四卷十一篇。原第十二篇迁入 `/notes/scaling_laws/`，按统计估计、谱、离散技能、特征学习、随机图和计算配置扩展成六篇。原五张 Scaling 图迁入独立专题，主系列现展示二十张图；历史数值输出仍保留。

基础实验 17–20 在 experiments.json 使用 topic/note_slug 标记归属，公开镜像指向相邻目录；其余实验的 note_number 保持不变。extended_checks.py 同步更新临界谱图路径。上文二十篇、十七篇和十二篇的说明为历史合并记录。

补充包现在包含 neural_network_theory/ 与 scaling_laws/ 两个相邻目录，使跨专题导出路径在解压后仍成立。旧网页与节锚点的跳转和当前页面校验见 verification.json。

## 2026-10-21 Grokking 独立专题

当前网页为四卷十篇。原第八篇迁入 `/notes/grokking/` 并扩展为五篇机制与实验笔记，保留三张旧图与全部已有训练记录，新增五张计算或函数分析图。核心系列现在展示十七张图；后续 SGD、Hopfield、RNN 编号顺延为 08、09、10。

历史基础实验 13 通过 experiments.json 指向 `../grokking/figures/linear_delay`。模加法训练的默认输出和 extended_checks.py 的权重读取、对称图发布都转向相邻 grokking/；旧原始输出文件保留。补充包加入第三个相邻目录，历史网页锚点转发到新系列对应内容。
