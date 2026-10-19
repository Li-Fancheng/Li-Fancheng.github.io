# 神经网络理论与可解模型

个人学习笔记与文献综述，五卷 20 篇，网站入口 `/notes/neural_network_theory/`。

正文与附件位于本目录，Hugo 页面元信息在 `content/zh/notes/neural_network_theory/`，沿用 `shared-note` 和 `cell-cycle` / `series` 布局。`series.json` 记录篇目与顺序；总目录 `body.md` 说明理论主线、统一记号及可解性层次。正文是完整推导笔记，包含模型假设、数值说明、适用边界和原始文献。

## 复现

```bash
python3 -m venv /tmp/neural-theory-venv
/tmp/neural-theory-venv/bin/pip install -r notes/neural_network_theory/requirements.txt
/tmp/neural-theory-venv/bin/python notes/neural_network_theory/experiments.py --all
/tmp/neural-theory-venv/bin/python notes/neural_network_theory/modular_grokking.py
/tmp/neural-theory-venv/bin/python notes/neural_network_theory/extended_checks.py
```

`--note 1` 到 `--note 20` 可单独重画。25 张模型、训练与函数分析图均提供 SVG、220 dpi PNG，源参数见脚本和图注，数值诊断见各篇 `figures/diagnostics.json`。不含论文原图，不需要图像生成工具、GPU 或下载数据集。

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

20 篇均增加了中间推导和解释：条件风险与噪声传播、逆 Wishart 与 MP resolvent、球面感知机 RS 容量、有限采样核风险、平均场有限矩闭合、非平衡 OU 概率流、Hopfield 零温分支，以及两副本 RNN 稳定性。核心篇进一步连接参数化几何、模加法函数对称性、临界谱的训练与测试指数和标度领先常数。

`extended_checks.py` 独立核对早停最优值、MP 积分、有符号对角模型、两层不平衡闭式解、OU 概率流、Hopfield 容量、时间与数据标度常数，并分析已保存模加法权重；生成四张补充图和 `extended_diagnostics.json`。补充计算不更新网络权重，保留原训练与测试划分。
