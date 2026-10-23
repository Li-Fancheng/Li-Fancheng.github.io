# Grokking 与延迟泛化

网站入口 `/notes/grokking/`，五篇个人学习综述与研究笔记，正文不展示代码、运行命令或练习。

- phenomenology：延迟定义、固定模态、随机设计、阈值和噪声。
- solution_competition：精确竞争模型、乘法动力学、电路效率、结构偏好和优化器。
- representation_circuits：Lazy-to-rich、二次风险、Fourier 算法、RFM 与函数对称投影。
- effective_theory：普通加法约束、零模、保持尺度的有效流及有限时间相图。
- experiments_diagnosis：原始训练与四个探索对照、阈值重算、末态函数干预及语言预训练类比。

## 图与记录

八张图：原有线性延迟、模加法训练、末态对称性三张图完整保留；新增五张独立模型或函数分析图。所有 SVG 与 PNG 并存，模型参数和数值结果见同名 JSON 及 numerical_checks.json。

figures/modular_history.csv、modular_run.json、modular_controls.json、modular_weights.npz 是原训练记录的字节一致副本。未重新运行四万步训练，未补造中间权重。

独立计算：安装 requirements.txt 后运行 `python generate_figures.py`。
公式与记录核查：`python review_checks.py`。

从相邻 neural_network_theory/ 的 modular_grokking.py 可以重新训练原模加法设置，其默认输出已指向本专题 figures/。基础 experiments.py 的历史实验 13 重画线性图，experiments.json 记录 topic/note_slug 和发布路径。extended_checks.py 读取本专题权重并镜像旧对称图。实验编号与网页编号分别维护。

保持尺度的约束模型使用明确推导的 Rayleigh 商流；它是独立机制模型，没有神经解码器。MP 图使用总体连续谱积分；Fourier 图使用显式算法函数。新增函数干预只作用于已保存末态 logits，不是训练时间中的因果干预。投影与主频选择使用完整有限域和已知模加法对称性，不作为新增模型的无偏测试性能。

## 链接与补充包

旧 `/notes/neural_network_theory/grokking/` 转向专题目录，原文章节锚点通过 grokking_routes.json 转向适当的新节。核心神经网络系列现为四卷十篇，编号及文献关联同步更新。

数值材料包保留 neural_network_theory/、scaling_laws/、grokking/ 三个相邻目录，维持跨专题输出路径。三个专题的下载链接提供同一份配套包；各目录 README.md 说明运行方式。

## 文字扩写记录

2026-10-21 在五篇专题完成后追加一轮解释与讨论，中文字符由 17004 增至 27140（增加约 59.6%）。新增内容穿插在原有章节中，展开慢模态、解选择、表示与读出、约束可辨识性和实验解释；五篇公式序列、章节锚点、图与实验数据保持不变。分篇统计见 prose_expansion_report.json；expansion_report.json 保留此前从单篇迁移为专题的历史统计。
