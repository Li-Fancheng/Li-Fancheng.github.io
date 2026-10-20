# Scaling Laws 机制与学习曲线

网站入口 `/notes/scaling_laws/`。六篇个人文献综述与研究笔记按机制组织，正文不展示代码、复现命令或练习。

1. statistical_geometry：统计估计、几何分辨率、宽度换算与数据/模型对偶。
2. spectral_learning：原联合谱模型、时间/数据/容量的完整推导。
3. discrete_structure：技能尾部、Poisson 计数、组合共享与数据选择。
4. feature_learning：时变度量、source 条件、动态加速与 SGD 暂态。
5. graph_language：2026 年 ER/BA 图与语言工作、熵率与计数基线。
6. compute_allocation：成本约束、经验曲面、损失底与协议敏感性。

## 图与计算

共十张自制计算图：五张由原神经网络理论系列迁移，五张为本专题新增。SVG 与 PNG 并存；绘图参数和诊断保存在 figures/ 与 numerical_checks.json。随机图使用独立生成的有限图，未进行 Transformer 训练；特征学习交叉图是有效时钟标度方程示意，不是完整 DMFT 复现。

新增图：安装 requirements.txt 后，在本目录运行 `python generate_figures.py`。
独立公式核查：`python review_checks.py`。覆盖时变度量、平滑计数交叉熵、Poisson 积分、几何偏差方差与计算最优条件。

原五图的基础计算位于相邻 neural_network_theory/：experiments.py 的历史编号 17–20 以及 extended_checks.py 的临界谱计算。历史实验编号稳定，experiments.json 的 topic、note_slug 和 figure 记录当前归属；输出将镜像至本专题 figures/。

## 文献与迁移

sources.json 记录十三篇原始文献及对应正文。注意预印本首发年与会议/期刊发表年不同；2025 年的特征学习模型初稿为 2024，Bahri 等的期刊版为 2024、初稿为 2021。

原 `/notes/neural_network_theory/scaling_laws/` 跳转本专题；time/data/model/compute 四个历史入口直接到对应章节。旧综合稿的六个主要节锚点由专题入口转发。原谱推导和五图全部保留，新增机制内容另行展开。

下载压缩包以 neural_network_theory/ 和 scaling_laws/ 两个相邻目录保存材料，以维持跨专题计算镜像路径。两个专题下载链接提供相同的配套材料包。解压后可只阅读或运行需要的目录。
