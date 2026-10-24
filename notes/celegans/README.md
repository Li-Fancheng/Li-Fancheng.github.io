# 线虫专题

以秀丽隐杆线虫（Caenorhabditis elegans）为模型的研究笔记与文献综述。

胚胎形态发生系列保留原入口 `/notes/celegans/embryonic_morphogenesis/`，分为七篇：

1. [细胞内移与力传递](embryonic_morphogenesis/gastrulation/body.md)
2. [表皮重排与腹侧闭合](embryonic_morphogenesis/epithelial_morphogenesis/body.md)
3. [早期伸长与各向异性力学](embryonic_morphogenesis/early_elongation/body.md)
4. [肌肉驱动与形态记忆](embryonic_morphogenesis/late_elongation/body.md)
5. [三维体形与模型辨识](embryonic_morphogenesis/geometry_and_inference/body.md)
6. [Active matter 软物质堆积与拓扑](embryonic_morphogenesis/active_matter_packing_topology/body.md)
7. [亚细胞骨架与流动](embryonic_morphogenesis/subcellular_cortex_and_flows/body.md)

正文通过 `shared-note` 与长篇布局展示，总览使用按 weight 排序的系列布局。各篇有独立目录、局部图号、文献和上一篇／下一篇链接。论文与计算图集中保存在 `embryonic_morphogenesis/figures/`，避免重复素材。

68 幅论文原图保留面板和标注，来源、权利信息、校验值与逐篇图号见 `embryonic_morphogenesis/media_sources.json`；91 条文献元数据见 `literature_sources.json`。显示宽度与原始尺寸见 `data/celegans_figures.json`，网页保留原图链接。

33 幅计算图使用文中明确给出的几何关系或无量纲模型。生成脚本为 `figures/generate_figures.py`、`figures/generate_model_figures.py` 与 `figures/generate_extended_models.py`、`figures/generate_active_packing_figures.py`、`figures/generate_worm_mechanics_figures.py`、`figures/generate_solid_fluid_figures.py`、`figures/generate_epithelial_spatial_models.py`、`figures/generate_subcellular_models.py`、`figures/generate_morphogenesis_network_models.py`，需要 NumPy、Matplotlib；扩展模型、表皮空间与亚细胞模型脚本另需 SciPy。论文原图不经过这些脚本修改。

各阶段开头嵌入 7 段原始显微影片，第二篇分别展示背侧插入、腹侧包裹和腹侧裂隙封闭；原有激光切割影片保留。影片直接使用研究机构或期刊的远程 MP4，保留播放控制、采样间隔和原始来源，记录见 media_sources.json 中的 embedded_stage_movies。竖向影片限制显示宽度，避免在长篇正文中过度放大。

总览补充胚胎结构、细胞谱系、荧光与 DIC 成像、RNAi、恢复实验和模型近似的说明。七篇正文在术语首次出现处解释含义，并展开从影像测量到模型变量的对应关系；具体假设和实验局限保留在相关段落。

模型部分按尺度与机制比较多条路线，新增三维细胞形状、相场、主动凝胶、连接受力失效、三维插入、非线性纤维、孔弹性、周期更新、Cosserat 杆与参数信息分析。新增图均有明确方程和无量纲参数；候选框架与已发表线虫模型、预印本结果分别说明。新增两幅论文图展示相场与皮层反馈模型。

第六篇比较受限软颗粒、主动 Brownian 粒子、顶点／Voronoi、三维细胞、主动皮层与向列模型。正文进一步讨论线虫四细胞排列、界面张力、内移的堆积与 clutch、手性皮层流、背侧定向插入、腹侧包裹与连接成熟、actin 各向异性、周期驱动与材料记忆、有效力反演；含 43 篇文献、15 幅论文图、9 幅计算或结构图，以及线虫皮层流和培养细胞取向运动的 2 段远程影片。新增算例核查双势阱失稳、邻居松弛的 Schur 补、clutch 力平衡、固定邻接的仿射形变和周期加载后的参考应变。

第六篇第 9 节集中讨论 solid–fluid transition：区分静态刚性、动力学解堵、黏弹性交叉与屈服，推导松弛模量、频率响应和蠕变，比较 SPV 相图、流动曲线及 Fielding 本构模型，并分别讨论线虫的组织响应与细胞质微流变。斑马鱼的组织转变与线虫的直接证据分别注明。新增解析图核查蠕变与松弛在 Laplace 域中的互易关系、长期极限与屈服前后响应。

第二篇补充表面与胞质极性耦合、四细胞几何导出的 T1 路径、双层相场穿插、rosette 受力与连接稳定、非均匀前缘及黏附成熟、负载共享失效与三维界面模型。现含 14 篇文献、11 幅论文原图和 5 幅模型图；新增三幅模型图含实际轮廓与形态序列。空间模型参数与数值核查见 figures/epithelial_spatial_model_parameters.json；相场面积误差、四细胞总面积与能垒、前缘相遇及再打开条件在生成脚本中核查。所有空间算例均未拟合胚胎数据。

第七篇从 actin 丝、马达与交联出发，讨论主动凝胶、屏蔽长度、分子周转与输运、Rho 可激发性及收缩失稳、流动排列、手性转矩、连接耦合、PAR 输运与胞质回流。含 33 篇文献、23 幅论文图、8 幅计算图与 3 段远程原始影片。新增论文图均完整保留并分别注明许可；2026 年周转与微流变工作明确标为预印本。计算脚本核查谱空间力平衡、标记群质量衰减、休息态稳定性、取向解析平均与椭圆边界无穿流条件，参数见 figures/subcellular_model_parameters.json。

第七篇以网络材料、流动与信号、细胞与邻居、全身形态、机制检验五条主线组织正文，仅在全身形态部分保留三处小节；沿内移、表皮插入与包裹、早期伸长、跨组织传力、周期重塑和整体弯扭展开网络的形态发生作用；明确区分合子皮层流与伸长期的方向性承载。新增计算核查局部活动到远场的线性力平衡、封闭圆柱体积约束、周期加载后可恢复与参考应变的分离，参数见 figures/morphogenesis_network_parameters.json。新增论文原图为 Zhang 等（2023）Fig. 1 与 Dai 和 Ben Amar（2023）Fig. 3、Fig. 7；其余阶段图复用系列现有原图。
