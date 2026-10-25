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

98 幅论文原图保留面板和标注，来源、权利信息、校验值与逐篇图号见 `embryonic_morphogenesis/media_sources.json`；105 条文献元数据见 `literature_sources.json`。显示宽度与原始尺寸见 `data/celegans_figures.json`，网页保留原图链接。

33 幅计算图使用文中明确给出的几何关系或无量纲模型。生成脚本为 `figures/generate_figures.py`、`figures/generate_model_figures.py` 与 `figures/generate_extended_models.py`、`figures/generate_active_packing_figures.py`、`figures/generate_worm_mechanics_figures.py`、`figures/generate_solid_fluid_figures.py`、`figures/generate_epithelial_spatial_models.py`、`figures/generate_subcellular_models.py`、`figures/generate_morphogenesis_network_models.py`，需要 NumPy、Matplotlib；扩展模型、表皮空间与亚细胞模型脚本另需 SciPy。论文原图不经过这些脚本修改。

各阶段开头嵌入 7 段原始显微影片，第二篇分别展示背侧插入、腹侧包裹和腹侧裂隙封闭；原有激光切割影片保留。影片直接使用研究机构或期刊的远程 MP4，保留播放控制、采样间隔和原始来源，记录见 media_sources.json 中的 embedded_stage_movies。竖向影片限制显示宽度，避免在长篇正文中过度放大。

总览补充胚胎结构、细胞谱系、荧光与 DIC 成像、RNAi、恢复实验和模型近似的说明。七篇正文在术语首次出现处解释含义，并展开从影像测量到模型变量的对应关系；具体假设和实验局限保留在相关段落。

模型部分按尺度与机制比较多条路线，新增三维细胞形状、相场、主动凝胶、连接受力失效、三维插入、非线性纤维、孔弹性、周期更新、Cosserat 杆与参数信息分析。新增图均有明确方程和无量纲参数；候选框架与已发表线虫模型、预印本结果分别说明。新增两幅论文图展示相场与皮层反馈模型。

第六篇比较受限软颗粒、主动 Brownian 粒子、顶点／Voronoi、三维细胞、主动皮层与向列模型。正文进一步讨论线虫四细胞排列、界面张力、内移的堆积与 clutch、手性皮层流、背侧定向插入、腹侧包裹与连接成熟、actin 各向异性、周期驱动与材料记忆、有效力反演；含 44 篇文献、17 幅论文图、9 幅计算或结构图，以及线虫皮层流和培养细胞取向运动的 2 段远程影片。新增算例核查双势阱失稳、邻居松弛的 Schur 补、clutch 力平衡、固定邻接的仿射形变和周期加载后的参考应变。

第六篇第二部分集中讨论 solid–fluid transition：区分静态刚性、动力学解堵、黏弹性交叉与屈服，推导松弛模量、频率响应和蠕变，比较 SPV 相图、流动曲线及 Fielding 本构模型，并分别讨论线虫的组织响应与细胞质微流变。斑马鱼的组织转变与线虫的直接证据分别注明。新增解析图核查蠕变与松弛在 Laplace 域中的互易关系、长期极限与屈服前后响应。

第二篇以已发表理论与实验为主线，展开 Brodland 的突起驱动、Weng 的独立皮层与爬行—收缩协作、Tetley 的身份依赖重排、Sknepnek 的主动 T1、Hutson 与 Layton 的组织受力、Wang 的细胞内网络及 Tah 的闭合刚性模型；每组分别讨论线虫可检验的对应。含 26 篇文献、27 幅论文原图与 5 幅原有计算图、3 段线虫阶段影片及 1 段已发表模型动画。极性、顶底耦合、重排能垒和闭合计算与相应机制合写，参数和生成脚本保留。新增图完整保留，来源、版本与许可记入媒体清单；检索与证据分工见 epithelial_theory_review_20261025.json。

第七篇从 actin 丝、马达与交联出发，讨论主动凝胶、屏蔽长度、分子周转与输运、Rho 可激发性及收缩失稳、流动排列、手性转矩、连接耦合、PAR 输运与胞质回流。含 33 篇文献、26 幅论文图、8 幅计算图与 3 段远程原始影片。新增论文图均完整保留并分别注明许可；2026 年周转与微流变工作明确标为预印本。计算脚本核查谱空间力平衡、标记群质量衰减、休息态稳定性、取向解析平均与椭圆边界无穿流条件，参数见 figures/subcellular_model_parameters.json。

第七篇按直接研究与可借鉴理论分成两部分，保留五条主要主题；沿内移、表皮插入与包裹、早期伸长、跨组织传力、周期重塑和整体弯扭展开网络的形态发生作用；明确区分合子皮层流与伸长期的方向性承载。新增计算核查局部活动到远场的线性力平衡、封闭圆柱体积约束、周期加载后可恢复与参考应变的分离，参数见 figures/morphogenesis_network_parameters.json。新增论文原图为 Zhang 等（2023）Fig. 1 与 Dai 和 Ben Amar（2024）Fig. 3、Fig. 7；其余阶段图复用系列现有原图。

2026-10-25 对其余六篇作理论综述复核：01 展开多细胞双稳态、噪声选择与分裂打断松弛；03 比较连续体、预应力与方向更新模型；04 以已发表参考状态、actin 重塑及连接周转模型重排主线；05 展开三维主动材料、能量转换假设与切割验证的力反演；06 补充持续驱动对软模态的选择及发育中的参考周长更新；07 展开产力与耗散的共同周转依赖，并以真实三维内部流场计算补充解析回流。通用推导与简化算例已按主题融入正文；07 保持原有五条主线。新增 14 幅论文原图，Miao 预印本两幅图从公开 PDF 中无损提取嵌入的原始 JPEG；所有原有图和阶段影片保留。逐篇复核记录见 series_theory_review_20261025.json。

按连续阅读重整七篇正文：14 个折叠块全部取消，方程、模型解释与计算图回到对应的机制讨论；重写各篇开头、重复小结及部分图题，删去撰写过程说明和重复提醒。显示方程、论文图、文献表和阶段影片保持完整，局部图号与媒体位置清单同步更新。

七篇现按研究对象分成两个正文部分：先集中介绍线虫的直接研究，再讨论其他体系、通用理论与候选计算。第一部分保留研究阶段和预印本状态；第二部分保留原体系的模型、检验及其与线虫的对应。图号、来源位置和旧主题链接同步调整，公式、论文图、文献表与视频保留。

编号沿用 Vertex 系列的层级：篇名为 01、02 等；正文大节为 1.、2.，小节为 1.1、1.2，更下一层为 1.1.1。两部分结构与旧主题锚点保留，参考文献不编号。
