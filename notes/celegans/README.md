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

100 幅论文原图保留面板和标注，来源、权利信息、校验值与逐篇图号见 `embryonic_morphogenesis/media_sources.json`；108 条文献元数据见 `literature_sources.json`。显示宽度与原始尺寸见 `data/celegans_figures.json`，网页保留原图链接。

26 幅当前展示的关系示意与计算图使用文中明确给出的几何关系或无量纲模型。生成脚本为 `figures/generate_figures.py`、`figures/generate_model_figures.py` 与 `figures/generate_extended_models.py`、`figures/generate_active_packing_figures.py`、`figures/generate_worm_mechanics_figures.py`、`figures/generate_solid_fluid_figures.py`、`figures/generate_epithelial_spatial_models.py`、`figures/generate_subcellular_models.py`、`figures/generate_morphogenesis_network_models.py` 与 `figures/generate_mechanism_synthesis.py`，需要 NumPy、Matplotlib；扩展模型、表皮空间、亚细胞与机制综合脚本另需 SciPy。论文原图不经过这些脚本修改。

各阶段开头嵌入 7 段原始显微影片，第二篇分别展示背侧插入、腹侧包裹和腹侧裂隙封闭；原有激光切割影片保留。影片直接使用研究机构或期刊的远程 MP4，保留播放控制、采样间隔和原始来源，记录见 media_sources.json 中的 embedded_stage_movies。竖向影片限制显示宽度，避免在长篇正文中过度放大。

总览补充胚胎结构、细胞谱系、荧光与 DIC 成像、RNAi、恢复实验和模型近似的说明。七篇正文在术语首次出现处解释含义，并展开从影像测量到模型变量的对应关系；具体假设和实验局限保留在相关段落。

模型部分按尺度与机制比较多条路线，新增三维细胞形状、相场、主动凝胶、连接受力失效、三维插入、非线性纤维、孔弹性、周期更新、Cosserat 杆与参数信息分析。新增图均有明确方程和无量纲参数；候选框架与已发表线虫模型、预印本结果分别说明。新增两幅论文图展示相场与皮层反馈模型。

第六篇以文献综述为主：前四节联系线虫的排列、内移、手性流与伸长，后半篇比较已发表的独立皮层、主动顶点、接触几何和分区身体模型。第8节的自编候选计算撤下，改为爪蟾插入与力传播、线虫骨架重塑和三维伸长的论文综述；原始模型图与实验图完整引用，并复用Tetley等的远程模型动画。第7节压缩重复判据，第9节综合文献建立的机制，第10节缩为向列与力网络拓扑的简短扩展。四细胞排列补充P2–EMS接触与分裂方向及内胚层诱导的联系。

第二篇以已发表理论与实验为主线，展开 Brodland 的突起驱动、Weng 的独立皮层与爬行—收缩协作、Tetley 的身份依赖重排、Sknepnek 的主动 T1、Hutson 与 Layton 的组织受力、Wang 的细胞内网络及 Tah 的闭合刚性模型；每组分别讨论线虫可检验的对应。含 26 篇文献、27 幅论文原图与 7 幅计算图、3 段线虫阶段影片及 1 段已发表模型动画。极性、顶底耦合、重排能垒和闭合计算与相应机制合写，参数和生成脚本保留。新增图完整保留，来源、版本与许可记入媒体清单；检索与证据分工见 epithelial_theory_review_20261025.json。

第七篇按五条机制主线组织：网络产力与耗散、反应和主动汇聚、取向与极性和手性、皮层向胞质及连接传力，以及网络在胚胎形态变化中的作用。线虫实验与对应理论相邻，其他体系的模型在相关段落说明原研究对象。含33篇文献、29幅论文图、5幅辅助原理示意与3段远程原始影片。新增Michaux（2018）Fig.10的实验约束反馈模型、Slabodnick（2023）Fig.5的配对运动测量流程，并复用Fig.6的zyxin扰动、Lardennois与Yang的原始模型和分子交换图。波传播、均匀圆柱及标量周期塑性演示撤下；已有屏蔽与网络传力示意紧邻相应实验，输运、仿射取向和解析回流图缩为辅助插图。完整组织重排与伸长模型衔接第六篇。

2026-10-25 对其余六篇作理论综述复核：01 展开多细胞双稳态、噪声选择与分裂打断松弛；03 比较连续体、预应力与方向更新模型；04 以已发表参考状态、actin 重塑及连接周转模型重排主线；05 展开三维主动材料、能量转换假设与切割验证的力反演；06 补充持续驱动对软模态的选择及发育中的参考周长更新；07 展开产力与耗散的共同周转依赖，并以真实三维内部流场计算补充解析回流。通用推导与简化算例已按主题融入正文；07 保持原有五条主线。新增 14 幅论文原图，Miao 预印本两幅图从公开 PDF 中无损提取嵌入的原始 JPEG；所有原有图和阶段影片保留。逐篇复核记录见 series_theory_review_20261025.json。

按连续阅读重整七篇正文：14 个折叠块全部取消，方程、模型解释与计算图回到对应的机制讨论；重写各篇开头、重复小结及部分图题，删去撰写过程说明和重复提醒。显示方程、论文图、文献表和阶段影片保持完整，局部图号与媒体位置清单同步更新。

除按机制整合的第七篇外，各篇按具体主题划分多个正文大节：前几节在标题中标明“线虫研究”，讨论直接针对线虫的模型与实验；后几节标明“理论借鉴”，讨论其他体系、通用理论与候选计算。相邻的短主题在同一大节中展开，避免目录过碎。研究阶段、预印本状态、模型的原研究体系及其与线虫的对应保留在正文中，公式、论文图、文献表、视频和旧主题链接保留。

编号沿用 Vertex 系列的层级：篇名为 01、02 等；正文大节为 1.、2.，小节为 1.1、1.2，更下一层为 1.1.1。正文主题大节连续编号，小节随所属大节编号；旧主题锚点保留，参考文献不编号。

2026-10-25 按机制综合意见重整：总览以有效形变、方向选择和形态保持开篇，01–05 为发育阅读路径，06–07 为横向专题；基础成像与遗传实验移至阅读指南，以模型—观测—竞争解释表收束。七篇开头分别给出当前判断与未决机制。02 末尾补充接触成熟与加载时序、闭合留下的内部状态，03 衔接应力、刚度与历史的辨识；04 区分 actin 与连接相关记忆的证据，并比较机械并联、串联及更新时序。06 保留邻接拓扑于主线，将主动向列和持续同调移至末尾扩展阅读。07 补充网络受力到组织材料参数的粗粒化关系；01 与05 分别补充内移及三维弯曲的竞争机制比较。

新增两幅候选计算图及可复现脚本、参数与校验：闭合后的相同外形、净应力和初始速率不唯一决定内部状态；线性负载增加与独立成熟下，接触稳定性可由明确的容量边界检验。未增加论文文献或替换论文原图，未将候选参考状态等同于已发表连接周转模型。复核记录见 `embryonic_morphogenesis/mechanism_synthesis_20261025.json`。

2026-10-26 第六篇按文献综述定位重写，未加入新的自建模型或计算验证。当前该篇包含29幅图，其中26幅为论文原图，3幅为约束、T1和标准材料关系的辅助示意。撤下的旧算例文件保留于图像库和历史记录，未删除或改写原论文图片。旧主题锚点、论文来源与远程影片保留，逐图位置重新同步。复核记录见 `embryonic_morphogenesis/active_matter_literature_revision_20261026.json`。

2026-10-26 第七篇按文献综述定位重排，图号1–34连续，论文图占29幅；未新增自建计算。新取得的两幅论文图保留原始面板，原有98幅论文图校验值不变，已撤下的三幅自制算例作为历史素材保留。旧主题锚点、三段影片与文献表保留。复核记录见 `embryonic_morphogenesis/subcellular_literature_revision_20261026.json`。
