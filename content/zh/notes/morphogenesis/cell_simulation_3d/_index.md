---
url: "/notes/cell_simulation_3d/"
weight: 3
title: "三维细胞模拟：从相场、网格到 Vertex"
date: 2026-10-25
summary: "九篇长文，结合 21 张论文图、18 幅计算图与 53 篇文献，从模型思想的来源出发，连贯讨论三维细胞的接触、材料、重排与形态发生，并分析各类模型的适用范围和局限。"
layout: cell-cycle
article_class: cell-simulation-3d-note
math: true
toc: true
---

两个细胞靠在一起，接触面为什么会变大？一层上皮怎样弯曲，一团细胞怎样在卵壳里重新排列，一个腔又怎样靠水通量长出来？三维细胞模拟试图把这些形状变化与具体的受力、材料和生物过程联系起来。

九篇笔记将论文中的三维形态与具体力学推导放在一起，逐步建立相场、表面网格、三维 Vertex、Cellular Potts、粒子与连续材料的方程。每篇以 5—6 个大节组织连续论述：开头说明模型思想的来源与变量选择，正文结合推导和论文案例，结尾讨论近似在什么条件下成立、何时需要扩展。细胞接触、形变和运动分别约束界面、材料与耗散；胚胎、视杯和类器官则展示这些作用怎样配合形成组织。

## 1. 专题目录：五部分、九篇笔记 {#chapters}

| 部分 | 笔记 | 主要问题 |
| --- | --- | --- |
| I · 共同框架 | [01 · 几何、能量与模型选择](/notes/cell_simulation_3d/foundations/) | 三维接触、压力、张力和耗散怎样进入同一框架？ |
| II · 细胞的离散表示 | [02 · 多相场：从弥散界面到细胞形变](/notes/cell_simulation_3d/phase_field/) | 不追踪膜节点，怎样演化细胞轮廓、保持体积与身份？ |
| II · 细胞的离散表示 | [03 · 表面网格：皮层、接触与重网格](/notes/cell_simulation_3d/surface_mesh/) | 三角形怎样承载张力、弹性与弯曲，如何避免穿透？ |
| II · 细胞的离散表示 | [04 · 三维 Vertex：共享多面体与拓扑重排](/notes/cell_simulation_3d/vertex_3d/) | 从顶基底上皮到体组织，怎样维护几何与邻接？ |
| II · 细胞的离散表示 | [05 · 格点与 Cellular Potts](/notes/cell_simulation_3d/lattice_potts/) | 身份复制如何形成运动，格点与随机时钟带来什么？ |
| II · 细胞的离散表示 | [06 · 中心粒子与亚细胞元](/notes/cell_simulation_3d/particles/) | 一个点何时够用，一组元又怎样产生有效材料？ |
| III · 胞内与环境 | [07 · 有限元与流固耦合](/notes/cell_simulation_3d/continuum_fluids/) | 何时需要体应变、流体、压力和渗透？ |
| IV · 生物过程 | [08 · 三维形态发生](/notes/cell_simulation_3d/morphogenesis/) | 分裂、极性、腔体、细胞核、基质和信号如何耦合？ |
| V · 定量检验 | [09 · 参数反演、数值验证与可复现模拟](/notes/cell_simulation_3d/validation/) | 哪些参数能识别，怎样证明结果不依赖数值选择？ |

## 2. 三维模型中的不同网格 {#mesh-types}

- **表面三角网格**：一层面片包住细胞，直接解析皮层与接触，见 [03](/notes/cell_simulation_3d/surface_mesh/)。
- **共享多面体网络**：多个细胞共用面、边和顶点，保存汇合组织的拓扑，见 [04](/notes/cell_simulation_3d/vertex_3d/)。
- **体素格点**：空间网格记录细胞身份，或用来求解相场，见 [02](/notes/cell_simulation_3d/phase_field/) 与 [05](/notes/cell_simulation_3d/lattice_potts/)。
- **体有限元网格**：填充细胞内部、计算材料应变和应力，见 [07](/notes/cell_simulation_3d/continuum_fluids/)。

三角形边通常是数值连接；细胞交线是生物几何。表面重网格与细胞换邻居是不同事件。曲面上的二维细胞网、有厚度上皮和完整三维体组织，也需要分别建模。

## 3. 研究对象与模型路线 {#reading-routes}

| 想研究什么 | 建议路线 | 首先核对 |
| --- | --- | --- |
| 接触面积、界面曲率与自由变形 | 01 → 02 / 03 → 09 | 体积、有效张力、界面与接触收敛 |
| 汇合组织重排与顶基底形状 | 01 → 04 → 08 → 09 | 厚度自由度、极性与重连规则 |
| 细胞分选、类型差异与随机运动 | 01 → 05 / 06 → 09 | 接触代价、涨落、连通性与时钟 |
| 大细胞群与营养、信号的反馈 | 01 → 06 → 08 → 09 | 源汇守恒、力学与周期时间尺度 |
| 胞内应变、细胞核和受限迁移 | 03 → 06 / 07 → 08 | 本构、核—皮层耦合与耗散 |
| 腔体膨胀、渗透与流体 | 01 → 07 → 08 → 09 | 压力、体积、水和溶质的记账 |

不同路线可以通过共同的单细胞、双细胞与小组织基准互相比较；更高分辨率是否必要，取决于新增形状或内部变量能否解释独立观测。

## 4. 论文图与三维形态 {#paper-figures}

正文收录 21 张论文图，分布在九篇讲解中。图旁说明了面板、颜色和视角的读法，并区分实验成像、分割重建、模拟形态和力学推断。可以先沿这条路线看形状，再读相应的方程。

| 想先看什么 | 论文图与读图位置 | 回到模型要问什么 |
| --- | --- | --- |
| 细胞铺展与高度改变 | [黏附条件下的三视图](/notes/cell_simulation_3d/foundations/#paper-adhesion) | 接触、体积和皮层张力怎样共同约束外形？ |
| 卵壳里的卵裂与挤压 | [相场与膜成像](/notes/cell_simulation_3d/phase_field/#paper-embryo) · [多视角接触检验](/notes/cell_simulation_3d/morphogenesis/#paper-embryo-views) | 哪些是分裂输入，哪些是力学产生的排列？ |
| 细胞从一层里挤出 | [三维表面与速度截面](/notes/cell_simulation_3d/phase_field/#paper-extrusion) | 平面内面积变化能否区分变形与离层？ |
| 片、管、囊与细长细胞 | [表面网格组成的组织](/notes/cell_simulation_3d/surface_mesh/#paper-mesh) · [细胞与核](/notes/cell_simulation_3d/surface_mesh/#paper-nuclei) | 表面自由度和核分别承担什么形变？ |
| 沿上皮厚度改变的邻居 | [scutoid 的成像与重建](/notes/cell_simulation_3d/vertex_3d/#paper-scutoid) · [三维重连](/notes/cell_simulation_3d/vertex_3d/#paper-reconnection) | 空间邻接改变与时间上的重排如何区分？ |
| 球状团侵袭与内部坏死 | [三维 CPM 侵袭形态](/notes/cell_simulation_3d/lattice_potts/#paper-invasion) · [PhysiCell 内部结构](/notes/cell_simulation_3d/particles/#paper-spheroid) | 边界粗糙、细胞散出和内部状态怎样分别检验？ |
| 黏性松弛与周围基质 | [分裂后的细胞团](/notes/cell_simulation_3d/continuum_fluids/#paper-active-surfaces) · [三维基质位移](/notes/cell_simulation_3d/continuum_fluids/#paper-tractions) | 时间尺度和材料选择怎样改变形态与力的解释？ |
| 内陷、出芽与器官形状 | [上皮壳形态图](/notes/cell_simulation_3d/morphogenesis/#paper-shell-shapes) · [视杯](/notes/cell_simulation_3d/morphogenesis/#paper-optic-cup) · [肠隐窝](/notes/cell_simulation_3d/morphogenesis/#paper-crypt) | 局部收缩、重排与腔体变化怎样配合？ |
| 从形状推断受力 | [foambryo](/notes/cell_simulation_3d/validation/#paper-foambryo) · [FIDES](/notes/cell_simulation_3d/validation/#paper-fides) | 图中的力学颜色经过了哪些假设？ |

论文图保留来源、图号和权利说明，点击可查看完整尺寸。配图用来理解现象、几何与证据；单张形态相似还需要接触、时间序列或干预来检验机制。

## 5. 教学计算与数值基准 {#figures}

正文有 18 幅可复现图，每篇两幅。新增计算包括双细胞接触角、相场薄界面动力学、等面积剪切、形状指数、体素能量换算、营养与持续运动、黏弹频率响应、多孔弹性长度标度、屈曲和曲率反演误差。

这 18 幅图是自制几何、解析曲线或明确方程的低维数值计算；图注分别说明计算内容与参数。它们用于解释模型与验证思路，不是实验数据，也不是完整三维组织模拟的结果。论文原图另按上述读图路线穿插在正文中。生成脚本、原图来源与校验方法见仓库中的 `notes/cell_simulation_3d/README.md`。

## 6. 综述与代表性研究 {#literature}

| 文献 | 适合解决的阅读问题 |
| --- | --- |
| [Van Liedekerke 等，2015](https://doi.org/10.1007/s40571-015-0082-3) | 格点、中心、可变形细胞、Vertex 与混合模型保留了哪些不同物理过程？ |
| [Fletcher、Cooper 与 Baker，2017](https://doi.org/10.1098/rstb.2015.0519) | 上皮形态发生中，不同细胞表示各自保留哪些力学过程？ |
| [Monfared、Ardaševa 与 Doostmohammadi，2026](https://doi.org/10.1146/annurev-conmatphys-060625-061354) | 多相场如何从单细胞形变发展到集体迁移与三维组织？ |
| [Okuda、Inoue 与 Adachi，2015](https://doi.org/10.2142/biophysico.12.0_13) | 三维 Vertex 怎样组织收缩、生长、重排与分裂？ |
| [Voß-Böhme，2012](https://doi.org/10.1371/journal.pone.0042852) | CPM 的随机更新与连续、多尺度解释之间有哪些条件？ |
| [Neagu 等，2026](https://doi.org/10.3390/cells15020177) | 类器官研究已经采用哪些细胞模型，分别用于什么问题？ |

专题共引用 53 篇论文、综述与方法分析。各篇篇末列出实际引用的完整题名、期刊信息和 DOI；正文在相关推导与案例旁给出链接。经典工作解释组织接触、弥散界面与格点模型的来源，综述建立不同表示的发展脉络，具体机制再结合原始研究的输入、观测与干预讨论。

## 7. 相关专题 {#related-topics}

[Vertex Model 专题](/notes/vertex_model/)详细讨论几何刚性、线性响应、T1、流变与主动涨落；这里将重点放在**三维表示之间的差别，以及从物理到可计算程序的选择**。关于从二维推广到三维，可先读[各向异性、曲面与三维组织](/notes/vertex_model/geometry_3d/)。

生物过程可与[线虫胚胎形态发生](/notes/celegans/embryonic_morphogenesis/)对照。具体软件实现及原始论文见 [09 · 软件实现与原始文献](/notes/cell_simulation_3d/validation/#software)。
