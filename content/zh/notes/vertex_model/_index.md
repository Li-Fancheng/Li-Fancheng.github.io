---
title: "Vertex Model：从细胞几何到组织力学"
date: 2026-10-18
summary: "以统计物理为主线的 7 篇理论笔记：局部细胞几何、力学相互作用与拓扑约束，怎样产生组织尺度的弹性、塑性、流动和形态变化。"
view: article-grid
layout: series
article_class: vertex-model-note
math: true
---

**局部的细胞几何、力学相互作用和拓扑约束，怎样产生整个组织的弹性、刚性、塑性与流动？** 这个专题围绕这一问题展开，从最小模型的推导一直走到生物形态发生中的机制判断。

全系列分为五部分，共 **7 篇讲义式 note**。模型基础合为一篇；刚性与动力学各用两篇；几何推广、生物应用各用一篇。正文从实验或具体力学问题出发，逐步解释模型假设、完整计算、论文结果和适用条件；较长的计算放在各篇技术补充中。

全系列包含 **34 张论文原图与 25 张自制教学图**。论文图附原文图号、来源链接和读图说明；自制图用于拆解几何、推导响应或复现明确给出的简化计算。图片可以点击查看原尺寸。

另配 **8 段论文原始视频**：7 段显微成像或成像叠加分析，1 段明确标注的模型动画。视频放在相应推导与实验旁边，附“先看哪里”和“怎样读物理”的说明。

| 部分 | Note | 要回答的问题 |
| --- | --- | --- |
| I · 模型建立 | [01 · 几何、能量与力学平衡](/notes/vertex_model/foundations/) | 为什么用多边形铺砌？面积、周长与黏附怎样进入受力？ |
| II · 刚性的几何起源 | [02 · 线性响应与非仿射弹性](/notes/vertex_model/linear_response/) | 怎样从顶点 Hessian 算出组织剪切模量？ |
| II · 刚性的几何起源 | [03 · 几何不相容、预应力与刚性转变](/notes/vertex_model/rigidity/) | 欠约束的网络为什么仍然能成为固体？ |
| III · 从固体到流动 | [04 · T1 重排、能量势垒与塑性](/notes/vertex_model/rearrangements/) | 邻接网络怎样改变，什么时候产生不可逆形变？ |
| III · 从固体到流动 | [05 · 组织流变与主动动力学](/notes/vertex_model/rheology_activity/) | 局部重排怎样变成应力松弛和持续流动？ |
| IV · 几何推广 | [06 · 各向异性、曲面与三维组织](/notes/vertex_model/geometry_3d/) | 换方向、加曲率或增加细胞高度，会改变哪些机制？ |
| V · 生物应用与边界 | [07 · 形态发生、实验检验与开放问题](/notes/vertex_model/morphogenesis/) | 模型复现形状以后，还需要哪些证据才能解释机制？ |

## 斑马鱼贯穿线：从局部界面到发育功能

斑马鱼在这里不仅是最后的应用例子。它贯穿模型建立、刚性、流变和几何推广，用实验检验“哪些局部变化真正决定组织状态”。同一物种的不同发育时期也需要不同模型。

| 从哪个问题读起 | 对应讲义 | 理论与实验的连接 |
| --- | --- | --- |
| 黏附为什么有时使组织变软、有时使它变硬？ | [01 · 界面能与接触角](/notes/vertex_model/foundations/#adhesion-interfaces) | 完全汇合的目标周长，与允许新接触的组织具有不同控制方向 |
| 接触网络仍然连通，为什么材料响应却会明显改变？ | [03 · 胚盘刚性渗流](/notes/vertex_model/rigidity/) | 独立约束、刚性簇、接触寿命与有限时间黏度 |
| 活细胞的涨落怎样参与流动，又由谁选择方向？ | [05 · 主动流变与液滴测量](/notes/vertex_model/rheology_activity/) | 自发涨落与受控响应不同；总事件活性与净塑性应变不同 |
| 细胞接触怎样反过来改变信号的空间范围？ | [06 · 细胞外通道与运输](/notes/vertex_model/geometry_3d/#porosity-transport) | 两张网络、孔隙中的物质守恒与 Nodal–Lefty 反馈 |
| 怎样判断材料变化真的解释了发育过程？ | [07 · 斑马鱼机制检验](/notes/vertex_model/morphogenesis/#zebrafish) | MPZ–PSM 材料梯度、胚盘铺展，以及独立调节密度和界面作用 |

论文线索从 Mongera（2018）、Petridou（2019、2021）与 Kim（2021），延伸到 Rustarazo-Calvo 和 Autorino 的两项 2026 年正式研究。每篇正文区分原论文的测量、模型结论、教学推导与进一步可检验的预测。

## 论文视频：把静态图还原成过程 {#paper-videos}

视频可在正文中点击播放，并使用播放器的全屏、进度条和速度菜单查看细节。延时成像的真实时间以画面时间戳和原图例为准；几秒的播放可能对应数十分钟发育。每段都保留原文及原始视频链接。

| 过程 | 观看重点 | 正文位置 |
| --- | --- | --- |
| 线虫表皮皮层切割 | 切口两侧的回缩与后续响应 | [01 · 激光切割](/notes/vertex_model/foundations/#video-vuong-laser) |
| 果蝇胚带交界收缩 | 两个端点是否同步，顶点相对轮廓怎样滑移 | [04 · 顶点滑移](/notes/vertex_model/rearrangements/#video-vanderleest-sliding) |
| 斑马鱼胚盘刚性图谱 | 黑色接触图与红色最大刚性簇的区别 | [03 · 成像与网络分析](/notes/vertex_model/rigidity/#video-autorino-rigidity) |
| 斑马鱼 A-PSM 细胞起伏 | 局部活动与长期接触更新是否相同 | [05 · 主动起伏](/notes/vertex_model/rheology_activity/#video-mongera-jiggling) |
| 主动泡沫形变后的慢过程 | 后续重排如何继续发生；这是模型动画 | [05 · 主动松弛](/notes/vertex_model/rheology_activity/#video-kim-active-relaxation) |
| 斑马鱼间隙与信号响应 | Smad2 核定位、间隙开闭与响应范围 | [06 · 孔隙与运输](/notes/vertex_model/geometry_3d/#video-autorino-nodal) |
| 斑马鱼胚盘铺展 | 发育动态与固定构型的三维旋转展示 | [07 · 全胚空间过程](/notes/vertex_model/morphogenesis/#video-rustarazo-epiboly) |
| 斑马鱼不同材料条件 | 同一条接触与周围液体空间如何重组 | [07 · 接触动力学](/notes/vertex_model/morphogenesis/#video-rustarazo-contacts) |

## 阅读时始终区分的三个量

- **刚性（rigidity）**：在指定边界、构型和松弛规则下，微小准静态应变是否需要有限的能量代价，例如剪切模量是否非零。
- **流动性（fluidity）**：在给定驱动和观察时间内，组织能否持续变形、松弛应力或发生不可逆重排。
- **几何指标**：例如实测的细胞形状指数 $p_i=P_i/\sqrt{A_i}$，可以与力学相关，却不是完整的力学状态。

目标形状指数 $p_0=P_0/\sqrt{A_0}$ 是模型参数，实测 $p_i$ 是构型的结果。经典二维无序 Vertex Model 中常见的 $p_0^*\simeq3.81$ 依赖模型与制备条件，不是所有组织共用的流动阈值；其来源与限制放在 [Note 03](/notes/vertex_model/rigidity/) 推导。

建议先读 01–05，建立“能量 → 刚性 → 势垒 → 重排 → 流变”的联系，再读 06–07。需要线性代数、微积分和基础统计物理；文中会定义组织力学所需的新量。读图时注意区分真实成像、论文数值模拟、解析计算与教学示意，各图的图注会说明其性质。
