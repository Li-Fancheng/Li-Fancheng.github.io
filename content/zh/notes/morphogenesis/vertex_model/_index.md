---
url: "/notes/vertex_model/"
weight: 1
title: "Vertex Model：从细胞几何到组织力学"
date: 2026-10-18
lastmod: 2026-10-25
summary: "从细胞几何、实验测量与生物过程出发，讨论 Vertex Model 的受力、刚性、重排和组织流动。"
view: article-grid
show_article_cards: false
layout: series
article_class: vertex-model-note
math: true
---

细胞边界的几何与能量决定局部恢复力；有限形变可以触发邻居交换；驱动、耗散与结构更新共同决定组织怎样松弛和流动。主线依次连接这些过程，最后回到果蝇、斑马鱼和培养上皮中的实验检验。

| 部分 | 笔记 | 核心问题 |
| --- | --- | --- |
| I · 几何与刚性 | [01 · 几何、能量与力学平衡](/notes/vertex_model/foundations/) | 面积、周长与黏附怎样进入受力？ |
| I · 几何与刚性 | [02 · 线性响应与非仿射弹性](/notes/vertex_model/linear_response/) | 内部调整怎样改变弹性？ |
| I · 几何与刚性 | [03 · 几何不相容、预应力与刚性转变](/notes/vertex_model/rigidity/) | 几何不相容与预应力怎样产生刚性？ |
| II · 重排与动力学 | [04 · T1 重排、能量势垒与塑性](/notes/vertex_model/rearrangements/) | 换邻居需要跨越怎样的能量路径？ |
| II · 重排与动力学 | [05 · 组织流变与主动动力学](/notes/vertex_model/rheology_activity/) | 局部形变与结构更新怎样形成黏弹性？ |
| II · 重排与动力学 | [06 · 涨落、相关与组织重塑](/notes/vertex_model/fluctuations/) | 涨落的持续时间与空间协同怎样影响重塑？ |
| III · 几何与生物推广 | [07 · 各向异性、曲面与三维组织](/notes/vertex_model/geometry_3d/) | 方向、曲率和厚度怎样改变力学？ |
| IV · 实验与机制判别 | [08 · 形态发生、实验检验与开放问题](/notes/vertex_model/morphogenesis/) | 怎样用实验区分竞争机制？ |
| 专题拓展 | [主动软物质的堆积与拓扑](/notes/vertex_model/active_packing/) | 用颗粒堆积、主动向列等体系比较几何约束与集体运动。 |

经典主线采用汇合二维 Vertex 网络的面积—周长能量。非汇合界面、细胞中心接触图、曲面与三维模型各有自己的自由度与力学假设；拓展篇中的 MIPS、二维熔化和向列缺陷用于对照，相关相图不直接作为 Vertex 刚性判据。

同一周期网络的[残差与弹性](/notes/vertex_model/rigidity/#unified-case)、[交界收缩路径](/notes/vertex_model/rearrangements/#unified-t1)和[驱动响应](/notes/vertex_model/rheology_activity/#unified-dynamics)使用同一组制备与参数定义。

三维表示与算法的进一步比较见[三维细胞模拟：从相场、网格到 Vertex](/notes/cell_simulation_3d/)，包括相场、表面网格、完整三维 Vertex、Cellular Potts、亚细胞元、连续材料与流固耦合。

<details>
<summary>生物案例与论文视频索引</summary>

## 生物案例：果蝇和斑马鱼为主，培养上皮与线虫作比较

果蝇提供了几组适合连续追踪的上皮过程：翅的细胞排列、胚带中的邻居交换、腹沟的顶端收缩，以及卵室的基底收缩。斑马鱼则用于研究胚盘与尾部的材料响应、细胞外液体空间和信号传播。

MDCK 和人支气管上皮的培养实验补充了可控几何、单层流变与运动停滞；线虫用于比较基底突出、边缘闭合和细胞骨架塑性。

| 生物过程 | 主要想弄明白什么 | 对应笔记 |
| --- | --- | --- |
| 果蝇翅上皮的排列与切割 | 一张稳定细胞图怎样与受力联系起来，外部连接起什么作用？ | [01 · 排列与切割](/notes/vertex_model/foundations/#fly-wing-packing)；[02 · 蛹翅回缩](/notes/vertex_model/linear_response/#fly-wing-ablation) |
| 果蝇背板的交界起伏与排列 | 平均张力、单条交界的涨落与交换速率怎样区分？ | [06 · 背板实验](/notes/vertex_model/fluctuations/#fly-notum-fluctuations) |
| 果蝇胚带延长 | 哪些接触先收缩，多细胞交汇怎样解开，方向从哪里来？ | [04 · 交界重塑](/notes/vertex_model/rearrangements/#fly-junction-remodelling)；[08 · 果蝇案例](/notes/vertex_model/morphogenesis/#drosophila) |
| 果蝇腹沟形成 | 为什么一阵阵收缩能累积，顶面缩小怎样走向内陷？ | [05 · 脉冲与稳定](/notes/vertex_model/rheology_activity/#fly-pulsed-constriction)；[07 · 面外运动](/notes/vertex_model/geometry_3d/#fly-ventral-furrow) |
| 果蝇卵室延长 | 收缩发生在基底侧时，为什么必须考虑曲面和基质？ | [07 · 基底驱动](/notes/vertex_model/geometry_3d/#fly-egg-chamber) |
| 斑马鱼胚盘铺展与尾部延长 | 怎样测量材料变化，接触更新与身体形状有什么关系？ | [03 · 刚性簇](/notes/vertex_model/rigidity/)；[05 · 液滴测量](/notes/vertex_model/rheology_activity/)；[08 · 发育比较](/notes/vertex_model/morphogenesis/#zebrafish) |
| 斑马鱼细胞外间隙与信号 | 接触重组怎样改变分子经过的通道与信号响应？ | [07 · 孔隙与运输](/notes/vertex_model/geometry_3d/#porosity-transport) |
| MDCK 培养上皮绕障碍物 | 细胞自身变形与换邻居，怎样检验 Maxwell 黏弹性？ | [05 · 单层流动实验](/notes/vertex_model/rheology_activity/#mdck-maxwell) |
| MDCK 培养上皮的细胞挤出 | 细胞长轴的缺陷、局部应力与挤出事件怎样联系？ | [拓展 · 取向与挤出](/notes/vertex_model/active_packing/#epithelial-defects) |
| 人支气管上皮培养 | 实测形状、集体运动和运动停滞怎样比较？ | [03 · 形状与停滞](/notes/vertex_model/rigidity/#airway-jamming) |
| 线虫表皮形态发生 | 突出、自由边缘和骨架更新为什么需要不同模型？ | [08 · 线虫比较](/notes/vertex_model/morphogenesis/#worm-comparison) |

## 论文视频 {#paper-videos}

视频可在正文中点击播放，并使用播放器的全屏、进度条和速度菜单查看细节。延时成像的真实时间以画面时间戳和原图例为准；几秒的播放可能对应数十分钟发育。

| 过程 | 观看重点 | 正文位置 |
| --- | --- | --- |
| 线虫表皮皮层切割 | 切口两侧的回缩与后续响应 | [01 · 激光切割](/notes/vertex_model/foundations/#video-vuong-laser) |
| 果蝇蛹翅外部连接切割 | 切割位置在哪里，翅边缘相对外壳怎样回缩 | [02 · 外部约束](/notes/vertex_model/linear_response/#video-etournay-cut) |
| 果蝇正常与突变蛹翅 | 整体轮廓、翅铰与翅片怎样逐渐改变 | [08 · 同步发育成像](/notes/vertex_model/morphogenesis/#video-etournay-wing) |
| 果蝇胚带交界收缩 | 两个端点是否同步，顶点相对轮廓怎样滑移 | [04 · 顶点滑移](/notes/vertex_model/rearrangements/#video-vanderleest-sliding) |
| 斑马鱼胚盘刚性图谱 | 黑色接触图与红色最大刚性簇的区别 | [03 · 成像与网络分析](/notes/vertex_model/rigidity/#video-autorino-rigidity) |
| 斑马鱼 A-PSM 细胞起伏 | 局部活动与长期接触更新是否相同 | [05 · 主动起伏](/notes/vertex_model/rheology_activity/#video-mongera-jiggling) |
| 主动泡沫形变后的慢过程 | 后续重排如何继续发生；这是模型动画 | [05 · 主动松弛](/notes/vertex_model/rheology_activity/#video-kim-active-relaxation) |
| 斑马鱼间隙与信号响应 | Smad2 核定位、间隙开闭与响应范围 | [07 · 孔隙与运输](/notes/vertex_model/geometry_3d/#video-autorino-nodal) |
| 斑马鱼胚盘铺展 | 发育动态与固定构型的三维旋转展示 | [08 · 全胚空间过程](/notes/vertex_model/morphogenesis/#video-rustarazo-epiboly) |
| 斑马鱼不同材料条件 | 同一条接触与周围液体空间如何重组 | [08 · 接触动力学](/notes/vertex_model/morphogenesis/#video-rustarazo-contacts) |
| MDCK 上皮缺陷与挤出 | 区分速度、长轴方向与缺陷跟踪，比较事件发生的先后 | [拓展 · 缺陷与细胞挤出](/notes/vertex_model/active_packing/#video-saw-extrusion) |
| 果蝇背板局部换邻居 | 交界的长短变化何时转为真实的接触更新 | [06 · 局部交界动态](/notes/vertex_model/fluctuations/#video-curran-neighbour-exchange) |

</details>

## 刚性、流动性与细胞形状

- **刚性（rigidity）**：轻轻使组织变形，它会不会产生恢复力？要用多大力才能让它变形一点？
- **流动性（fluidity）**：持续受力或等待一段时间后，组织能不能重新排列、继续变形，并释放原来的应力？
- **几何指标**：从图像量出的面积、周长和细长程度。它们可能帮助预测软硬，但还要结合接触、受力和时间来理解。

文中的 $p_0=P_0/\sqrt{A_0}$ 是模型设定的偏好；从真实轮廓量出的 $p_i=P_i/\sqrt{A_i}$ 则是实际形状。它们类似弹簧的原长与当前长度，并不总相等。论文常出现的 3.81 是特定二维模型的目标形状指数，[Note 03](/notes/vertex_model/rigidity/) 讨论它的几何线索与适用范围。
