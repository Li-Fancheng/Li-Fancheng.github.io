---
title: "从生物现象到可解问题：涨落、几何与集体行为"
date: 2026-11-02
lastmod: 2026-11-02
summary: "通过八篇深入案例，追踪生物观察如何经过变量选择、近似、求解与检验，成为可理解的研究问题。"
layout: series
show_article_cards: false
article_class: neural-theory-note
math: true
tags:
  - 定量生物学
  - 统计物理
  - 文献阅读
---

一个陌生的生物现象，最难的地方有时在计算之前：为什么应沿细胞谱系追踪大小？为什么一个二维菌落可以先降为一维前沿，再降为少数边界？为什么间歇转录能用一个隐藏开关描述，而静态分布仍可能无法识别这个开关？

这里用八篇深入案例重构这一步研究逻辑。先呈现令人困惑的观察，再识别直接相关的自由度；尝试最小解释，实际完成一次求解，检查预测能否回到最初的问题。每次增加模型复杂性，都应有一个原先没有解释的观测作为理由。

各篇按论文证据重构有物理依据的研究路线，不将其当作作者真实发现历史。原论文的实验结果、本文的教学推导和待执行的计算方案在相关段落分别说明。

## 八篇案例，从亚细胞到生态群落

| 尺度与案例 | 跟随的研究路线 | 最值得学习的选择 |
| --- | --- | --- |
| 亚细胞：[布朗棘轮](/notes/solvable_biology/brownian_ratchet/) | 膜与纤维的等待 → 间隙变量 → 快膜消元 → 快反应首次到达 → 化学做功 | 根据谁在限速，选择消去哪一类自由度 |
| 单细胞：[转录爆发](/notes/solvable_biology/transcription_bursts/) | Poisson 基线不足 → 隐藏开关 → 矩闭合 → 完整分布 → 时间辨识 | 用涨落发现隐藏状态，再检验静态等价机制 |
| 细胞谱系：[大小控制](/notes/solvable_biology/cell_size/) | 条件分布识别 Adder → 代际收缩与稳态 → 阈值累积 | 分开观察规律、统计后果和候选机制 |
| 组织：[刚性转变](/notes/solvable_biology/tissue_rigidity/) | 密度解释的困难 → T1 重排 → 几何约束 → 能垒与软模 | 沿实际事件寻找宏观力学的控制变量 |
| 胚胎：[位置编码](/notes/solvable_biology/positional_information/) | 联合表达测量 → 位置后验 → 局部精度 → 全局歧义 → 实际读取 | 先定义任务，再判断局部公式是否足够 |
| 个体行为：[最优采食](/notes/solvable_biology/foraging/) | 重复周期 → 机会成本 → 离开阈值 → 随机捕获 → 信念状态 | 选择长期目标与足以预测未来的状态 |
| 种群演化：[扩张前沿](/notes/solvable_biology/range_expansion/) | 大种群遗传分离 → 畴边界 → 扩散不足 → 粗糙界面 → 径向存活 | 怎样根据观测修正有效描述 |
| 生态群落：[中性生态](/notes/solvable_biology/neutral_ecology/) | 个体等价 → 出生死亡闭合 → 联合稳态 → 丰度谱 → 动态检验 | 让零假设算出一组共享参数的预测 |

## 可以怎样穿起来读？

想看一次近似怎样被观测推动着修正，可以先读[种群扩张](/notes/solvable_biology/range_expansion/)，再读[布朗棘轮](/notes/solvable_biology/brownian_ratchet/)中两种不同的限速极限。两篇都强调，有效描述的选择需要与观测和时间尺度相配。

想看统计怎样连接机制，可以从[细胞大小](/notes/solvable_biology/cell_size/)进入，再对照[转录爆发](/notes/solvable_biology/transcription_bursts/)与[中性生态](/notes/solvable_biology/neutral_ecology/)。三者分别用代际关系、分子计数和群落替代获得可解结构，再追问什么额外数据可以辨识机制。

[组织刚性](/notes/solvable_biology/tissue_rigidity/)和[位置编码](/notes/solvable_biology/positional_information/)呈现两种几何：前者是实际形状的相容性，后者是表达空间中的统计距离。[最优采食](/notes/solvable_biology/foraging/)则将问题转向目标和状态的选择，展示已有优化工具如何在生物假设明确以后接上来。

每篇在结尾回看关键转译。阅读时可以不断追问：换一种自然的简化，哪项观测会无法解释？保留下来的自由度，是否足以把计算结果带回最初的生物问题？

相关站内笔记：[神经网络中的问题重构与物理近似](/notes/statistical_physics_of_learning/)、[定量生物学课程](/notes/qbio_course/)和[形态发生系列](/notes/morphogenesis/)。
