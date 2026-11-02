---
title: "从生物现象到可解问题：涨落、几何与集体行为"
date: 2026-11-02
summary: "追踪生物问题怎样转化为可计算的随机过程、几何力学、统计推断和优化模型。"
layout: series
show_article_cards: false
article_class: neural-theory-note
math: true
tags:
  - 定量生物学
  - 统计物理
  - 文献阅读
---

生物现象往往同时涉及分子、细胞和环境。理论建模的第一步不是把所有细节写进方程，而是选择一个问题：什么变量足以表达现象？去掉哪些机制后，仍能得到可以检验的预测？

本专题沿着八个案例，追踪从观察到可计算对象的转译。每篇区分原论文的生物背景、便于推导的最小模型与可做的对照实验；推导和实验方案是理解机制的工具，不等于对原论文全部结果的复现。静态拟合、信息上限和最优性假设都不能单独证明真实机制。

## 八个问题，八种转译

| 尺度 | 问题与最小模型 | 理论入口 |
| --- | --- | --- |
| 亚细胞 | [细胞怎样用随机涨落产生定向力？](/notes/solvable_biology/brownian_ratchet/)：布朗棘轮：把力学运动和化学反应接起来 | 受限布朗运动、化学反应 |
| 单细胞 | [相同基因的细胞为什么表达不同？](/notes/solvable_biology/transcription_bursts/)：转录爆发：从平均值转向分布与时间 | 随机开关、分子数分布 |
| 单细胞与谱系 | [随机生长和分裂下，细胞大小为何稳定？](/notes/solvable_biology/cell_size/)：Adder：细胞周期问题变成代际随机映射 | 随机映射、稳态与相关 |
| 组织 | [组织密度不变，为什么仍能变得容易流动？](/notes/solvable_biology/tissue_rigidity/)：顶点模型：把重排写成几何能垒 | 几何约束、刚性和能垒 |
| 胚胎 | [细胞怎样从带噪声的信号推断位置？](/notes/solvable_biology/positional_information/)：位置编码：把发育问题变成统计解码 | 贝叶斯推断、信息论 |
| 个体行为 | [动物什么时候应该离开食物区域？](/notes/solvable_biology/foraging/)：边际价值定理：先选定目标函数 | 收益率优化、最优停止 |
| 种群演化 | [大种群为什么会被随机漂变分成色块？](/notes/solvable_biology/range_expansion/)：扩张前沿：局部抽样放大成宏观扇区 | 随机生长、谱系与界面 |
| 生态群落 | [为什么少数物种常见，而许多物种稀有？](/notes/solvable_biology/neutral_ecology/)：中性生态：用随机过程建立可反驳基线 | 中性过程、抽样与模型比较 |

## 怎样阅读？

可以从[细胞大小控制](/notes/solvable_biology/cell_size/)入手：一个随机递推就能算出稳态方差与代际相关。接着读[组织刚性](/notes/solvable_biology/tissue_rigidity/)和[扩张前沿](/notes/solvable_biology/range_expansion/)，比较随机映射、能量景观与非平衡增长三种建模方式。

每篇最后都提出可以区分解释的测量或数值对照。它们是待执行的实验设计，文中的玩具推导及预测都注明了假设；不能把建议的计算当成已经得到的数值结果。

相关站内笔记：[神经网络中的问题重构与物理近似](/notes/statistical_physics_of_learning/)、[定量生物学课程](/notes/qbio_course/)和[形态发生系列](/notes/morphogenesis/)。
