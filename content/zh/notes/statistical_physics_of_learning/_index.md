---
title: "从神经网络问题到可解理论：问题重构与物理近似"
date: 2026-10-30
lastmod: 2026-11-02
summary: "通过问题转译和自己的最小计算，理解神经网络研究怎样与统计学习、动力系统及统计物理建立联系。"
layout: series
show_article_cards: false
article_class: neural-theory-note
math: true
tags:
  - 统计物理
  - 神经网络理论
  - 文献阅读
---

做理论研究时，困难有时在工具之前：即使知道怎样求一个自由能、怎样分析随机矩阵，也未必知道如何把眼前的神经网络问题变成它们能处理的形式。

这里记录的是中间那一步：原来的问题为什么难？换成什么变量、模型或极限后，哪项计算突然变得可行？得到的答案又能解释原现象的什么？有些转化是精确重写，有些是渐近闭合，也有些只是值得检验的机制模型。把这些关系弄清楚，才知道什么时候可以借用已有理论。

## 问题怎样转过去，工具怎样接上来？

| 案例 | 最值得追踪的转化 | 接上的理论 |
| --- | --- | --- |
| [Scaling laws](/notes/statistical_physics_of_learning/spectral_scaling/) | 网络资源与泛化 → 固定特征的随机设计回归 | 统计估计、核谱与随机矩阵 |
| [深度线性网络](/notes/statistical_physics_of_learning/deep_linear_model/) | 多层训练 → 矩阵分解动力学 → 特殊条件下的独立模式 | SVD 与非线性常微分方程 |
| [NTK 与平均场](/notes/statistical_physics_of_learning/changing_variables/) | 参数运动 → 预测函数或神经元分布的演化 | 核方法、谱松弛与平均场 |
| [在线学习的重叠变量](/notes/statistical_physics_of_learning/online_order_parameters/) | 权重更新 → 局部场分布 → 协方差闭合 | 高斯概率与自平均动力学 |
| [RHM](/notes/statistical_physics_of_learning/designing_hierarchy/) | 层次表示 → 随机规则 → 局部信号的可检测性 | 组合概率与统计涨落 |
| [扩散模型](/notes/statistical_physics_of_learning/diffusion_score/) | 经验 score → 对数配分函数梯度 → 样本权重凝聚 | 高维概率与玻璃统计物理 |
| [SGD 随机动力学](/notes/statistical_physics_of_learning/stochastic_dynamics/) | 局部轨迹或随机高维系统 → 对应的有效随机过程 | OU 与非平衡 DMFT |
| [层间现象论](/notes/statistical_physics_of_learning/phenomenological_model/) | 表示的集体变化 → 可测宏观变量的候选耦合 | 阈值动力学与弹簧—滑块模型 |

各篇开头用一张简短的“问题转译卡”标明入口，正文保留不同的展开方式。对原论文逻辑的重构不等于作者真实的发现历史；自己的玩具模型与实验则另外注明。尤其是 [RHM 的研究练习](/notes/statistical_physics_of_learning/hierarchy_exercise/)，完整记录了两层规则计算、固定世界的检测与代码。

[文献地图](/notes/statistical_physics_of_learning/literature_map/)保留按主题组织的阅读资料。相关站内材料还有[神经网络理论与可解模型](/notes/neural_network_theory/)、[Scaling Laws](/notes/scaling_laws/)、[层次结构数据](/notes/hierarchical_data/)和[扩散模型的动力学](/notes/diffusion_dynamics/)。
