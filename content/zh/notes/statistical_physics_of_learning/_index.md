---
title: "从神经网络问题到可解理论：问题重构与物理近似"
date: 2026-10-30
lastmod: 2026-11-02
summary: "从学习现象出发，分析自由度选择、物理近似与可计算理论怎样共同构成机制解释。"
layout: series
show_article_cards: false
article_class: neural-theory-note
math: true
tags:
  - 统计物理
  - 神经网络理论
  - 文献阅读
---

神经网络的理论问题常有多个相互纠缠的来源：训练改变表示，有限数据限制可估计的方向，随机更新带来历史相关，而数据本身又含多层结构。研究首先需要确定要解释哪一项现象，以及哪些自由度控制它。

这个专题分析具体工作怎样选择变量、构造模型或取极限，让原来难以计算的问题获得新的数学结构。随机协方差的逆矩阵、独立模式的动力学、局部场的自洽过程、样本能量的配分函数，分别计算不同的观测。最终还要回到测试误差、训练时间、表示几何或生成记忆，判断理论解释了什么，以及什么实验能区分剩余机制。

## 从现象到理论对象

| 研究笔记 | 要解释的现象与保留的自由度 | 转化后的计算 |
| --- | --- | --- |
| [Scaling laws](/notes/statistical_physics_of_learning/spectral_scaling/) | 有限数据与有限表示怎样改善泛化；协方差、随机投影与目标谱 | 固定特征回归、Wishart 逆矩阵、未解析方向的风险 |
| [深度线性网络](/notes/statistical_physics_of_learning/deep_linear_model/) | 平台和快速学习；数据奇异模态与分层增益 | 特殊对齐轨道的非线性 ODE、初始化时间尺度 |
| [NTK 与平均场](/notes/statistical_physics_of_learning/changing_variables/) | 预测怎样演化，表示是否移动；切核或神经元分布 | 核谱松弛、测试预测、分布输运 |
| [在线学习的重叠变量](/notes/statistical_physics_of_learning/online_order_parameters/) | 随机输入流中的对齐和泛化；局部场协方差 | 高斯积分、宏观闭合、观测量所需的最小变量 |
| [RHM 与层级学习](/notes/statistical_physics_of_learning/designing_hierarchy/) | 组合结构怎样改变数据需求；固定规则与局部相关 | 概率传播、信号与样本涨落、检测与学习的比较 |
| [扩散模型](/notes/statistical_physics_of_learning/diffusion_score/) | 类别分化与样本记忆；经验 score、能量和权重 | 对数配分函数、典型能量鞍点与凝聚 |
| [SGD 随机动力学](/notes/statistical_physics_of_learning/stochastic_dynamics/) | 局部涨落或全程学习；曲率、噪声及两时响应 | OU 协方差、带记忆的自洽有效过程 |
| [层间现象论](/notes/statistical_physics_of_learning/phenomenological_model/) | 表示怎样联动；可测分离与候选耦合 | 阈值动力学、启动时间及跨条件预测 |
| [两层随机规则计算](/notes/statistical_physics_of_learning/hierarchy_exercise/) | 规则平均能否代表固定世界的可检测性 | 超几何矩、信号分布、有限样本实验与复现代码 |

精确变量重写、受控极限和机制模型支持的结论不同。函数动力学的链式法则本身不冻结核，SVD 本身不消去模式竞争，统计信号可检测也不保证网络学会它。各篇沿现象、自由度、近似、计算、检验和解释缺口展开，并标明原文结果与自构例子的关系。

## 延伸阅读

[文献地图](/notes/statistical_physics_of_learning/literature_map/)按科学问题比较 27 篇工作，帮助查找某个解释缺口对应的变量和理论。相关站内材料包括[神经网络理论与可解模型](/notes/neural_network_theory/)、[Scaling Laws](/notes/scaling_laws/)、[层级学习](/notes/hierarchical_data/)和[扩散模型动力学](/notes/diffusion_dynamics/)。
