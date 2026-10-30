---
title: "神经网络的统计物理：好工作的近似方法"
date: 2026-10-30
summary: "以近似和可解极限为线索，阅读 scaling laws、扩散模型、层次数据与特征学习中的统计物理工作。"
layout: series
show_article_cards: false
article_class: neural-theory-note
math: true
tags:
  - 统计物理
  - 神经网络理论
  - 文献阅读
---

做理论研究时，我常觉得最困难的不是计算，而是不知道应该计算什么。面对复杂系统，有很多因素可以舍弃，也有很多极限可以取：哪些结构必须留下？一个简化究竟能回答什么问题？

这里想从几篇有启发性的工作出发，练习理解研究者怎样把复杂问题变成有物理内容的简单问题。不同做法不一定是同一种意义上的近似：有的是受控极限，有的是隔离机制的最小模型，还有的是人为设计数据。它们能支持的结论也不同。

## 从问题出发的案例

- [深度线性网络：去掉非线性，为什么还能看到复杂的学习动力学？](/notes/statistical_physics_of_learning/deep_linear_model/) 从平台期出发，保留层间乘法，分离深度参数化造成的动力学。
- [NTK 与平均场：先换变量，再决定怎样闭合](/notes/statistical_physics_of_learning/changing_variables/) 从函数动力学的恒等式出发，比较固定切核与移动的神经元分布。
- [RHM：如果想研究层次结构，应该怎样造数据？](/notes/statistical_physics_of_learning/designing_hierarchy/) 从高斯混合开始，构造包含组合与等价替换的数据世界。

- [扩散模型：如果先假设 score 已经学好，还剩下什么问题？](/notes/statistical_physics_of_learning/diffusion_score/) 用可写出的加噪分布分离生成动力学与网络训练。
- [SGD 的噪声什么时候能当成热噪声？](/notes/statistical_physics_of_learning/stochastic_dynamics/) 从离散局部模型到 OU 近似，再问何时需要记忆和自洽反馈。
- [Scaling laws：一条幂律背后，究竟是哪种资源在受限？](/notes/statistical_physics_of_learning/spectral_scaling/) 分别计算样本涨落与未解析模态的最小模型。
- [微观方程太复杂时，怎样构造现象论？](/notes/statistical_physics_of_learning/phenomenological_model/) 从可测层间变化出发，试探阈值与耦合如何形成可检验模型。

[文献地图](/notes/statistical_physics_of_learning/literature_map/)保留了原来的论文资料，可在需要时按主题继续查阅。也可参看[神经网络理论与可解模型](/notes/neural_network_theory/)、[Scaling Laws](/notes/scaling_laws/)、[层次结构数据](/notes/hierarchical_data/) 和 [扩散模型的动力学](/notes/diffusion_dynamics/)。
