---
title: "06 · NTK：非线性网络怎样产生线性训练动力学"
subtitle: "切线特征、冻结核与训练集外预测"
date: 2026-10-19
weight: 6
summary: "由参数梯度流推导时变 NTK，再求解冻结核下的残差、测试预测和 RKHS 插值，并给出失效诊断。"
layout: cell-cycle
article_class: neural-theory-note
math: true
toc: true
reading_time: 17
tags:
  - 神经网络理论
  - 统计物理
  - 可解模型
---

{{% shared-note "neural_network_theory/ntk" %}}
