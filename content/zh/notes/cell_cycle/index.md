---
title: "细胞周期：从单细胞噪声到群体增长"
subtitle: "Cell-Cycle Regulation, Temporal Correlations, and Stochastic Dynamics"
date: 2026-10-14
summary: "从独立分裂计时模型出发，推导 Euler–Lotka 方程、时间关联、adder 尺寸稳态和涨落功率谱，理解单细胞噪声如何影响群体增长。"
translationKey: cell-cycle-correlations
layout: cell-cycle
math: true
toc: true
tags:
  - Cell Cycle
  - Stochastic Dynamics
  - Statistical Physics
---

这份英文学习笔记围绕一个问题展开：分裂时间的随机性为什么能加速群体增长，而维持尺寸稳态的细胞又如何通过代际关联抵消这一效应？内容包括模型假设、逐步推导、实验图示、技术附录和练习题。

**阅读路线：** 独立计时模型 → Euler–Lotka 方程 → 时间关联 → adder 尺寸调控 → 功率谱与振荡。建议先读第 0–2 节了解问题，再读第 6 节的核心推导。

{{% cell-cycle-note %}}
