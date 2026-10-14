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

**阅读路线：** 第 0–2 节提出独立噪声的增长优势与尺寸稳态之间的矛盾；第 3 节把 Euler–Lotka 方程推广到相关的分裂历史，第 4 节说明累积时间涨落如何影响增长率；第 5 节用 Markov 模型展示相同的单代分裂时间分布如何因时间关联不同而产生不同增长率；第 6 节推导 adder 如何产生关联并抵消累积时间涨落；第 7–8 节讨论功率谱与振荡，第 9 节总结。若主要关注尺寸调控的核心推导，建议按第 0–4 节 → 第 6 节的顺序阅读；第 5 节可作为具体例子选读。

{{% cell-cycle-note %}}
