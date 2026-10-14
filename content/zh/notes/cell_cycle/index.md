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

**阅读路线：** 第 0–2 节提出独立噪声的增长优势与尺寸稳态之间的矛盾；第 3–4 节建立时间关联与群体增长的联系；第 5 节用简短对照展示相同的分裂时间分布如何产生不同增长率；第 6 节推导 adder 如何通过尺寸反馈产生关联并抵消累积时间涨落；第 7–8 节讨论功率谱与振荡，第 9 节总结。核心推导可按第 0–4 节 → 第 6 节阅读。第 5 节提供直观例子，完整的 Markov 模型计算放在附录 E。

{{% cell-cycle-note %}}
