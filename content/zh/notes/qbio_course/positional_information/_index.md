---
title: "基因表达如何编码空间位置"
subtitle: "从 Bicoid 梯度和 French flag model 到位置解码与网络优化"
date: 2026-10-17
summary: "从 Bicoid 梯度出发，介绍基因表达中的位置信息与位置解码，讨论噪声对定位精度的限制，以及多基因联合读取和调控网络优化的作用。"
layout: qbio-series
article_class: positional-information-note
math: true
toc: false
tags:
  - 位置信息
  - 基因调控
  - 统计物理
qbio_topic: positional_information
weight: 2
---

细胞核接触的是分子，却需要在合适的位置作出表达选择。第一篇从梯度的形成讲到有噪声的表达边界；第二篇讨论多个基因怎样共同定位，以及局部精度为何不等于全局可辨认；第三篇再问，哪些分子网络能够产生可靠的位置图案。

各篇沿用原有章节、公式和图号。补充推导可按需要展开，跨篇引用可直接跳转。
