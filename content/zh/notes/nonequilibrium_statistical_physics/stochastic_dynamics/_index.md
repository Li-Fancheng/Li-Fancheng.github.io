---
title: "随机动力学：从单条轨迹到概率与历史"
subtitle: "三个物理视角，一份随机过程与工具指南"
date: 2026-10-28
lastmod: 2026-10-30
summary: "从 Brownian 与 OU 理解粒子运动、概率流和历史统计，再按物理机制认识常见随机过程。"
layout: series
show_article_cards: false
weight: 1
math: true
tags: [非平衡统计物理, 随机动力学, 随机过程]
---

同一种随机运动，可以看一次记录中的轨迹，也可以看反复准备后形成的概率云，还可以比较连接相同起终点的所有历史。三种描述各自突出不同信息，却使用同一套局部运动规律。

**I 轨迹 · II 概率流 · III 随机历史 · IV 三种描述的联系 · V 常见随机过程与工具。** 前三篇沿着连续的物理问题展开，第四篇整理它们的相互推导；第五篇是一份可以单独查阅的过程指南。自由 Brownian 和线性 OU 是主线，常漂移与驱动圆环提供最简单的输运和非平衡稳态例子。

{{% shared-note "nonequilibrium_statistical_physics/stochastic_dynamics" %}}

## 从一次运动到许多可能的历史

1. **[Langevin Dynamics：从牛顿方程到随机轨迹](langevin/)**\
   为什么确定性力学仍需要随机力？从环境作用到阻尼与噪声，实际计算速度记忆、弹道到扩散的跨越，以及谐势中持续的涨落和回复。
2. **[Fokker–Planck Dynamics：从随机轨迹到概率流](fokker_planck/)**\
   许多粒子的概率云如何移动、展宽并趋于稳态？从概率守恒和无偏游走建立漂移流、扩散流，再看 OU 的宽度、圆环持续流与边界。
3. **[Path Integrals：从终点分布到随机历史](path_integrals/)**\
   同一个终点由哪些历史抵达？用时间切片构造权重，计算 Brownian bridge、偏离典型运动的终点和轨迹面积。
4. **[三种描述的联系：从随机方程到概率与路径](relations/)**\
   从短时转移连接 Langevin、Fokker–Planck 和路径积分，双向推导并用 OU 过程对照三种计算。
5. **[Stochastic Processes：随机过程与随机微积分导览](stochastic_processes/)**\
   当运动出现跳跃、持续直行、重尾位移、长等待或相关增量，应怎样选模型？对照 Wiener、OU、Poisson、持续游走、Lévy、CTRW 与 fBm，并查阅 Itô/Stratonovich 规则。

前三篇适合顺次阅读；第四篇用于复习与整理，第五篇既可接着读，也可在需要具体随机过程或微积分工具时查阅。

## 主线模型与统一记号

第一篇从牛顿受力出发：势能产生的力为 $F=-dV/dx$，环境提供阻力与随机力。在过阻尼尺度上，同一个力造成的平均漂移速度为 $\mu F(x)$，其中 $\mu=1/\zeta$ 为迁移率。

为了也能讨论一般有效模型，后面把平均漂移速度统记为 $b(x)$；平衡热浴模型中 $b(x)=\mu F(x)$。约化到这一尺度后，主线采用一维、时间齐次、加性 Gaussian 扰动：

$$
\boxed{dx=b(x)dt+\sqrt{2D}\,dW_t}.
$$

扩散系数 $D$ 为正，单位是位置平方／时间。标准 Wiener 过程满足 $\langle(W_{t+h}-W_t)^2\rangle=h$，形式白噪声 $\xi=\dot W$ 通过冲量或时间积分理解。

概率密度记为 $p(x,t)$，给定初态的转移密度记为 $G(x,t\mid x_0,t_0)$，概率流记为 $J(x,t)$。位置方差用 $\sigma^2$ 表示，下标 $\mathrm{st}$ 表示稳态。系综平均用 $\langle\cdot\rangle$，条件平均也可用 $\mathbb E[\cdot\mid\cdot]$。第五篇采用标准随机过程记号 $X_t$。

| 模型 | 运动规律 | 每次出现的新任务 |
| --- | --- | --- |
| 自由 Brownian | $b=0$ | 扩散尺度、Gaussian 核、时间切片、bridge 与轨迹面积 |
| 线性回复 OU | $b=-\gamma x$ | 旧扰动的遗忘、概率云宽度、关联与过程分类 |
| 常漂移与驱动圆环 | $b=v$；必要时对周长 $L$ 取模 | 移动的云、罕见终点、持续稳态流和绕圈方向 |

**平衡热浴**中才有 $b=-\mu\,dV/dx$、$D=\mu k_{\mathrm B}T$。对谐势 $V=\kappa x^2/2$，OU 回复率 $\gamma=\mu\kappa$。一般有效随机模型的参数不必具有温度解释。

阅读起点是牛顿力学、微积分和基础概率。
