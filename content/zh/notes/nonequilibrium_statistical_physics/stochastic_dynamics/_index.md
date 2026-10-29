---
title: "随机动力学：从单条轨迹到概率与历史"
subtitle: "三个物理视角，一份随机过程与工具指南"
date: 2026-10-28
lastmod: 2026-10-29
summary: "从 Brownian 与 OU 理解粒子运动、概率流和历史统计，再按物理机制认识常见随机过程。"
layout: series
show_article_cards: false
weight: 1
math: true
tags: [非平衡统计物理, 随机动力学, 随机过程]
---

同一种随机运动，可以看一次记录中的轨迹，也可以看反复准备后形成的概率云，还可以比较连接相同起终点的所有历史。三种描述各自突出不同信息，却使用同一套局部运动规律。

**I 轨迹 · II 概率流 · III 随机历史 · IV 常见随机过程与工具。** 前三篇沿着连续的物理问题展开；第四篇是一份可以单独查阅的过程指南。自由 Brownian 和线性 OU 是主线，常漂移与驱动圆环提供最简单的输运和非平衡稳态例子。

{{% shared-note "nonequilibrium_statistical_physics/stochastic_dynamics" %}}

## 从一次运动到许多可能的历史

1. **[Langevin Dynamics：从牛顿方程到随机轨迹](langevin/)**\
   为什么确定性力学仍需要随机力？从环境作用到阻尼与噪声，实际计算速度记忆、弹道到扩散的跨越，以及谐势中持续的涨落和回复。
2. **[Fokker–Planck Dynamics：从随机轨迹到概率流](fokker_planck/)**\
   许多粒子的概率云如何移动、展宽并趋于稳态？从概率守恒和无偏游走建立漂移流、扩散流，再看 OU 的宽度、圆环持续流与边界。
3. **[Path Integrals：从终点分布到随机历史](path_integrals/)**\
   同一个终点由哪些历史抵达？用时间切片构造权重，计算 Brownian bridge、偏离典型运动的终点、轨迹面积和历史的方向性。
4. **[Stochastic Processes：随机过程与随机微积分导览](stochastic_processes/)**\
   当运动出现跳跃、持续直行、重尾位移、长等待或相关增量，应怎样选模型？对照 Wiener、OU、Poisson、持续游走、Lévy、CTRW 与 fBm，并查阅 Itô/Stratonovich 规则。

前三篇不要求预先掌握生成元、半群、Feynman–Kac 或响应场。第四篇既可接着读，也适合读文献时按需返回；其中 Itô 等距、CTRW 的 Laplace 推导和时间平均 MSD 的细算放在可展开的选读框，首次阅读可以跳过。

## 主线模型与统一记号

前三篇用 $x(t)$ 表示粒子位置，$v(t)=\dot x(t)$ 表示保留惯性时的速度。第一篇先从牛顿受力出发：外力 $F(x)=-V'(x)$，环境提供阻力与随机力。在过阻尼尺度上，同一个外力造成的平均漂移速度为 $\mu F(x)$，其中 $\mu=1/\zeta$。

为了也能讨论一般有效模型，后面把平均漂移速度统记为 $b(x)$；平衡热浴模型中 $b(x)=\mu F(x)$。约化到这一尺度后，主线采用一维、时间齐次、加性 Gaussian 扰动：

$$
\boxed{dx(t)=b(x(t))dt+\sqrt{2D}\,dW_t}.
$$

扩散系数 $D$ 为正，单位是位置平方／时间。$dx(t)$ 表示随机位置增量，$dW_t$ 保留 Wiener 过程的标准记号，满足 $\langle(W_{t+h}-W_t)^2\rangle=h$。形式白噪声 $\xi=\dot W$ 通过冲量或时间积分理解。

概率密度记为 $p(x,t)$，给定初态的转移密度记为 $G(x,t\mid x_0,t_0)$，概率流记为 $J(x,t)$。这些函数中的 $x$ 是空间坐标，$x(t)$ 才是运动记录中的位置。位置方差用 $\sigma^2$ 表示，以免与速度 $v$ 混淆；下标 $\mathrm{st}$ 统一表示稳态。普通系综平均用 $\langle\cdot\rangle$，条件平均也可用 $\mathbb E[\cdot\mid\cdot]$。第四篇及过程论选读采用标准随机过程记号 $X_t$，第四篇开头会解释它与 $x(t)$ 的对应。

| 模型 | 运动规律 | 每次出现的新任务 |
| --- | --- | --- |
| 自由 Brownian | $b=0$ | 扩散尺度、Gaussian 核、时间切片、bridge 与轨迹面积 |
| 线性回复 OU | $b=-\gamma x$ | 旧扰动的遗忘、概率云宽度、关联与过程分类 |
| 常漂移与驱动圆环 | $b=v$；必要时对周长 $L$ 取模 | 移动的云、罕见终点、持续稳态流和绕圈方向 |

$V(x)$ 统一表示真实势能，$F=-V'$；$\zeta$ 为摩擦系数，$\mu=1/\zeta$ 为迁移率。**平衡热浴**中才有 $b=-\mu V'$、$D=\mu k_{\mathrm B}T$。对谐势 $V=\kappa x^2/2$，OU 回复率 $\gamma=\mu\kappa$。一般有效随机模型的参数不必具有温度解释。

## 进阶选读：按问题进入

以下材料保存原有深入计算，独立于四篇主线；不需要从头再读一遍。

| 想进一步研究什么 | 选读材料 | 保留的内容 |
| --- | --- | --- |
| 为什么换阱极罕见 | [双稳亚稳态与罕见历史](rare_events/) | 跨阱输运、两态约化、上坡路径与势垒指数 |
| 到达多久、先到哪里、为何弛豫慢 | [生成元、首次到达与弛豫谱](markov_generators/) | 首达积分、committor、Kramers 前因子、谱、主方程生成元与 Feynman–Kac |
| 怎样系统计算多时刻关联与响应 | [Gaussian 泛函与 MSRJD 响应场](path_integral_methods/) | OU Gaussian 计算、响应场离散推导、因果性与 OM Jacobian |

阅读起点是牛顿力学、微积分和基础概率。需要从离散反应事件或生命系统观测出发时，可对照 QBIO 的[生命系统中的随机动力学](/notes/qbio_course/stochastic_dynamics/)。
