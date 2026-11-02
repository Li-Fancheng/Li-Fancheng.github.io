---
title: "随机生长和分裂下，细胞大小为何稳定？"
subtitle: "Adder：细胞周期问题变成代际随机映射"
date: 2026-11-02
summary: "Adder：细胞周期问题变成代际随机映射"
layout: cell-cycle
article_class: neural-theory-note
math: true
toc: true
weight: 3
---

## 问题转译卡

| 要问的事 | 本篇的选择 |
| --- | --- |
| 原始问题 | 每代生长、分裂都不完全一样，体积分布为什么不发散？ |
| 关键转译 | 沿细胞谱系追踪出生体积，用随机映射描述下一代。 |
| 理论入口 | 自回归过程、稳态分布和代际相关。 |
| 保留与舍弃 | 保留加入体积的随机性；忽略不对称分裂、测量误差及细胞周期机制。 |
| 逻辑地位 | Adder 是操作性规律，不直接说明分子控制机制。 |

[Taheri-Araghi 等 (2015)](https://pubmed.ncbi.nlm.nih.gov/25544609/)在所研究的细菌和条件下发现，出生到分裂之间增加的体积近似与出生体积无关。这使复杂的生长调控能先变成一个代际问题。

## 一个可算到底的随机递推

设出生体积为 $V_n$，该代增加 $\Delta_n$，对称分裂后

$$V_{n+1}=\frac{V_n+\Delta_n}{2}.$$

若 $\Delta_n$ 独立同分布，均值 $\mu$、方差 $\sigma_\Delta^2$，且与当前体积独立，稳态满足

$$\langle V\rangle=\mu,\qquad \operatorname{Var}(V)=\frac{\sigma_\Delta^2}{3},\qquad \operatorname{Corr}(V_n,V_{n+k})=2^{-k}.$$

推导只需取均值、方差并用独立性消去交叉项。偏大的细胞下一代只保留一半偏差，因此存在恢复；但上述数值依赖理想的独立噪声和对称分裂。

## Sizer、Adder 和 Timer 怎样比较？

用统一的局部模型 $V_{n+1}=aV_n+b+\eta_n$，分别取 $a=0,1/2,1$ 表示理想 sizer、adder 和固定倍数生长后对称分裂的 timer。前两者在 $|a|<1$ 时可有稳态；$a=1$ 加非零独立加性噪声时方差随代数增长。不过真实 timer 的噪声常是乘性的，不能直接套用这条加性结论。

模拟多个谱系，匹配平均出生体积和噪声，再比较 $V_{\rm div}$ 对 $V_{\rm birth}$ 的斜率、出生体积自相关与方差随代数的变化。理想 adder 的分裂体积斜率为 1、增加体积对出生体积斜率为 0。注意同一个出生体积既作为横轴又被用于计算增加量，测量误差会制造负相关；应模拟测量误差并使用独立校准。谱系采样与随机抽取总体细胞也有不同的生长偏倚。

[← 相同基因的细胞为什么表达不同？](/notes/solvable_biology/transcription_bursts/) · [专题目录](/notes/solvable_biology/) · [组织密度不变，为什么仍能变得容易流动？ →](/notes/solvable_biology/tissue_rigidity/)
