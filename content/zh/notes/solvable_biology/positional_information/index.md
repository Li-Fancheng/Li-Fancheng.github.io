---
title: "细胞怎样从带噪声的信号推断位置？"
subtitle: "位置编码：把发育问题变成统计解码"
date: 2026-11-02
summary: "位置编码：把发育问题变成统计解码"
layout: cell-cycle
article_class: neural-theory-note
math: true
toc: true
weight: 5
---

## 问题转译卡

| 要问的事 | 本篇的选择 |
| --- | --- |
| 原始问题 | 相同基因组的细胞如何获得位置相关的命运信息？ |
| 关键转译 | 用表达分布 $P(\boldsymbol g\mid x)$ 构造后验 $P(x\mid\boldsymbol g)$。 |
| 理论入口 | 贝叶斯推断、Fisher 信息和互信息。 |
| 保留与舍弃 | 保留表达均值、噪声与相关性；不预设细胞实际解码电路。 |
| 逻辑地位 | 最优统计解码给出信息上的参照，而非具体命运机制。 |

[Tkačik 等 (2015)](https://pubmed.ncbi.nlm.nih.gov/25361898/)将位置信息、误差和读取精度放在同一数学框架下。对局部表达 $\boldsymbol g$，给定先验 $P(x)$，贝叶斯公式是

$$P(x\mid\boldsymbol g)=\frac{P(\boldsymbol g\mid x)P(x)}{\int P(\boldsymbol g\mid x')P(x')\,dx'}.$$

模型先回答“数据允许多准确地推断位置”，再问实际细胞是否能够达到这种精度。信息存在不等于下游网络真的读取了它。

## 一个局部可算的例子

假设 $\boldsymbol g\mid x$ 服从均值 $\boldsymbol\mu(x)$、与位置无关的协方差 $\Sigma$ 的高斯分布。局部 Fisher 信息为

$$J(x)=\boldsymbol\mu'(x)^\mathsf{T}\Sigma^{-1}\boldsymbol\mu'(x).$$

在噪声较小、均值可局部区分且远离边界时，位置误差可近似为 $\sigma_x\simeq J^{-1/2}$。若噪声也随位置变化，Fisher 信息还有协方差导数项；多峰后验或边界附近也不能把这个局部公式当成全局误差。互信息则依赖位置先验，衡量整个表达向量对位置的不确定性减少，并不等同某一点的误差。

两个基因斜率都为 $s$、噪声方差均为 $\sigma^2$、相关系数为 $\rho$ 时，$J=2s^2/[\sigma^2(1+\rho)]$。同向的共同噪声越强，第二个基因的增益越小；若斜率方向改变，相关性的作用也会改变。

## 计算和检验

在有限区间上生成含噪表达，以数值积分求完整后验，比较单基因与双基因的平均解码误差；保持各基因方差和均值曲线不变，扫描 $\rho$。再用不同先验、边界位置和非单调均值检查局部公式何时失效。实验比较需要分离技术噪声、胚胎间差异和位置配准误差；更高的最优解码精度本身不能证明命运决定使用了该信号。

[← 组织密度不变，为什么仍能变得容易流动？](/notes/solvable_biology/tissue_rigidity/) · [专题目录](/notes/solvable_biology/) · [动物什么时候应该离开食物区域？ →](/notes/solvable_biology/foraging/)
