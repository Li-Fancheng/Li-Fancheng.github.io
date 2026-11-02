---
title: "Li₂ 与 grokking 的特征涌现动力学"
subtitle: "从输出层消元到单神经元能量景观"
date: 2026-11-02
lastmod: 2026-11-02
summary: "分析输出层怎样构造标签相关的反传信号，以及独立特征近似怎样将 grokking 连接到能量极值、群表示和有限样本稳定性。"
layout: cell-cycle
article_class: neural-theory-note
math: true
toc: true
---

## 1. 现象与需要解释的问题

Grokking 的张力在于，训练准确率已经很高，测试准确率却可能在很久以后才上升。训练集拟合完成，并不表示表示学习已经完成；但若只说网络由记忆转向理解，仍没有解释梯度为何继续改变隐藏层，以及它会选择什么特征。

田渊栋的 [Li₂ 框架](https://arxiv.org/html/2509.21519v5)把两层非线性网络的学习分为 lazy、independent 和 interactive 三个阶段。论文的具体分析见 *Provable Scaling Laws of Feature Emergence from Learning Dynamics of Grokking*。这里关心三个相连的问题：拟合输出层怎样产生有用的反传信号，信号怎样选择特征，以及一部分特征学会以后，网络怎样继续寻找剩余结构。

需要观测的不只有训练与测试误差，还包括反传梯度的方向、神经元之间的相关，以及隐藏权重在任务相关子空间中的投影。相似的误差平台可能对应完全不同的内部过程。

## 2. 关键自由度：隐藏激活与反传信号

取 $n$ 个样本，数据矩阵 $X\in\mathbb R^{n\times d}$、one-hot 标签 $Y\in\mathbb R^{n\times M}$，网络为

$$
F=\sigma(XW),\qquad \hat Y=FV,
\qquad W\in\mathbb R^{d\times K},\quad V\in\mathbb R^{K\times M}.
$$

令 $P=I-\mathbf1\mathbf1^{\mathsf T}/n$，$\tilde F=PF$、$\tilde Y=PY$，使用沿样本维去均值的平方损失。隐藏层收到的负梯度精确为

$$
G_F=-\frac{\partial J}{\partial F}
=P(Y-FV)V^{\mathsf T},\qquad
\dot w_j=X^{\mathsf T}D_j[G_F]_j,
\quad D_j=\operatorname{diag}\sigma'(Xw_j).
$$

因此残差大小本身不够：残差还要经过 $V^{\mathsf T}$，才能决定每个隐藏节点怎样移动。$\tilde Y\tilde Y^{\mathsf T}$ 比较样本的标签关系，$\tilde F^{\mathsf T}\tilde F$ 比较神经元激活，$\tilde F\tilde F^{\mathsf T}$ 则描述当前特征已覆盖的样本方向。这些矩阵直接进入反馈，因而比所有权重坐标更适合组织问题。

## 3. 关键近似与模型选择

第一步是暂时固定隐藏特征，求输出层的线性学习。这保留了输出层如何解释当前表示的过程，并允许随后把它的反馈代回隐藏层。把固定特征解用于移动的 $W$，还需要输出层弛豫比表示变化快的准静态假设；它在固定 $F$ 下的精确性，不自动成为全程训练的精确性。

第二步才是独立特征近似：将反馈中的样本空间滤波和神经元空间耦合近似为标量作用，使 $[G_F]_j\simeq\kappa\tilde Y\tilde Y^{\mathsf T}\tilde f_j$。保留下来的标签相关算子仍能推动 $w_j$ 发生有限变化。这比冻结整个隐藏层多保留了一项机制，同时比一般耦合训练少追踪了神经元间的反馈。

原文用随机初始化、输入几何和大宽度讨论这种结构。应用时应检查矩阵对实际信号的作用：中心化会产生零模态，矩形 $F$ 的两个 Gram 矩阵还受秩约束，不能在任意宽高比下把二者都当作全空间的非零单位阵。

最后，为研究特征方向的极值，需要固定范数。下文采用 $\|w_j\|=1$ 的球面模型；普通权重衰减会影响半径和时间尺度，却不与这一约束严格等价。群算术和平方激活则进一步选择了一个可以分析极值结构的数据模型。

## 4. 从神经网络问题到可计算的理论问题

### 输出层消元留下了什么

先只对输出层加入 $\lambda\|V\|_F^2/2$、$\lambda>0$，损失按样本求和，不除以 $n$。固定 $F$ 时，

$$
\dot V=\tilde F^{\mathsf T}(\tilde Y-\tilde FV)-\lambda V,
\qquad V_*=B\tilde F^{\mathsf T}\tilde Y,
\quad B=(\tilde F^{\mathsf T}\tilde F+\lambda I)^{-1}.
$$

定义软残差算子

$$
A_\lambda=I-\tilde FB\tilde F^{\mathsf T}
=\lambda(\tilde F\tilde F^{\mathsf T}+\lambda I)^{-1}.
$$

将 $V_*$ 代回原始反传式，得到固定特征输出层平衡处的恒等式

$$
G_F=A_\lambda\tilde Y\tilde Y^{\mathsf T}\tilde FB.
$$

输出层的消去产生两个不同的反馈：左边 $A_\lambda$ 强调当前特征难以拟合的样本方向，右边 $B$ 调整冗余神经元的组合。独立近似 $G_F\simeq\kappa\tilde Y\tilde Y^{\mathsf T}\tilde F$ 恰好删除了这两类结构的方向差异。原文用 $\eta$ 表示这里的衰减系数；本页以 $\lambda$ 与学习率区分。

有用信号也可以在输出层到达平衡以前出现。若 $V(0)=0$，固定 $F$ 的短时展开直接给出

$$
V(t)=t\tilde F^{\mathsf T}\tilde Y+O(t^2),\qquad
G_F(t)=t\tilde Y\tilde Y^{\mathsf T}\tilde F+O(t^2).
$$

这一项不要求 $\lambda>0$。衰减可以维持后期反馈，初始瞬态也可以提供特征学习信号；这两条路径在[原文 v5 的第 3 节](https://arxiv.org/html/2509.21519v5#S3)中已明确区分。若无衰减且已精确插值，平衡残差才会为零，不能把这个末态结论反推到全部早期动力学。

### 独立反传成为单神经元的能量上升

在 $G_F\simeq\kappa\tilde Y\tilde Y^{\mathsf T}\tilde F$、$\kappa>0$ 的窗口内，定义

$$
\mathcal E(w)=\tfrac12\|\tilde Y^{\mathsf T}\sigma(Xw)\|^2,
\qquad
\nabla_w\mathcal E=X^{\mathsf T}D(w)\tilde Y\tilde Y^{\mathsf T}\sigma(Xw).
$$

链式法则说明隐藏节点沿 $\mathcal E$ 上升。球面模型中，

$$
\dot w=\kappa(I-ww^{\mathsf T})\nabla\mathcal E,
\qquad
\frac{d\mathcal E}{dt}=\kappa\|(I-ww^{\mathsf T})\nabla\mathcal E\|^2\ge0.
$$

这个能量度量一个非线性特征与标签的相关强度。它由反馈结构导出；在该近似内，稳定的特征方向对应球面上的局部极大值。原来要解整个网络的梯度轨迹，现在可以先问一个节点的能量有哪些极值、吸引域怎样分布、抽样会怎样扰动这些极值。

### 模加法为什么选择 Fourier 特征

考虑完整的 $\mathbb Z_M$ 模加法表：输入为 $(e_a,e_b)$，标签为 $h=a+b\pmod M$，$w=(u;v)$，取 $\sigma(z)=z^2$。以下用全部 $M^2$ 个样本的未归一化求和，直接计算系数。限制到 $\sum_a u_a=\sum_b v_b=0$ 的非常数方向，则

$$
[\tilde Y^{\mathsf T}f]_h
=\sum_a(u_a+v_{h-a})^2
-\frac1M\sum_{a,b}(u_a+v_b)^2
=2(u*v)_h.
$$

线性激活下同一个中心化标签和为零：各类别见到的单个输入元素完全均衡。平方激活保留的交叉项 $u_av_b$ 才能携带和的结构。因此这里去掉非线性会删掉目标机制，与[深度线性案例](/notes/statistical_physics_of_learning/deep_linear_model/)的取舍不同。

定义单位化 Fourier 变换 $\hat u_k=M^{-1/2}\sum_a u_ae^{-2\pi\mathrm i ka/M}$。卷积被对角化，Parseval 恒等式给出

$$
\mathcal E=2\sum_h|(u*v)_h|^2
=2M\sum_k|\hat u_k|^2|\hat v_k|^2,
\qquad \|u\|^2+\|v\|^2=1.
$$

能量奖励两半权重集中到相同频率。对实权重、奇数 $M$，频率以 $k,-k$ 配对；只比较两组不同频率对，并令两半的范数各为 $1/\sqrt2$、频率质量均为 $r,1-r$，有

$$
\mathcal E(r)=\frac M4\bigl[r^2+(1-r)^2\bigr].
$$

均匀混合是谷底，集中于任一频率对提高能量。单频内部的相位仍可连续变化，所以是带平坦方向的极大值集合。这里 Fourier 方法能用，来自标签运算把类条件相关变成循环卷积；原文对一般有限群进一步使用不可约表示分块。

### 把被删掉的相互作用放回来

当多个节点学到相似特征，$B$ 的非对角元重新重要。只取两个等范数的中心化激活，令其 Gram 矩阵为 $\left(\begin{smallmatrix}a&c\\c&a\end{smallmatrix}\right)$，则精确有

$$
[\tilde FB]_1=
\frac{(a+\lambda)\tilde f_1-c\tilde f_2}{(a+\lambda)^2-c^2}.
$$

正相关使第二项相减，削弱重复方向。与此同时，$A_\lambda$ 在样本核特征值 $s$ 的方向上作用为 $\lambda/(s+\lambda)$：已经充分覆盖的方向被压低，尚未解释的标签分量继续反馈。这为后期的特征多样性和缺失结构提供了具体来源。它们在独立能量中被省略，在 interactive 阶段重新控制学习。

## 5. 物理机制、理论预测与检验

物理图像是一个会改变反馈的学习过程：输出层先利用现有表示，标签相关信号随后推动节点选择任务结构，已学会的结构又改变后续节点收到的力。单节点偏向集中频率，整个网络却需要覆盖足够多的特征；后期的竞争与残差反馈协调这两种要求。

样本量也会改变特征景观。[原文定理 4](https://arxiv.org/html/2509.21519v5#Thmtheorem4)通过经验能量的扰动分析，给出保持 $d_k$ 维不可约表示相关极大值的充分样本尺度

$$
n\gtrsim d_k^2M\log(M/\delta).
$$

对循环群，复不可约表示维数为 1，因而得到 $M\log M$ 的量级。它控制特征极值的稳定性；完整网络还需要到达这些极值并覆盖所需特征，才能推出泛化。这样，样本阈值有了比测试准确率跳变更具体的中间对象。

按上述求和规范，对 $M=11$ 的完整表直接计算：单频对得到 $\mathcal E=2.75$，两频对等质量混合得到 $1.375$，与公式一致。有限差分也核对了能量梯度与独立反传式；这些是本页方程的数值核对，并未重跑原论文的端到端 grokking 实验。作者的[实验代码](https://github.com/yuandong-tian/understanding/tree/main/ssl/real-dataset/cogo)提供了进一步比较的起点。

一个有判别力的检查是随训练测量

$$
\varepsilon_{\rm ind}(t)=
\min_{\kappa\ge0}\frac{\|G_F-\kappa\tilde Y\tilde Y^{\mathsf T}\tilde F\|_F}{\|G_F\|_F},
\qquad \|G_F\|_F>0,
$$

再与频谱集中、节点冗余和测试误差一起比较。它直接检验独立近似是否存在有效窗口。只观察测试误差变好，不能区分节点先独立选频再互相协调，还是从一开始就由强耦合反馈选择特征。

## 6. 尚未解决的矛盾与进一步的问题

能量景观解释了特征的候选方向，却还没有单独给出等待多久才泛化。等待时间还涉及信号强度随时间的积累、初始投影、输出层弛豫和后期特征覆盖。值得先比较相同数据与隐藏初始化下的零输出初始化、小随机输出初始化及不同输出学习率，同时测 $G_F$ 的结构。这样能分辨平台主要受输出层生成信号的速度控制，还是受已有信号驱动特征移动的速度控制。

另一个问题来自单频选择与全网覆盖之间的张力。可在输出层准静态的对照系统中，分别使用完整反馈、只保留 $A_\lambda$、只保留 $B$ 和独立近似，比较学到的频率种类与测试误差。这能把单节点的景观选择、节点去冗余和缺失特征反馈分开检验，也能判断独立近似究竟解释了 grokking 的哪一段。
