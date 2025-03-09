---
title: 细菌活性物质的统计力学研究
date: 2024-05-03
math: true
---

论文Statistical Mechanics of Interacting Run-and-Tumble Bacteria 的阅读笔记

## 1 文献概述
### 1.1  Run-and-Tumble
Run-and-Tumble 运动是一些细菌（如 Escherichia coli，大肠杆菌）典型的自推进运动模式。这类细菌的运动方式包括：Run（奔跑）：细菌沿几乎直线的轨迹以恒定速度 \( v \) 运动。Tumble（翻滚）：细菌随机改变方向，发生的概率由 翻滚速率 \( \lambda \) 决定。

这一类运动本质上是一种随机游走，在大尺度下，类似于扩散过程，其扩散系数为：
\[
D = \frac{v^2}{\lambda}
\]
这个扩散系数通常比同温度下的布朗运动扩散系数大得多，因此在研究活性物质（Active Matter）动力学和自驱动粒子统计力学时，Run-and-Tumble 运动是一个重要的模型。Run-and-Tumble 运动是非平衡过程，其稳态分布无法用标准的玻尔兹曼因子 \( e^{-H/k_BT} \) 直接描述。细菌的集体行为（如群体运动、自聚集等）不能简单地用热力学模型解释，而需要专门的非平衡统计物理分析。过去的研究主要集中在单个Run-and-Tumble粒子的扩散模型，但对多体相互作用下的集体行为缺乏深入研究。


### 1.2 主要研究问题
本文探讨了多个 Run-and-Tumble 细菌相互作用时的集体动力学，主要关注以下问题：
- 如何从单粒子的动力学推广到多体系统？过去的研究多关注单个粒子的运动，而忽略了粒子间的相互作用及其对整体动力学的影响。本文通过统计力学的方法，建立描述 Run-and-Tumble 细菌集体行为的模型。

- 如何在统计力学框架下描述 Run-and-Tumble 运动？由于 Run-and-Tumble 运动不满足详细平衡，无法直接套用传统的热力学方法。本文推导了描述细菌密度演化的漂移-扩散方程，并利用 Langevin 形式分析其噪声特性。

- 如何理解细菌的集体行为，如自聚集（self-trapping）和相分离？研究在何种条件下，细菌群体会自发形成高密度和低密度区域，而不是均匀分布。研究 Run-and-Tumble 运动中的“自聚集”现象是否类似于平衡系统中的相分离。

- 在 Run-and-Tumble 体系中，是否可以找到等效的“有效温度”来描述系统？研究当粒子相互作用导致速度变化时，系统能否等效地映射到一个有效温度的热力学系统。探讨何种情况下该映射会失效。

## 2 Run-and-Tumble 运动的基本理论

### 2.1 细菌运动模式（Run 与 Tumble）

Run-and-Tumble 运动的基本特性：

细菌的翻滚发生是随机的，受泊松分布控制，翻滚之间的时间间隔服从指数分布：
  \[
  P(t) = \lambda e^{-\lambda t}
  \]
  其中 \( t \) 是连续 Run 的时间，\( \lambda \) 是翻滚率（即单位时间内翻滚的概率）。每次翻滚后，细菌会在所有方向上随机选择新的运动方向。这一模式在大尺度下呈现出扩散行为，但扩散系数远大于普通布朗运动。这种运动模式可以使细菌在复杂环境中有效搜索食物源。在趋化过程中，细菌可以调节翻滚率 \( \lambda \) 以趋向高营养浓度区域。


### 2.2 单个粒子的随机游走动力学

在一维情况下，Run-and-Tumble 运动可以用两种状态来描述：右移状态（Right-moving）： 速度 \( +v \)，概率分布记作 \( R(x,t) \)。左移状态（Left-moving）： 速度 \( -v \)，概率分布记作 \( L(x,t) \)。设翻滚率为 \( \lambda \)，则单个粒子的演化方程为：

\[
\frac{\partial R}{\partial t} = -v \frac{\partial R}{\partial x} - \frac{\lambda}{2} R + \frac{\lambda}{2} L
\]
\[
\frac{\partial L}{\partial t} = v \frac{\partial L}{\partial x} - \frac{\lambda}{2} L + \frac{\lambda}{2} R
\]
第一项表示粒子以速度 \( v \) 或 \( -v \) 移动。第二项表示粒子以速率 \( \lambda \) 翻滚并失去当前状态。第三项表示粒子从另一状态翻滚至当前状态。将两式相加，定义总概率密度：
\[
p(x,t) = R(x,t) + L(x,t)
\]
求和可得：
\[
\frac{\partial p}{\partial t} = -\frac{\partial J}{\partial x}
\]
其中概率流 \( J(x,t) \) 定义为：
\[
J(x,t) = v R - v L
\]
概率守恒方程表明，Run-and-Tumble 运动可以用漂移-扩散形式来描述。

### 2.3 扩散极限与漂移扩散方程

当粒子的翻滚率 \( \lambda \) 较高，相较于 Run 时间尺度，翻滚发生得很快，可以用连续介质近似，即扩散极限（diffusive limit）。  
对 \( J(x,t) \) 求时间导数，并使用上面的概率流方程：
\[
\frac{\partial J}{\partial t} = v \frac{\partial R}{\partial t} - v \frac{\partial L}{\partial t}
\]
代入 \( R, L \) 的演化方程，并假设在长时间尺度上 \( \frac{\partial J}{\partial t} \approx 0 \)，可以得到稳态近似：
\[
J \approx -D \frac{\partial p}{\partial x}
\]
其中：
\[
D = \frac{v^2}{\lambda}
\]
即 Run-and-Tumble 运动的有效扩散系数。这个结果表明：在大尺度上，Run-and-Tumble 运动表现为扩散行为，扩散系数 \( D \) 取决于粒子的速度 \( v \) 和翻滚率 \( \lambda \)。当 \( v \) 增加或 \( \lambda \) 减小时，扩散性增强，细菌的运动范围扩大。

Run-and-Tumble 运动在大尺度上的动力学可以写成：
\[
\frac{\partial p}{\partial t} = D \frac{\partial^2 p}{\partial x^2} - \frac{\partial (V p)}{\partial x}
\]
其中：扩散项 \( D \frac{\partial^2 p}{\partial x^2} \) 代表无方向偏向的运动。漂移项 \( -\frac{\partial (V p)}{\partial x} \) 代表由于密度变化引起的额外速度贡献。

如果 Run-and-Tumble 粒子的速度 \( v \) 依赖于局部密度 \( \rho(x) \)，即：
\[
v = v(\rho)
\]
那么漂移项 \( V \) 可以表示为：
\[
V = -\frac{1}{\rho} \frac{\partial}{\partial x} \left( v^2(\rho) \rho \right)
\]
这说明：如果 \( v(\rho) \) 随密度 \( \rho \) 下降（即高密度区域速度较低），则可能发生自聚集现象。该方程在高维度 (\( d > 1 \)) 下可以推广，描述更多复杂的细菌集体行为。

### 2.4 Run-and-Tumble 运动的统计力学描述

Run-and-Tumble 运动是一个非平衡过程，其稳态分布不能用玻尔兹曼因子 \( p_s(x) \sim e^{-H/k_BT} \) 描述：由于细菌的运动并非受热噪声控制，无细致平衡（detailed balance）。在无外力情况下，粒子密度 \( p(x) \) 不一定是均匀分布，而可能受速度依赖影响，形成非均匀稳态。

某些情况下，可以将 Run-and-Tumble 体系等效映射到有效热力学系统：设定一个有效自由能：
  \[
  F_{\text{ex}}[\rho] = \int f_{\text{ex}}(\rho) dx
  \]
若存在等效势场 \( V(x) \)，则可以定义一个有效温度：
  \[
  T_{\text{eff}} = \frac{D}{\mu}
  \]
  其中 \( \mu \) 是热力学意义下的摩擦系数。但当相互作用强烈时，这种等效映射会失效。

## 3 多体作用与噪声影响

### 3.1 由单体动力学推广到多体系统

在前面的讨论中，我们研究了单个 Run-and-Tumble 细菌的运动，其行为在大尺度上可以描述为漂移-扩散过程。然而，在真实的细菌群体中，粒子间存在相互作用，例如：排斥作用：由于物理体积排斥，细菌不能无限靠近。碰撞诱导翻滚（collision-induced tumbling）：高密度区域的碰撞增加翻滚率 \( \lambda \)。营养消耗：细菌可能倾向于离开高密度区域，以减少营养竞争。

### 3.2 细菌密度与相互作用对运动的影响

假设细菌的速度 \( v \) 和翻滚率 \( \lambda \) 依赖于局部密度 \( \rho(x) \)，即：
\[
v = v(\rho), \quad \lambda = \lambda(\rho)
\]

在这种情况下，单体的漂移-扩散方程：
\[
\frac{\partial p}{\partial t} = -\frac{\partial J}{\partial x}
\]

其中，概率流修正为：
\[
J = -D(\rho) \frac{\partial p}{\partial x} + V(\rho) p
\]

扩散系数：
\[
D(\rho) = \frac{v^2(\rho)}{\lambda(\rho)}
\]

漂移速度：
\[
V(\rho) = -\frac{1}{\rho} \frac{\partial}{\partial x} \left( D(\rho) \rho \right)
\]

如果 \( v(\rho) \) 下降（例如高密度区域速度降低），则漂移项 \( V(\rho) \) 可能导致细菌在高密度区域“自聚集（self-trapping）”，类似于相分离现象。

#### 相分离条件
设定一个自由能密度：
\[
f_{\text{ex}}(\rho) = \int_{\rho_0}^{\rho} s_1(u) du
\]

其中：
\[
V/D = \frac{d}{dx} \left( \frac{\delta F_{\text{ex}}}{\delta \rho} \right)
\]

如果自由能密度的二阶导数：
\[
\frac{d^2 f_{\text{ex}}}{d\rho^2} < 0
\]
则系统在密度涨落下可能进入相分离状态，即局部高密度区域和低密度区域共存。

### 3.3 统计力学中的噪声项推导

在多体系统中，细菌密度 \( \rho(x,t) \) 是一个随机变量，其演化方程可以用Langevin 形式描述：
\[
\frac{\partial \rho}{\partial t} = -\frac{\partial J}{\partial x} + \eta(x,t)
\]

其中噪声项 \( \eta(x,t) \) 表示密度涨落，它的均值为零，且满足高斯白噪声性质：
\[
\langle \eta(x,t) \eta(x',t') \rangle = 2 D(\rho) \delta(x-x') \delta(t-t')
\]

进一步，我们可以推导出密度演化的 Fokker-Planck 方程：
\[
\frac{\partial P[\rho]}{\partial t} = -\int dx \frac{\delta}{\delta \rho(x)} \left[ \left( -\frac{\partial}{\partial x} J \right) P[\rho] \right] + \frac{1}{2} \int dx dx' \frac{\delta^2}{\delta \rho(x) \delta \rho(x')} \left[ \langle \eta(x) \eta(x') \rangle P[\rho] \right]
\]

这表明：噪声项 \( \eta(x,t) \) 由扩散系数 \( D(\rho) \) 调制，高密度区域的噪声可能更强。该系统类似于一个有效的非平衡涨落动力学系统。

#### Ito-Langevin 过程
密度场 \( \rho(x,t) \) 也可以用 Langevin 过程描述：
\[
\frac{dx_i}{dt} = A(x_i) + C(x_i) \xi_i(t)
\]

其中：
- 漂移项 \( A(x_i) = V(\rho) + \partial_x D(\rho) \)
- 噪声项 \( C(x_i) = \sqrt{2 D(\rho)} \)
- 白噪声 \( \xi_i(t) \) 满足：
  \[
  \langle \xi_i(t) \xi_j(t') \rangle = \delta_{ij} \delta(t-t')
  \]

即使没有外部驱动力，相互作用本身也会导致密度涨落。


### 3.4 Run-and-Tumble 细菌系统的集体行为

多体相互作用导致 Run-and-Tumble 细菌系统出现多种集体现象：

（1）自聚集（Self-trapping）

当 \( v(\rho) \) 下降（即 \( dv/d\rho < 0 \)）时，粒子可能被困在高密度区域，形成稳定的“团簇”。这种现象类似于热平衡系统中的相分离，但它是由于运动学效应（kinetic effect）而非热力学驱动。

（2）相分离（Phase Separation）

在适当的密度范围内，Run-and-Tumble 细菌可分离为高密度相和低密度相。自由能分析表明，当 \( d^2 f_{\text{ex}} / d\rho^2 < 0 \) 时，系统局部不稳定。计算机模拟显示，Run-and-Tumble 细菌可形成明显的高密度和低密度区，类似于液-气相分离。

（3）有效温度控制的动力学

在弱相互作用情况下，系统可以用一个有效温度 \( T_{\text{eff}} \) 描述：
\[
T_{\text{eff}} = \frac{D}{\mu}
\]
但在强相互作用下，该映射失效，因为系统不再满足平衡态的涨落-耗散定理。

（4）外场影响

如果细菌处于外场（如重力场），漂移项 \( V(x) \) 受到修正：
\[
V(x) = -\frac{1}{\rho} \frac{\partial}{\partial x} \left( D(\rho) \rho \right) - mg
\]
当 \( mg \) 足够大时，系统密度可能呈指数衰减，类似于沉降平衡（sedimentation equilibrium）。但由于 Run-and-Tumble 运动不同于布朗运动，密度衰减曲线可能偏离传统玻尔兹曼分布。

## 4 总结
### 4.1 单个粒子的运动方程
#### 4.1.1 Run-and-Tumble 运动的基本描述

在一维情况下，我们定义：\( R(x,t) \)：时间 \( t \) 处于 \( x \) 位置并向右运动的粒子概率密度。\( L(x,t) \)：时间 \( t \) 处于 \( x \) 位置并向左运动的粒子概率密度。Run-and-Tumble 细菌的运动包括两个过程：以速度 \( \pm v \) 直线运动。以速率 \( \lambda \) 翻滚，改变方向。

根据此模型，\( R(x,t) \) 和 \( L(x,t) \) 满足：
\[
\frac{\partial R}{\partial t} + v \frac{\partial R}{\partial x} = -\frac{\lambda}{2} R + \frac{\lambda}{2} L
\]
\[
\frac{\partial L}{\partial t} - v \frac{\partial L}{\partial x} = -\frac{\lambda}{2} L + \frac{\lambda}{2} R
\]
对流项：\( \pm v \frac{\partial R}{\partial x} \) 表示细菌沿 \( x \) 方向移动。翻滚项：\( -\lambda R/2 + \lambda L/2 \) 表示细菌从右移变成左移，反之亦然。

#### 4.1.2 概率流方程
定义总概率密度：
\[
p(x,t) = R(x,t) + L(x,t)
\]
定义概率流：
\[
J(x,t) = v R(x,t) - v L(x,t)
\]
求和可得：
\[
\frac{\partial p}{\partial t} = -\frac{\partial J}{\partial x}
\]
这表明系统满足概率守恒，即密度随时间变化由概率流 \( J \) 决定。

求差可得：
\[
\frac{\partial J}{\partial t} + v \frac{\partial p}{\partial x} = -\lambda J
\]
该方程描述了概率流的时间演化，其中 \( -\lambda J \) 表示翻滚导致的流衰减。

### 4.2 扩散极限

#### 4.2.1 Langevin 方程的构造
Run-and-Tumble 运动可以用 Langevin 形式表示：
\[
\frac{dx}{dt} = v \sigma(t)
\]
其中：\( \sigma(t) \) 为二值随机变量，取值 \( \pm 1 \)。\( \sigma(t) \) 服从泊松过程，以速率 \( \lambda \) 发生翻滚：
  \[
  P(\sigma \to -\sigma) = \lambda dt
  \]

在长时间尺度下，\(\sigma(t)\) 近似为Ornstein-Uhlenbeck 过程：
\[
\frac{d\sigma}{dt} = -\lambda \sigma + \eta(t)
\]
其中 \( \eta(t) \) 是均值为零、方差为 \( 2\lambda \) 的白噪声。

#### 4.2.2 漂移扩散方程的推导
在长时间尺度 \( t \gg \lambda^{-1} \) 下，假设概率流达到准稳态：
\[
\frac{\partial J}{\partial t} \approx 0
\]
解得：
\[
J \approx -\frac{v^2}{\lambda} \frac{\partial p}{\partial x}
\]
代入概率守恒方程：
\[
\frac{\partial p}{\partial t} = \frac{v^2}{\lambda} \frac{\partial^2 p}{\partial x^2}
\]
我们得到了扩散方程：
\[
\frac{\partial p}{\partial t} = D \frac{\partial^2 p}{\partial x^2}, \quad D = \frac{v^2}{\lambda}
\]
该方程表明：在长时间尺度上，Run-and-Tumble 运动表现为扩散过程。扩散系数 \( D \) 由速度 \( v \) 和翻滚率 \( \lambda \) 决定。

### 4.3 多体作用粒子的统计力学描述

#### 4.3.1 由单粒子描述推广到密度场方程
如果粒子间存在相互作用，例如速度 \( v \) 依赖于局部密度 \( \rho(x) \)：
\[
v = v(\rho)
\]
则总密度满足：
\[
\frac{\partial \rho}{\partial t} = -\frac{\partial}{\partial x} J
\]
其中：
\[
J = -D(\rho) \frac{\partial \rho}{\partial x} + V(\rho) \rho
\]

#### 4.3.2 Ito-Langevin 过程推导
密度场的随机动力学可写为：
\[
\frac{\partial \rho}{\partial t} = -\frac{\partial}{\partial x} \left( A(\rho) \rho - D(\rho) \frac{\partial \rho}{\partial x} + \eta(x,t) \right)
\]
其中：漂移项：
  \[
  A(\rho) = V(\rho) + \frac{\partial D}{\partial x}
  \]
噪声项：
  \[
  \langle \eta(x,t) \eta(x',t') \rangle = 2D(\rho) \delta(x-x') \delta(t-t')
  \]

#### 4.3.3 Fokker-Planck 方程的建立
密度分布的概率满足 Fokker-Planck 方程：
\[
\frac{\partial P[\rho]}{\partial t} = -\int dx \frac{\delta}{\delta \rho(x)} \left[ \left( -\frac{\partial}{\partial x} J \right) P[\rho] \right] + \frac{1}{2} \int dx dx' \frac{\delta^2}{\delta \rho(x) \delta \rho(x')} \left[ \langle \eta(x) \eta(x') \rangle P[\rho] \right]
\]

### 4.4 详细平衡条件与有效自由能

#### 4.4.1 何种情况下系统可以映射到热力学平衡
如果 Run-and-Tumble 体系可以等效映射到一个热平衡系统，则必须满足：
\[
V(\rho) / D(\rho) = -\frac{d}{dx} \left( \frac{\delta F_{\text{ex}}}{\delta \rho} \right)
\]
其中有效自由能：
\[
F_{\text{ex}}[\rho] = \int f_{\text{ex}}(\rho) dx
\]

如果 \( f_{\text{ex}}(\rho) \) 的二阶导数：
\[
\frac{d^2 f_{\text{ex}}}{d\rho^2} < 0
\]
则系统可能发生自聚集或相分离。

#### 4.4.2 过渡到高维情况下的修正
在高维 (\( d > 1 \))，Run-and-Tumble 动力学可推广为：
\[
\frac{\partial p}{\partial t} = D \nabla^2 p - \nabla \cdot (V p)
\]
如果粒子之间存在相互作用，则相分离行为可以用 Cahn-Hilliard 方程描述：
\[
\frac{\partial \rho}{\partial t} = \nabla \cdot (M(\rho) \nabla \mu)
\]
其中 \( \mu = \delta F_{\text{ex}} / \delta \rho \) 是化学势。
