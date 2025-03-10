---
title: 重整化群流最优传输与扩散模型
date: 2024-12-17
math: true
---

关于重整化流、最优传输与扩散模型的笔记，一个由重整化群思想启发的扩散模型设计和优化。

---

## 1 背景

### 1.1 重整化群（Renormalization Group, RG）

重整化群（RG）是理论物理中的一个重要框架，用于研究系统在不同尺度上的行为。在统计物理和量子场论中，RG描述了如何在不同的能量或长度尺度上调整系统的自由度，使其行为保持一致。

RG的两种主要形式

Wilson 重整化群：主要思想是积分出短程自由度，只保留低能有效自由度。通过逐步消除高频自由度，得到系统在更大尺度下的有效描述。主要用于统计力学和凝聚态物理，例如：Ising 模型的相变，超导和量子相变。

Exact Renormalization Group (ERG, 精确重整化群)以连续方程的形式描述RG流，如 Polchinski 方程。通过一个微分方程，描述系统的有效作用量如何随着尺度变化。


RG 解释了为什么在临界点附近不同物理系统具有相同的临界指数。可以用于分析相变，例如铁磁相变、液-气临界点等。在量子场论中，RG 用于处理发散问题，如 QED 和 QCD 的重整化。

### 1.2 最优传输（Optimal Transport）
最优传输（Optimal Transport, OT）源于 Gaspard Monge (1781) 的研究，目标是在给定的运输成本下，寻找将一个质量分布转换到另一个质量分布的最优方式。Kantorovich (1942) 对其进行了推广，提出了一种基于线性规划的方法，使得问题更容易求解。

1. Monge 形式：给定两个概率分布 \(\mu_X, \mu_Y\)，寻找一个映射 \(T: X \to Y\)，使得质量守恒：
     \[
     T_* \mu_X = \mu_Y
     \]
目标是最小化运输成本：
     \[
     \inf_T \int_X c(x, T(x)) d\mu_X
     \]

Kantorovich 形式：允许质量在多个目标点之间分配，引入联合分布 \(\pi(x, y)\)，最小化成本：
     \[
     \inf_{\pi \in \Gamma(\mu_X, \mu_Y)} \int_{X \times Y} c(x, y) d\pi(x, y)
     \]

Wasserstein 距离：最优传输的一个核心概念是 Wasserstein 距离（OT-距离），它测量两个概率分布之间的“最小工作量”：
     \[
     W_p(\mu_X, \mu_Y) = \left( \inf_{\pi} \int |x - y|^p d\pi(x, y) \right)^{\frac{1}{p}}
     \]
其中 \(W_2\) Wasserstein 度量在统计物理和机器学习中应用广泛。

Wasserstein 距离可用于刻画系统的熵变化，RG 流可视为一个熵增加过程。最优传输可以用于描述系统的相关性如何随着尺度变化。OT 框架下的 Benamou-Brenier 公式描述了扩散过程。

- 生成对抗网络（GANs）：传统GAN 使用 Jensen-Shannon 散度，而 Wasserstein GAN（WGAN）利用 Wasserstein 距离，提高了生成质量。
- 图神经网络（GNNs）：在图匹配、节点嵌入等任务中，OT 提供了一种测量不同图结构的方法。
- 扩散模型（Diffusion Models）：通过 Wasserstein 距离优化扩散过程，提高生成质量。

### 1.3 生成式模型（Generative Models）与扩散模型（Diffusion Models）
- 概率生成模型：变分自编码器（VAE）玻尔兹曼机（Boltzmann Machines）
- 非概率生成模型：GANs，基于自回归的方法（如 PixelCNN）
- 基于扩散的生成模型（Diffusion Models）：受热力学扩散过程启发，采用逐步添加噪声并逆向去噪的方式进行生成。

## 2 重整化群（RG）与最优传输
### 2.1 重整化群（RG）基础

重整化群（Renormalization Group, RG） 研究物理系统在不同尺度上的行为变化，广泛用于统计物理、凝聚态物理和量子场论。在RG框架下，一个物理系统在不同能量尺度上的行为可以通过“流动”描述，即随着尺度变化，系统的参数（如耦合常数）会发生演化。基本思想：去除短程自由度（高能自由度），获得有效的低能理论。

### 2.2 Wilson 重整化群 vs. Exact Renormalization Group (ERG)

#### (1) Wilson 重整化群
通过积分出高能自由度来获得低能有效理论：设系统的分区函数为：
    \[
    Z = \int \mathcal{D} \phi \, e^{-S[\phi]}
    \]
通过积分掉短程自由度（如高动量模式 \(\phi_{>}\)），得到新的有效作用量 \(S'[\phi_{<}]\)：
    \[
    e^{-S'[\phi_{<}]} = \int \mathcal{D} \phi_{>} e^{-S[\phi_{<}, \phi_{>}]}
    \]
这样，低能自由度的作用量发生变化，从而定义了RG流。

#### (2) Exact Renormalization Group (ERG)
以连续方程描述RG流：Polchinski方程：
    \[
    \Lambda \frac{d S_{\Lambda}}{d\Lambda} = \frac{1}{2} \int d^d k \, \frac{d R_{\Lambda}(k)}{d\Lambda} \left[ \frac{\delta S_{\Lambda}}{\delta \phi_k} \frac{\delta S_{\Lambda}}{\delta \phi_{-k}} - \frac{\delta^2 S_{\Lambda}}{\delta \phi_k \delta \phi_{-k}} \right]
    \]
Wegner-Morris流方程：
    \[
    \frac{\partial S_k}{\partial k} = \mathcal{F}[S_k]
    \]
其中 \(\mathcal{F}\) 是一个微分算子，描述了作用量随尺度变化的流动。

### 2.3 RG流的单调量（RG Monotones）
定义：一个单调量（Monotone）是在RG流中单调变化的量，通常用于描述熵或信息损失。新提出的RG单调量：
  \[
  \mathcal{M}(t) = \int p(x) \log p(x) dx
  \]
其中 \(p(x)\) 是场论中的概率密度函数。该量在RG流中单调减少，类似于信息熵的增长。

### 2.4. 最优传输问题
目标：寻找将一个概率分布转换到另一个概率分布的最优方式。形式化描述：设 \(X, Y\) 为测度空间，\(\mu_X, \mu_Y\) 为两个概率分布，目标是寻找最优的传输方式。

#### (1) Monge 公式
Monge 公式定义了一个映射 \(T: X \to Y\)，满足：
  \[
  T_* \mu_X = \mu_Y
  \]
目标是最小化传输成本：
  \[
  \inf_T \int_X c(x, T(x)) d\mu_X
  \]
其中 \(c(x, y)\) 是运输成本。

#### (2) Kantorovich 公式
允许概率质量分布在多个目标点之间，使用联合分布 \(\pi(x, y)\)：
  \[
  \inf_{\pi \in \Gamma(\mu_X, \mu_Y)} \int_{X \times Y} c(x, y) d\pi(x, y)
  \]
其中 \(\Gamma(\mu_X, \mu_Y)\) 是所有满足边际约束的联合分布：
    \[
    \int_Y \pi(x, y) dy = \mu_X(x), \quad \int_X \pi(x, y) dx = \mu_Y(y)
    \]

### 2.5 Wasserstein度量
#### (1) Wasserstein-2 度量
Wasserstein-2 度量（OT-距离）定义为：
  \[
  W_2(\mu_X, \mu_Y) = \left( \inf_{\pi} \int |x - y|^2 d\pi(x, y) \right)^{\frac{1}{2}}
  \]
  - 描述了从 \(\mu_X\) 到 \(\mu_Y\) 进行最优传输的“最小工作量”。
  - 在许多物理系统中，该度量可以用于描述相变或分布演化。

#### (2) Otto计算（Otto Calculus）与信息熵
Wasserstein空间可以定义几何结构，使其成为黎曼流形，并引入“Wasserstein梯度流”。设自由能泛函为：
  \[
  \mathcal{F}(\mu) = \int V(x) d\mu(x) + \int d\mu(x) \log d\mu(x)
  \]
其 Wasserstein 梯度流为：
  \[
  \frac{\partial \mu}{\partial t} = \nabla_W \mathcal{F}(\mu)
  \]
这类方程描述了熵的增长过程，与RG流有直接联系。

### 2.6 Benamou-Brenier 公式
该公式提供了一种将最优传输问题转化为流体力学问题的方法。最优传输距离可以通过求解如下优化问题得到：
  \[
  W_2^2(\mu_0, \mu_1) = \inf_{\rho, v} \int_0^1 \int_{\mathbb{R}^d} \rho(x, t) |v(x, t)|^2 dx dt
  \]
  其中：\(\rho(x, t)\) 表示时间 \(t\) 时刻的概率密度。\(v(x, t)\) 表示速度场，满足连续性方程：
    \[
    \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho v) = 0
    \]
该变分问题的解描述了从初始分布 \(\mu_0\) 到最终分布 \(\mu_1\) 的最优流动。



## 3. RG流的最优传输解释
### 3.1 证明Polchinski方程可以重写为最优传输梯度流
Polchinski 方程描述了作用量 \( S_\Lambda[\phi] \) 在不同尺度 \(\Lambda\) 下的演化：
\[
\Lambda \frac{d S_{\Lambda}}{d\Lambda} = \frac{1}{2} \int d^d k \, \frac{d R_{\Lambda}(k)}{d\Lambda} \left[ \frac{\delta S_{\Lambda}}{\delta \phi_k} \frac{\delta S_{\Lambda}}{\delta \phi_{-k}} - \frac{\delta^2 S_{\Lambda}}{\delta \phi_k \delta \phi_{-k}} \right]
\]
其中：\( R_{\Lambda}(k) \) 是一个调节函数（cutoff function），用于控制不同动量模式的贡献。第一项描述了作用量对场的梯度项，第二项对应于涨落修正。这一方程表明，在RG流中，短程自由度逐步被积分掉，并导致作用量的演化。

在最优传输理论中，Wasserstein-2 度量用于衡量两个概率分布之间的最优传输距离：
\[
W_2^2(\mu_0, \mu_1) = \inf_{\rho, v} \int_0^1 \int_{\mathbb{R}^d} \rho(x, t) |v(x, t)|^2 dx dt
\]
其中：\(\rho(x,t)\) 为时间 \( t \) 时刻的概率密度函数。\( v(x,t) \) 为速度场，满足连续性方程：
  \[
  \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho v) = 0
  \]
该方程确保质量守恒。为了将RG流重写为Wasserstein梯度流，我们需要找到适当的自由能泛函 \(\mathcal{F}(\mu)\)，使得其梯度流能描述RG流。

我们引入一个自由能泛函：
\[
\mathcal{F}(\mu) = \int \left[ \frac{1}{2} \phi (-\nabla^2) \phi + V(\phi) \right] d^dx
\]
其中：第一项表示自由场的动能项（kinetic term）。第二项 \(V(\phi)\) 代表相互作用势能。Wasserstein梯度流的通用形式：
\[
\frac{\partial \mu}{\partial t} = - \nabla_W \mathcal{F}(\mu)
\]
其中 \(\nabla_W\) 代表Wasserstein空间上的梯度运算。

如果我们选择合适的度量结构 \( g^{-1} = \frac{\partial \Sigma}{\partial t} \)，并设定：
\[
\frac{\partial S}{\partial t} = \nabla_W \mathcal{F}(\mu)
\]
那么RG流可以重写为最优传输梯度流的形式。

这一变换说明：
1. RG流是一个熵增加过程，与最优传输中的熵增长原理一致。
2. RG流的不可逆性可以通过Wasserstein度量来度量。

### 3.2 信息熵的增长与RG流的不逆性
在RG流中，高能自由度逐步被消去，导致信息熵增加：
\[
S(\mu) = - \int \mu(x) \log \mu(x) dx
\]
根据RG流的单调性，信息熵在RG流中是单调增加的：
\[
\frac{dS}{d\Lambda} \geq 0
\]
这与热力学第二定律类似，表明RG流是一个不可逆过程。进一步，RG单调量（RG monotone）可以通过相对熵（Kullback-Leibler 散度）定义：
\[
\mathcal{M}(\mu) = \int \mu(x) \log \frac{\mu(x)}{\mu_*(x)} dx
\]
其中 \(\mu_*\) 是一个固定点分布。在RG流中，该单调量随着\(\Lambda\) 的演化而单调减少：
\[
\frac{d\mathcal{M}}{d\Lambda} \leq 0
\]
这表明RG流在信息度量下是单向演化的，进一步支持其最优传输梯度流的解释。

### 3.3 计算自由场的RG流：连接RG单调量与信息熵
考虑自由标量场的作用量：
\[
S = \int d^dx \left[ \frac{1}{2} (\nabla \phi)^2 + \frac{m^2}{2} \phi^2 \right]
\]
其分区函数为：
\[
Z = \int \mathcal{D} \phi \, e^{-S}
\]
RG流的单调量可以计算为：
\[
\mathcal{M} = \int d^d k \, \log \left( \frac{k^2 + m^2}{\Lambda^2 + m^2} \right)
\]
随着\(\Lambda\) 下降，该单调量单调减少。说明自由场的RG流确实对应于信息熵的增长过程。Wasserstein度量可以用于测量RG流的熵变。

## 4. 扩散模型
### 4.1 扩散生成模型（Diffusion Models）基础

扩散模型（Diffusion Models）是一类基于马尔可夫过程的生成模型，其核心思想是：
1. 前向扩散（Forward Process）：逐步向数据添加噪声，使其变为一个标准高斯分布。
2. 逆向去噪（Reverse Process）：通过学习逆向扩散过程，从高斯噪声中逐步恢复原始数据。

逐步添加噪声（Forward Process）

前向扩散是一个离散时间马尔可夫过程，从数据分布 \( x_0 \sim q(x_0) \) 开始，逐步添加噪声：
\[
q(x_t | x_{t-1}) = \mathcal{N}(x_t; \sqrt{1 - \beta_t} x_{t-1}, \beta_t I)
\]
其中：\( x_t \) 是第 \( t \) 步的噪声版本。\( \beta_t \) 是一个时间相关的噪声调度参数，通常是递增的。通过重复应用该过程，数据最终收敛到一个各向同性的高斯分布。

前向扩散的 闭式表达式：
\[
q(x_t | x_0) = \mathcal{N}(x_t; \sqrt{\bar{\alpha}_t} x_0, (1 - \bar{\alpha}_t) I)
\]
其中：\( \alpha_t = 1 - \beta_t \)\( \bar{\alpha}_t = \prod_{s=1}^{t} \alpha_s \)

这表明，可以直接在第 \( t \) 步从数据 \( x_0 \) 生成 \( x_t \)：
\[
x_t = \sqrt{\bar{\alpha}_t} x_0 + \sqrt{1 - \bar{\alpha}_t} \epsilon, \quad \epsilon \sim \mathcal{N}(0, I)
\]

逆向去噪（Reverse Process）

逆向扩散过程旨在学习从噪声恢复数据的马尔可夫链：
\[
p_\theta(x_{t-1} | x_t) = \mathcal{N}(x_{t-1}; \mu_\theta(x_t, t), \Sigma_\theta(x_t, t))
\]
其中：\( \mu_\theta(x_t, t) \) 是神经网络（通常为 U-Net）预测的去噪均值。\( \Sigma_\theta(x_t, t) \) 是协方差矩阵，通常设为对角矩阵。

一种常见的估计方法是预测噪声：
\[
\epsilon_\theta(x_t, t) = \frac{x_t - \sqrt{\bar{\alpha}_t} x_0}{\sqrt{1 - \bar{\alpha}_t}}
\]
则去噪均值可由以下公式估计：
\[
\mu_\theta(x_t, t) = \frac{1}{\sqrt{\alpha_t}} \left( x_t - \frac{\beta_t}{\sqrt{1 - \bar{\alpha}_t}} \epsilon_\theta(x_t, t) \right)
\]

### 4.2 连接信息论与扩散模型
Kullback-Leibler 散度（KL散度）

扩散模型训练的核心是最大化数据对数似然：
\[
\log p_\theta(x_0)
\]
但直接优化它通常很困难，因此采用变分下界（Variational Lower Bound, VLB） 进行优化：
\[
\mathbb{E}_{q(x_{1:T} | x_0)} \left[ \sum_{t=1}^{T} D_{KL}(q(x_t | x_{t-1}) || p_\theta(x_{t-1} | x_t)) - \log p(x_T) \right]
\]
其中：KL 散度 \( D_{KL} \) 衡量模型分布与真实数据分布之间的差异。

变分推断（Variational Inference）

VLB 目标可拆解为：
\[
L = L_T + \sum_{t=2}^{T} L_t + L_1
\]
其中：\( L_T = D_{KL}(q(x_T | x_0) || p(x_T)) \) 使得最终噪声逼近标准正态分布。\( L_t = D_{KL}(q(x_t | x_{t-1}) || p_\theta(x_{t-1} | x_t)) \) 训练模型学习逆扩散分布。\( L_1 \) 是数据重构损失项。

优化这个目标可训练扩散模型，使其学习如何从高斯噪声恢复数据。

### 4.3 频域扩散模型（FDDM）

FDDM 通过离散余弦变换（DCT）将数据转换到频域：
\[
\tilde{x} = F(x), \quad x = F^{-1}(\tilde{x})
\]
其中 \( F \) 表示 DCT 变换，\( F^{-1} \) 为逆变换。

在频域上执行扩散：
\[
q(\tilde{x}_t | \tilde{x}_{t-1}) = \mathcal{N}(\tilde{x}_t; \sqrt{1 - \beta_t} \tilde{x}_{t-1}, \beta_t I)
\]
与时域扩散类似，但能更有效地区分信号和噪声。低频部分对应全局结构，高频部分对应细节。通过分阶段去噪，可以优先恢复全局结构，再修复细节，提高生成质量。在频域中，噪声可以更高效地去除，使得模型训练更加稳定。计算复杂度更低，适用于大规模数据生成。

#### 前向扩散（Forward Diffusion）
FDDM 在频域执行前向扩散：
\[
q(\tilde{x}_t | \tilde{x}_0) = \mathcal{N}(\tilde{x}_t; \sqrt{\bar{\alpha}_t} \tilde{x}_0, (1 - \bar{\alpha}_t) I)
\]
其中 \( \tilde{x} \) 为 DCT 变换后的数据。

#### 逆向扩散（Reverse Diffusion）
在逆向过程中，FDDM 估计去噪步骤：
\[
p_\theta(\tilde{x}_{t-1} | \tilde{x}_t) = \mathcal{N}(\tilde{x}_{t-1}; \mu_\theta(\tilde{x}_t, t), \Sigma_\theta)
\]
其中：
\[
\mu_\theta(\tilde{x}_t, t) = \frac{1}{\sqrt{\alpha_t}} \left( \tilde{x}_t - \frac{\beta_t}{\sqrt{1 - \bar{\alpha}_t}} \epsilon_\theta(\tilde{x}_t, t) \right)
\]

#### 训练目标与损失函数
训练目标是最小化预测噪声的均方误差（MSE）：
\[
L_\theta = \mathbb{E}_{x_0, \epsilon} \left[ ||\epsilon - \epsilon_\theta(F(x_t), t)||^2 \right]
\]
其中：\( \epsilon \) 是真实噪声。\( \epsilon_\theta \) 是模型预测的噪声。这种损失函数确保模型学习到最佳的去噪策略，从而实现高质量的生成。

## 5. RG流、Wasserstein梯度流与扩散模型的联系
### 5.1 证明扩散模型的反向去噪过程等价于RG的逆向流

重整化群（RG）描述了系统在不同尺度上的演化。RG流的基本方程（Polchinski方程）为：
\[
\Lambda \frac{d S_{\Lambda}}{d\Lambda} = \frac{1}{2} \int d^d k \, \frac{d R_{\Lambda}(k)}{d\Lambda} \left[ \frac{\delta S_{\Lambda}}{\delta \phi_k} \frac{\delta S_{\Lambda}}{\delta \phi_{-k}} - \frac{\delta^2 S_{\Lambda}}{\delta \phi_k \delta \phi_{-k}} \right]
\]
其中：\(\Lambda\) 是 RG 变换的尺度参数。\( S_{\Lambda} \) 是作用量，描述了低能有效理论。第一项表示作用量的非线性贡献，第二项表示涨落修正。

RG 逆向流的思想：在 RG 过程中，低能有效理论逐步丢失高能自由度，相当于增加熵。逆向 RG 流（Inverse RG Flow）可以视为一个去噪过程，试图从低能有效理论恢复更精细的细节。

扩散模型的逆向去噪方程：
\[
p_\theta(x_{t-1} | x_t) = \mathcal{N}(x_{t-1}; \mu_\theta(x_t, t), \Sigma_\theta(x_t, t))
\]
其中：
\[
\mu_\theta(x_t, t) = \frac{1}{\sqrt{\alpha_t}} \left( x_t - \frac{\beta_t}{\sqrt{1 - \bar{\alpha}_t}} \epsilon_\theta(x_t, t) \right)
\]
这个公式描述了如何从高斯噪声逐步恢复原始数据。

两者的相似性：
- 去除噪声 vs. 恢复高能自由度：逆向扩散试图去除高斯噪声，使数据回归到真实分布。逆向 RG 试图恢复被丢弃的高能自由度，重建原始理论。

- 梯度流形式：RG 逆向流可写为 Wasserstein 梯度流：
     \[
     \frac{\partial S}{\partial t} = \nabla_W \mathcal{F}(\mu)
     \]
逆向扩散的均值计算公式也可以表示为 Wasserstein 梯度下降：
     \[
     x_{t-1} = x_t - \eta \nabla_W \mathcal{F}(\mu)
     \]
这表明 去噪过程等价于最优传输梯度下降，即扩散模型的逆向过程等价于 RG 的逆向流。

### 5.2 Wasserstein梯度流如何影响扩散模型
Wasserstein梯度流的定义:Wasserstein-2 流是一类优化问题：
\[
\frac{\partial \mu}{\partial t} = - \nabla_W \mathcal{F}(\mu)
\]
其中自由能泛函 \(\mathcal{F}(\mu)\) 由数据分布的熵与势能项构成：
\[
\mathcal{F}(\mu) = \int \left[ \frac{1}{2} |\nabla \phi|^2 + V(\phi) \right] d^dx
\]

扩散模型的优化目标是最小化 Wasserstein-2 变分问题：
\[
\min_\theta \mathbb{E}_{x_0, \epsilon} \left[ \| \epsilon - \epsilon_\theta(x_t, t) \|^2 \right]
\]
这相当于沿 Wasserstein 路径执行梯度下降，以最优方式从噪声恢复数据。

FDDM的噪声调度（Noise Schedule）设计: FDDM 采用基于 Wasserstein 流的噪声调度方案，在频域定义：
\[
\tilde{x}_t = \sqrt{\bar{\alpha}_t} \tilde{x}_0 + \sqrt{1 - \bar{\alpha}_t} \epsilon
\]
噪声调度函数 \(\alpha_t\) 采用动量依赖形式：
\[
\alpha_t(k) = \frac{1}{e^{(\epsilon_k - \mu_t)/T} + 1}
\]
其中：\(\epsilon_k\) 是频域动量分布。\(\mu_t\) 控制信噪比（SNR）。\(T\) 控制噪声扩散的平滑程度。在低频区域（\(k\) 小）保持较高的 \(\alpha_t\)，确保全局结构先恢复。在高频区域（\(k\) 大）降低 \(\alpha_t\)，确保去噪过程更平稳。

在扩散过程中，不同频率成分具有不同的噪声扩散模式：低频部分（大尺度结构） 受到较少的噪声污染，先进行去噪。高频部分（细节信息） 受到较强的噪声污染，后进行恢复。通过 Wasserstein 梯度流优化的噪声调度，可以使得低频部分先恢复，提高整体生成质量。

## 参考文献
1. Renormalization Group Flow as Optimal Transport
   
2. Renormalization Group Flow, Optimal Transport and Diffusion-Based Generative Model