---
title: 生成模型的几何结构
date: 2024-10-20
---

关于生成模型最近又看到一篇理论性方面进行研究的文章 Deep Generative Models through the Lens of the Manifold Hypothesis，这篇文章主要探讨了深度生成模型（DGMs）在流形假设下的几何结构，并分析了不同生成模型在流形上的表现和局限性。

## 1 流形学习（Manifold Learning）

### 1.1 高维数据往往存在低维结构  

在许多实际应用中，数据通常具有高维特性。例如：图像数据集（如ImageNet、MNIST）的像素维度通常非常高，但图像本身可能由较少的关键因素决定，如光照、视角、物体类别等。语音信号虽然是时间序列数据，包含大量采样点，但其核心特征（如音素、语调）可能受低维因素控制。蛋白质三维结构虽然存在于高维空间（如PDB格式表示的三维坐标），但其形态通常由较低维的生物化学约束决定。  

流形假设（Manifold Hypothesis） 认为：高维数据通常分布在一个比原始维度低得多的流形（Manifold）上，而不是均匀填充整个高维空间。这一假设在多个领域得到了实验验证，例如：  

- 研究发现，图像数据的内在维度（Intrinsic Dimension）远小于像素空间维度。例如，在CIFAR-10数据集中，图像的内在维度可能只有几十到几百，而不是数千个像素点的维度。  
- 语音数据的低维特征可以通过降维技术提取，如主成分分析（PCA）、自编码器等。  
- 在物理系统（如分子动力学模拟）中，粒子的运动受到低维约束，而不是随机填充整个高维状态空间。  

### 1.2 常见的流形学习方法  

由于数据分布在低维流形上，许多机器学习方法试图学习这些低维结构。常见的方法包括：  

谱方法（Spectral Methods）：线性降维：主成分分析（PCA）、多维缩放（MDS）。非线性降维：t-SNE（t-distributed Stochastic Neighbor Embedding）、Isomap、LLE（Locally Linear Embedding）。  

自编码器（Autoencoders, AEs）通过神经网络构造一个编码器（Encoder）和解码器（Decoder），将高维数据映射到低维潜在空间（Latent Space），再从潜在空间重构数据。变种包括变分自编码器（VAE）、对抗自编码器（AAE）等。  

流形学习算法（Manifold Learning Algorithms）旨在学习数据的流形结构，如流形正则化（Manifold Regularization）、流形流（Manifold Flows）。在生成模型中，一些方法尝试通过学习流形结构来提高生成数据的质量，如基于流形的扩散模型。  

流形学习的核心目标是获取数据的低维本质结构，而不是仅仅关注高维空间中的数据分布。  

### 1.3 深度生成模型（Deep Generative Models）  

深度生成模型（DGMs）利用神经网络从数据中学习概率分布，并生成新样本。核心目标是建模真实数据分布 \( p^*(x) \) ，并从中采样新的数据点 \( x \sim p_{\theta}(x) \) 。  

变分自编码器（Variational Autoencoders, VAE）通过最大化证据下界（ELBO），近似真实数据的后验分布。使用一个潜在变量 \( z \) ，假设数据由低维流形上的变量生成，并通过解码器 \( p(x|z) \) 进行数据重建。由于采用高斯分布假设，VAE 容易产生模糊样本，不适用于高精度图像生成。虽然学习了潜在流形，但其高斯假设限制了其表达能力。

归一化流（Normalizing Flows, NF）通过可逆变换构建显式概率密度模型，使数据从复杂分布映射到简单分布（如高斯分布）。训练时利用变换的雅可比行列式（Jacobian Determinant）计算精确的对数似然。由于需要计算高维矩阵行列式，NF 在高维数据上训练成本较高。无法直接建模低维流形上的数据，而是倾向于覆盖整个高维空间。  

生成对抗网络（Generative Adversarial Networks, GANs）采用博弈论框架，训练一个生成器（Generator）和判别器（Discriminator）。标是让生成器学习真实数据分布，并欺骗判别器，使其无法区分真实数据和生成数据。训练时存在模式崩溃（Mode Collapse）、不稳定性等问题。由于采用隐式分布，往往能够学习流形上的数据分布，但训练不稳定。

扩散模型（Diffusion Models, DM）通过向数据添加噪声并学习去噪过程来生成数据。训练过程类似于马尔可夫链，最终可以生成高质量的样本。目前在高分辨率图像合成（如DALL-E、Stable Diffusion）方面表现突出。通过去噪过程逐步逼近数据流形，适用于复杂数据分布的学习。

大多数生成模型未能显式考虑数据的流形结构，导致数值不稳定或生成质量下降。  

### 1.4 常用的度量方法  

KL散度是计算两个概率分布 \( P \) 和 \( Q \) 之间差异的常见方法，其定义为：  
\[
KL(P \| Q) = \int P(x) \log \frac{P(x)}{Q(x)} dx
\]  
然而，在流形学习中，KL 散度存在两个主要问题：1. KL散度可能趋于无穷大：如果 \( Q(x) \) 的支撑集（support）不包含 \( P(x) \) 的支撑集，KL散度会变成无穷大，导致优化问题不稳定。2. KL散度对概率密度的绝对值敏感：这意味着即使 \( P \) 和 \( Q \) 仅在流形上不同，KL 散度可能会过度惩罚微小的差异，而不考虑整体结构的相似性。  

Wasserstein 距离（Wasserstein Distance, WD）是一种更适合流形学习的度量，它衡量了两个分布之间的“最优传输代价”：
\[
W(P, Q) = \inf_{\gamma \in \Pi(P, Q)} \mathbb{E}_{(x, y) \sim \gamma} [c(x, y)]
\]  
即使分布支撑集不同，Wasserstein 距离仍然是有限的，适用于流形上的数据建模。可解释性强，直观上等价于“移动质量”的最优成本，适用于生成模型训练。避免KL散度导致的模式崩溃，因此WGAN（Wasserstein GAN）优于传统GAN。  

## 2. 流形感知/非流形感知的深度生成模型

### 2.1 非流形感知的深度生成模型（Manifold-Unaware DGMs）

非流形感知的深度生成模型指的是未显式建模数据流形结构的生成方法。这些模型通常假设数据分布在整个高维空间中，而非流形上，因此在高维环境下可能会出现流形过拟合（Manifold Overfitting） 和 数值不稳定性 的问题。

#### 2.1.1 基于似然的方法（Likelihood-Based Approaches）

在基于似然的模型（如 VAE 和 Normalizing Flows）中，模型试图最大化对数似然：
\[
\max_{\theta} \mathbb{E}_{x \sim p^*(x)} \log p_{\theta}(x)
\]
然而，如果数据 \( x \) 仅分布在某个低维流形 \( \mathcal{M} \subset \mathbb{R}^D \) 上，而模型 \( p_{\theta}(x) \) 是定义在整个 \( \mathbb{R}^D \) 的概率密度函数，则：由于流形 \( \mathcal{M} \) 具有零测度（\( \lambda_D(\mathcal{M}) = 0 \)），因此模型需要让 \( p_{\theta}(x) \) 在 \( \mathcal{M} \) 上趋近于无穷大，才能匹配数据分布：
  \[
  p_{\theta}(x) \to \infty, \quad \forall x \in \mathcal{M}
  \]
这种现象称为 流形过拟合（Manifold Overfitting），意味着模型可能在流形上形成奇异点，而在流形外的区域密度极低或为零，从而影响泛化能力。

在高维情况下，计算对数似然需要求解雅可比行列式：
  \[
  \log p_{\theta}(x) = \log p_{z}(f^{-1}(x)) + \log \left| \det \frac{\partial f^{-1}(x)}{\partial x} \right|
  \]
  其中 \( f^{-1}(x) \) 是从数据空间到潜在空间的变换。当流形维度 \( d^* \ll D \) 时，\(\frac{\partial f^{-1}(x)}{\partial x}\) 的秩可能不足，使得雅可比行列式的计算变得不稳定，甚至出现数值溢出问题。


### 2.1.2 具体方法分析

VAE 采用变分推断，将数据 \( x \) 映射到潜在空间 \( z \)，其目标是最大化证据下界（ELBO）：
\[
\mathcal{L}_{\text{ELBO}} = \mathbb{E}_{q_{\phi}(z|x)} [\log p_{\theta}(x|z)] - KL(q_{\phi}(z|x) \| p(z))
\]
存在的问题：

- 模式崩溃（Mode Collapse）：由于 VAE 强制 \( q_{\phi}(z|x) \) 服从标准正态分布 \( \mathcal{N}(0, I) \)，如果数据的真实分布具有多个模式（modes），VAE 可能会忽略部分模式，导致生成质量下降。

- 生成模糊样本：由于最大似然目标等价于最小化均方误差（MSE），导致生成图像模糊：
     \[
     D_{\text{MSE}}(x, x') = \| x - x' \|_2^2
     \]
这意味着 VAE 可能无法准确重建数据流形上的复杂结构。

NF（归一化流） 通过构造可逆变换 \( f \) 来建模数据分布：
\[
p_{\theta}(x) = p_z(f^{-1}(x)) \left| \det J_f^{-1}(x) \right|
\]
其中 \( J_f^{-1}(x) \) 是变换的雅可比矩阵。

存在的问题：流形学习能力受限： NF 需要高维可逆映射，而数据通常分布在低维流形上，因此 NF 可能会拟合错误的概率密度。数值不稳定： 计算雅可比行列式的成本高，可能导致数值不稳定。

EBMs（基于能量的模型） 定义能量函数：
\[
p_{\theta}(x) = \frac{\exp(-E_{\theta}(x))}{Z(\theta)}
\]
其中 \( Z(\theta) \) 是配分函数。

存在的问题：计算 \( Z(\theta) \) 的成本极高，通常需要 MCMC 近似采样，优化过程难以收敛。

GANs 通过博弈过程优化生成器 \( G \) 和判别器 \( D \)：
\[
\min_G \max_D \mathbb{E}_{x \sim p_{\text{data}}} [\log D(x)] + \mathbb{E}_{z \sim p(z)} [\log (1 - D(G(z)))]
\]
存在的问题：模式崩溃（Mode Collapse）： 由于对抗训练，GANs 可能只生成部分模式的数据，而无法覆盖整个数据流形。梯度消失（Gradient Vanishing）： 当判别器过强时，生成器的梯度可能趋于零，导致训练难以继续优化。

### 2.2 流形感知的深度生成模型（Manifold-Aware DGMs）

流形感知的深度生成模型旨在明确建模数据的流形结构，从而提高生成质量和稳定性。

#### 2.2.1 通过噪声增强流形感知

去噪分数匹配（Denoising Score Matching, DSM）目标：学习数据分布的分数函数：
  \[
  s_{\theta}(x) = \nabla_x \log p_{\theta}(x)
  \]
通过向数据添加噪声，并优化去噪目标：
  \[
  \mathbb{E}_{x \sim p(x)} \mathbb{E}_{\epsilon \sim \mathcal{N}(0, \sigma^2 I)} \|\nabla_x \log p_{\theta}(x + \epsilon) - s_{\theta}(x + \epsilon)\|^2
  \]
这样，模型可以在数据流形上学习到平滑的概率梯度，从而提高生成质量。