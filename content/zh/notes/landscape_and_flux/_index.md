---
title: 基于景观和流理论描述生物系统的涌现行为
date: 2025-03-07
math: true
---

近期在进行保研的准备，计划把各个领域都先了解一下，《Perspectives on the landscape and flux theory for describing emergent behaviors of the biological systems》，这篇文章是关于生物系统涌现行为的研究综述，计划通过阅读这篇文章来了解这个领域的前沿动态。

## 1. 基本理论 

生物系统是由大量分子和细胞相互作用形成的复杂系统，具有跨尺度的层次结构，并且呈现出涌现行为（emergent behavior）。传统的还原论方法难以完整描述这些复杂系统，因此需要采用全局性的理论框架，如景观（landscape）理论和流（flux）理论，来刻画生物系统的行为和稳定性。  

### 1.1 复杂系统的涌现行为与多尺度特性
  
在物理、生物和化学系统中，大量个体元素（如分子、细胞或生态个体）的相互作用可能导致整体行为与个体特性截然不同，即涌现行为（emergence）。这些系统通常涉及不同的时空尺度，例如：分子尺度（纳米）上的蛋白质折叠、分子识别；细胞尺度（微米）上的基因表达、信号转导；组织和生态尺度（毫米-米）上的：发育、进化、生态系统动态

许多复杂系统可以用动力系统表示：
  \[
  \frac{d \mathbf{C}}{dt} = \mathbf{F}(\mathbf{C})
  \]
其中，\(\mathbf{C} = (C_1, C_2, ..., C_n)\) 是系统状态变量（如蛋白质浓度、基因表达水平），\(\mathbf{F}(\mathbf{C})\) 表示系统的驱动力。在存在随机性（噪声）时，系统可由Langevin 方程描述：
  \[
  \frac{d \mathbf{C}}{dt} = \mathbf{F}(\mathbf{C}) + \mathbf{\eta}(t)
  \]
其中 \(\mathbf{\eta}(t)\) 是噪声项，服从高斯白噪声分布：
  \[
  \langle \eta_i(t) \eta_j(t') \rangle = 2 D_{ij} \delta(t - t')
  \]

这些数学模型可用于描述基因调控、信号传导和分子马达等生物过程。

### 1.2 景观（landscape）理论与流（flux）理论的基本概念

景观理论的核心思想是，系统的状态可以用势能函数（potential function） 或 稳态概率分布 \(P_{ss}(\mathbf{C})\) 描述：
  \[
  U(\mathbf{C}) = -\ln P_{ss}(\mathbf{C})
  \]
这个势能景观 \(U(\mathbf{C})\) 形状决定了系统的稳定性、跃迁路径及动力学行为。例如，蛋白质折叠中的漏斗景观（funnel landscape）可以用来描述折叠动力学。细胞分化的Waddington 景观用于刻画基因调控网络的不同稳态。

当系统处于非平衡稳态（nonequilibrium steady state, NESS）时，除了势能梯度外，还存在非零的稳态概率流 \(\mathbf{J}_{ss}(\mathbf{C})\)：
  \[
  \mathbf{J}_{ss}(\mathbf{C}) = \mathbf{F}(\mathbf{C}) P_{ss}(\mathbf{C}) - D \nabla P_{ss}(\mathbf{C})
  \]
如果 \(\mathbf{J}_{ss}(\mathbf{C}) \neq 0\)，表明系统存在稳定的概率循环流，破坏了详细平衡条件（detailed balance），即：
  \[
  P_{ss}(\mathbf{C}_i \to \mathbf{C}_j) \neq P_{ss}(\mathbf{C}_j \to \mathbf{C}_i)
  \]
这意味着系统在稳态下仍然具有方向性的流动，例如：细胞周期振荡（G1 → S → G2 → M → G1）；神经网络中的振荡行为；进化系统中的 Red Queen 现象。

### 1.3 生物系统中的平衡与非平衡态

在热力学平衡下，系统满足玻尔兹曼分布：
  \[
  P_{eq}(\mathbf{C}) \propto e^{-U(\mathbf{C})/k_B T}
  \]
在此情况下：系统不存在稳态流动，\(\mathbf{J}_{ss} = 0\)可以用势能景观 \(U(\mathbf{C})\) 完全描述系统行为典型例子包括蛋白质折叠、DNA 结合、受体-配体结合

在非平衡系统中，存在外部驱动（如 ATP/GTP 水解）维持系统远离平衡：
  \[
  \frac{d P(\mathbf{C}, t)}{dt} = -\nabla \cdot \mathbf{J}(\mathbf{C}, t)
  \]
在稳态情况下：
  \[
  \nabla \cdot \mathbf{J}_{ss} = 0, \quad \mathbf{J}_{ss} \neq 0
  \]
这意味着系统在稳态下仍然存在非平衡概率流，如：细胞周期（G1-S-G2-M）；细胞分化过程；神经网络中的循环活动。流来源于非平衡系统的驱动力通常来自耗散过程，如：ATP/GTP 水解（细胞生长、分裂）；电子传递链（能量代谢）；受体-配体结合与信号转导。这些过程通过非平衡流维持系统功能，表现为方向性动力学行为。

### 1.4 研究生物系统的物理与数学框架
  
要研究景观与流理论，需要结合统计物理、随机过程和非线性动力学。核心数学工具包括：

- Fokker-Planck 方程：描述系统状态的概率演化
   \[
   \frac{\partial P(\mathbf{C}, t)}{\partial t} = -\nabla \cdot \mathbf{J}(\mathbf{C}, t)
   \]
- 主方程（Master Equation）：描述离散状态系统的跃迁
   \[
   \frac{d P_i}{dt} = \sum_j (W_{ji} P_j - W_{ij} P_i)
   \]
- 路径积分方法（Path Integral）：用于计算最可能跃迁路径
   \[
   P(\mathbf{C}_f, t | \mathbf{C}_i, 0) = \int D[\mathbf{C}(t)] e^{-S[\mathbf{C}(t)]}
   \]
- 非平衡热力学：计算熵生产
   \[
   \sigma = \int \mathbf{J} \cdot \nabla (\frac{\mathbf{J}}{P})
   \]

## 2 平衡态景观

景观理论的核心思想是使用自由能势（Free Energy Potential） \(U(C)\) 描述系统状态的稳定性，低势能区域对应稳定状态，高势能区域对应过渡态或不稳定态。

### 2.1 蛋白质折叠

#### 2.1.1 Levinthal悖论与折叠路径
Levinthal 在 1969 年提出，蛋白质的天然结构是在极短时间内形成的，但如果依靠随机搜索，找到正确构象的时间远超宇宙年龄。这被称为Levinthal悖论：

设一个 100 残基蛋白，每个残基有 \(10\) 个可能构象，则构象总数为：
  \[
  10^{100}
  \]
若每个构象搜索时间为 \(10^{-13} s\)，则总搜索时间为：
  \[
  10^{100} \times 10^{-13} s \gg 10^{17} s \quad (\text{宇宙年龄} \approx 10^{17} s)
  \]
但实验表明蛋白质折叠通常在毫秒到分钟内完成。

Levinthal 悖论的解决方案：1. 能量漏斗景观（Energy Funnel Landscape）2. 分级折叠（Hierarchical Folding）3. 折叠路径上的关键中间态（Key Intermediate States）

#### 2.1.2 漏斗景观理论
蛋白质折叠并非随机搜索，而是受到自由能景观（free energy landscape）的驱动：自由能定义：
  \[
  F(C) = -k_B T \ln P_{eq}(C)
  \]
其中：\( P_{eq}(C) \) 是构象 \(C\) 的平衡概率，\( k_B \) 是玻尔兹曼常数，\( T \) 是温度。

漏斗形自由能景观 使蛋白质沿最优路径折叠：
  \[
  F(C) = F_{\text{native}} + \sum_{i} E_i e^{-\frac{\Delta E_i}{k_B T}}
  \]
其中：\( E_i \) 是构象 \(i\) 的局部自由能\( \Delta E_i \) 是其相对于本征状态的能垒。折叠路径受本征吸引（Native Basin Attraction）主导，避免随机搜索。

#### 2.1.3 能量景观的量化
使用分子动力学模拟（Molecular Dynamics, MD）计算能量景观：
  \[
  F(C) = -k_B T \ln P(C)
  \]
通过势能面分析（Potential Energy Surface, PES），识别稳定态与过渡态：
  \[
  \tau_f = \tau_0 e^{\frac{\Delta G^\dagger}{k_B T}}
  \]
  其中：\( \tau_f \) 是折叠时间，\( \Delta G^\dagger \) 是折叠过渡态自由能垒。

### 2.2 特异性与药物发现

#### 2.2.1 亲和力与特异性的概念
亲和力（Affinity）：配体与靶点结合的稳定性，通常用结合自由能衡量：
  \[
  \Delta G_{\text{binding}} = -k_B T \ln K_d
  \]
  其中：\( K_d \) 是解离常数，越小表示结合越稳定。

特异性（Specificity）：区分不同靶点的能力：
  \[
  S = \frac{\Delta E}{\Delta E_{\text{roughness}}}
  \]
  其中：\( \Delta E \) 是目标结合状态与错误结合状态的能量差，\( \Delta E_{\text{roughness}} \) 是结合景观的局部粗糙度。

#### 2.2.2 结合景观及其在药物设计中的应用
结合能景观呈漏斗形：
  \[
  F(C) = F_{\text{binding}} + \sum_{i} E_i e^{-\frac{\Delta E_i}{k_B T}}
  \]
目标是优化配体，使其具有更陡峭的漏斗，使目标靶点的结合自由能 \( \Delta G \) 最小：
  \[
  \Delta G_{\text{optimized}} = \min_{\text{ligand}} \Delta G_{\text{binding}}
  \]

### 2.3 蛋白质进化与设计

#### 2.3.1 漏斗景观的进化基础
自然选择趋向于产生具有最优折叠景观的序列，进化景观的自由能表达式：
  \[
  F(\text{sequence}) = -k_B T \ln P(\text{sequence})
  \]
进化选择优化折叠速度和稳定性：
  \[
  k_{\text{folding}} = k_0 e^{-\frac{\Delta G^\dagger}{k_B T}}
  \]

#### 2.3.2 序列空间与结构空间的关联
序列-结构映射函数：
  \[
  S = f(\text{sequence})
  \]
漏斗景观确保突变不会显著影响折叠路径。

### 2.4 生物分子结合

#### 2.4.1 结合-折叠景观
IDP（内在无序蛋白）在结合前是无序的，在结合后折叠：
  \[
  P_{\text{binding-folding}}(C) = P_{\text{binding}}(C) P_{\text{folding}}(C)
  \]
结合景观受折叠景观影响，导致结合路径多样化。

#### 2.4.2 IDP（内在无序蛋白）识别机制
IDP 采用诱导契合（Induced Fit）和构象选择（Conformational Selection）：
  \[
  P_{\text{binding}}(C) = P_{\text{unfolded}}(C) e^{-\frac{\Delta G_{\text{binding}}}{k_B T}}
  \]
  
#### 2.4.3 结合特异性的物理基础
结合自由能决定特异性：
  \[
  \Delta G_{\text{specificity}} = \Delta G_{\text{target}} - \Delta G_{\text{non-target}}
  \]
特异性取决于结合漏斗的陡峭程度。

## 3 非平衡系统的景观和流

生物系统往往处于非平衡态（nonequilibrium state），如细胞周期、基因表达调控、神经信号传递、进化动力学等过程。这些系统通常受到外部能量输入（如ATP/GTP水解）的驱动，导致概率流的形成，打破平衡态的详细平衡（detailed balance）。  非平衡景观与流理论结合势能景观（potential landscape）与稳态概率流（steady-state probability flux），刻画系统的稳定性与动力学特性。

### 3.1 非平衡动力学的基本方程

#### 3.1.1 Fokker-Planck 方程
在非平衡系统中，状态分布 \( P(C, t) \) 随时间变化，可由Fokker-Planck方程（FP方程）描述：
\[
\frac{\partial P(C,t)}{\partial t} = -\nabla \cdot J(C, t)
\]
其中，概率流 \( J(C,t) \) 由漂移项（drift）和扩散项（diffusion）组成：
\[
J(C, t) = F(C) P(C,t) - D \nabla P(C,t)
\]
\( F(C) \) 是系统的驱动力（通常来自势能梯度与非平衡流），\( D \) 是扩散系数，表示系统内随机扰动的强度。

当系统达到非平衡稳态（NESS, Nonequilibrium Steady State）时，概率分布不随时间变化：
\[
\frac{\partial P_{ss}(C)}{\partial t} = 0
\]
从而概率流满足：
\[
\nabla \cdot J_{ss}(C) = 0, \quad J_{ss}(C) \neq 0
\]
此时，系统存在一个非零的稳态概率流，导致循环流动，形成非平衡行为。

#### 3.1.2 非平衡稳态流的形成
稳态概率流 \( J_{ss}(C) \) 可进一步分解为势能梯度项与旋转流项：
\[
J_{ss}(C) = -D P_{ss}(C) \nabla U(C) + J_{\text{curl}}(C)
\]
其中：\( U(C) = -\ln P_{ss}(C) \) 是自由能景观；非平衡流 \( J_{\text{curl}}(C) \) 不能由势能景观描述，是由外部驱动维持的旋转流：
  \[
  \nabla \times J_{\text{curl}}(C) \neq 0
  \]

例如：细胞周期振荡（G1 → S → G2 → M）；神经网络振荡；生态系统的进化 Red Queen 动力学。

### 3.2 流的驱动力与生物系统的稳定性

#### 3.2.1 景观势能与流的相互作用
在非平衡系统中，动力学受两种力的影响：
\[
F(C) = -D \nabla U(C) + J_{\text{curl}}(C) / P_{ss}(C)
\]
势能梯度项 \( -D \nabla U(C) \)：吸引系统向低自由能区域移动；旋转流项 \( J_{\text{curl}}(C) \)：驱动系统绕特定轨道循环。两者的相互作用决定了系统的稳态行为。例如：细胞周期：势能景观决定稳定相位（G1/S/G2/M），旋转流驱动周期性变化，进化系统：适应度景观决定种群稳定性，旋转流驱动持续进化。

#### 3.2.2 旋转流如何驱动生物过程
非平衡流在生物系统中起着至关重要的作用：

细胞周期：
  - 景观：G1、S/G2、M 三个稳态
  - 旋转流：驱动细胞周期循环
神经网络：
  - 景观：离散的神经活动模式
  - 旋转流：维持神经振荡和决策过程
生物进化：
  - 景观：适应度高的基因型
  - 旋转流：Red Queen 动力学，维持持续适应

旋转流的计算：
\[
J_{\text{curl}}(C) = \nabla \times F(C)
\]
如果 \( J_{\text{curl}}(C) \neq 0 \)，则系统处于非平衡状态。

### 3.3 主导路径与状态间动力学

#### 3.3.1 最优路径的计算
在非平衡系统中，状态转换通常遵循最优路径，其计算可使用路径积分方法：
\[
P(C_f, t | C_i, 0) = \int D[C(t)] e^{-S[C(t)]}
\]
其中作用量 \( S[C(t)] \) 为：
\[
S[C(t)] = \int_0^t \left[ \frac{1}{4} \dot{C}^T D^{-1} \dot{C} - \frac{1}{2} J_{\text{curl}}^T D^{-1} \dot{C} + V_{\text{eff}} \right] dt
\]

最优路径满足最小作用量原理：
\[
\delta S[C(t)] = 0
\]


#### 3.3.2 速率计算与能垒分析
在非平衡系统中，跃迁速率由 Kramers 公式修正：
\[
k_{A \to B} = A e^{-\frac{\Delta G^\dagger - W}{k_B T}}
\]
其中：\( \Delta G^\dagger \) 是跃迁能垒，\( W \) 是非平衡流的贡献。如果 \( W > 0 \)，则流降低跃迁能垒，增加跃迁速率。

### 3.4 非平衡热力学

#### 3.4.1 内在能量、熵与自由能
自由能：
  \[
  F = E - TS
  \]
  其中：\( E = \int P(C) U(C) dC \) 是内在能量，\( S = -\int P(C) \ln P(C) dC \) 是熵。

非平衡自由能的演化：
  \[
  \frac{dF}{dt} = - \sigma + J_{\text{curl}} \cdot \nabla U
  \]
其中 \( \sigma \) 是熵产生率。

#### 3.4.2 自由能最小化原理
在非平衡系统中，自由能随时间减少：
\[
\frac{dF}{dt} \leq 0
\]
但不会达到最小值，而是在稳态时趋于一个非零值。

#### 3.4.3 非平衡系统的熵产生
熵产生是非平衡系统的重要特征：
\[
\sigma = \int J_{\text{curl}} \cdot D^{-1} J_{\text{curl}} dC
\]
熵产生衡量非平衡系统的能量耗散，如：细胞周期的 ATP 耗散，进化中的基因变异速率，神经系统的信息处理效率。


## 4. 非平衡景观与流的生物应用
非平衡景观与流理论在许多生物系统中具有广泛的应用，如细胞周期、细胞分化、癌症、神经网络、生物进化和跨尺度动力学。本部分通过数学模型和理论框架分析这些系统的非平衡特性，并提供相应的参考文献。

### 4.1 细胞周期

#### 4.1.1 细胞周期的调控网络
细胞周期由一系列周期性激活和抑制的基因网络调控，涉及主要的Cyclin/Cdk复合物：
- Cyclin D/Cdk4-6
- Cyclin E/Cdk2
- Cyclin A/Cdk2
- Cyclin B/Cdk1

动力学方程可用非线性微分方程表示：
\[
\frac{dC_i}{dt} = F_i(\mathbf{C})
\]
其中：\( C_i \) 表示细胞周期相关蛋白的浓度，\( F_i(\mathbf{C}) \) 由激活和抑制关系决定。

#### 4.1.2 景观与流如何驱动细胞周期振荡
细胞周期的非平衡态表现在：
1. 自由能景观：G1、S、G2、M 对应不同稳态
2. 稳态流（J_ss）：在细胞周期环路上驱动振荡

细胞周期的概率流满足：
\[
J_{\text{ss}}(C) = -D P_{\text{ss}}(C) \nabla U(C) + J_{\text{curl}}(C)
\]
势能梯度 \( -D \nabla U(C) \) 促使系统进入最低能量的稳定状态（如 G1、S、G2、M）。旋转流 \( J_{\text{curl}}(C) \) 维持细胞周期振荡

细胞周期的振荡条件：
\[
\oint J_{\text{ss}}(C) \cdot dC \neq 0
\]
意味着系统不会停留在最低自由能，而是被流驱动沿着细胞周期循环。

### 4.2 细胞分化与发育

#### 4.2.1 Waddington景观的量化
Waddington 景观用于描述细胞命运决定：低势能区域（basins） 代表不同的分化状态，高势能屏障（barriers） 阻止随机跃迁，外部信号 影响景观形状。

自由能景观：
\[
U(C) = -\ln P_{\text{ss}}(C)
\]
其中 \( C \) 表示细胞的基因表达状态。

#### 4.2.2 细胞命运决策的非平衡动力学
细胞分化遵循非平衡动力学：
\[
\frac{dC}{dt} = -D \nabla U(C) + J_{\text{curl}}(C)
\]
势能项 \( -D \nabla U(C) \) 控制细胞向稳定分化状态移动，旋转流 \( J_{\text{curl}}(C) \) 影响分化路径。

关键性状变化的预测：
\[
\frac{dS}{dt} = \frac{dS_{\text{total}}}{dt} - \frac{dS_{\text{env}}}{dt}
\]
熵的剧烈变化往往意味着分化或再生的发生。

### 4.3 癌症

#### 4.3.1 癌症形成的景观模型
癌细胞的动力学可用多稳态景观描述：
- 正常细胞态（N）
- 炎症态（G）
- 癌细胞态（C）

稳态概率：
\[
P_{\text{ss}}(C) = e^{-U(C)/k_B T}
\]

癌症形成路径：
\[
N \to G \to C
\]
癌症最可能通过炎症状态（G）形成，而非直接突变。

#### 4.3.2 关键基因和调控网络的识别
癌症调控网络可表示为非线性动力学方程：
\[
\frac{dC}{dt} = F(C) + \eta(t)
\]
- \( F(C) \) 是基因调控网络的决定性动力学
- \( \eta(t) \) 是随机噪声

### 4.4 神经网络与大脑功能

#### 4.4.1 Hopfield神经网络模型
Hopfield 模型用能量景观表示记忆存储：
\[
H = -\sum_{i<j} w_{ij} S_i S_j
\]
- \( w_{ij} \) 是神经元间的突触强度
- \( S_i \) 是神经元状态（+1 或 -1）

低能量态代表存储的记忆。

#### 4.4.2 景观与流在决策、记忆和神经振荡中的作用
大脑的决策过程可以用非平衡势能景观描述：
\[
\frac{dC}{dt} = -\nabla U(C) + J_{\text{curl}}(C)
\]
- 决策点 是能量屏障
- 旋转流 \( J_{\text{curl}} \) 影响决策方向

### 4.5 生物进化

#### 4.5.1 Wright适应度景观
适应度景观 \( W(C) \) 控制种群进化：
\[
\frac{dP(C)}{dt} = D \nabla^2 P(C) + \nabla (P(C) \nabla W(C))
\]

适应度景观的特征：
- 局部最大值 代表稳定种群
- 旋转流 解释非平衡进化

### 4.5.2 Red Queen假说的物理基础
Red Queen 假说描述物种在相互竞争下的持续进化：
\[
\oint J_{\text{ss}}(C) \cdot dC \neq 0
\]
即系统不会停留在固定适应度，而是不断变化。

### 4.6 跨尺度问题

#### 4.6.1 基因组动力学
基因组的三维结构影响基因表达：
\[
\frac{dP}{dt} = -\nabla \cdot J_{\text{ss}}(P)
\]
Hi-C 数据用于量化基因组景观。

#### 4.6.2 从分子到细胞的层次结构
非平衡景观决定跨尺度行为：
- 分子水平：蛋白质折叠
- 细胞水平：基因表达
- 组织水平：发育调控
  
#### 4.6.3 景观和流如何影响跨尺度行为
跨尺度相互作用可用分层建模：
\[
U_{\text{total}} = U_{\text{molecule}} + U_{\text{cell}} + U_{\text{tissue}}
\]
能量景观决定信息如何在不同尺度上传递。


