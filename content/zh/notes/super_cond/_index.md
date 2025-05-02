---
title: 转角体系超导理论入门
date: 2025-05-02
math: true
---

论文：Coulomb interaction, phonons, and superconductivity in twisted bilayer graphene 的附录文件，包含对超导理论的简单介绍。


## 1 配对机制总览

超导配对机制的本质是电子之间通过某种方式形成库珀对，进而在低温下凝聚为相干的量子态。配对机制主要分为传统BCS机制与非常规配对机制两类。

### 1.1 传统BCS机制简介

BCS理论认为电子通过晶格振动（声子）间接相互吸引而形成库珀对。这一吸引来自电子-声子耦合：第一个电子吸引正离子形成局域势能畸变，第二个电子随后被该势能井吸引。

这种配对在动量空间表现为零总动量对，波函数为s波（无角动量）对称性。BCS能隙方程如下：

\[
\Delta = \hbar \omega_D \, \exp\left( -\frac{1}{N(0)V} \right)
\]

其中：\( \Delta \) 是超导能隙；\( \omega_D \) 是德拜频率（声子的最大频率）；\( N(0) \) 是费米能级处的态密度；\( V \) 是有效吸引势强度。

能隙的存在导致低温下电子无法以单粒子激发形式存在，从而出现零电阻。

### 1.2 非常规配对概览

当电子-电子相互作用为主导（如强库仑排斥），传统BCS机制失效，但系统仍可能通过高阶过程或特定对称性通道形成配对。

Kohn-Luttinger机制是一类典型机制，其核心思想是：

尽管一阶库仑排斥是排斥性的，但二阶交换过程可在某些角动量通道诱导出吸引成分。这种“次级诱导吸引”来源于费米面附近的粒子-空穴激发的贡献。

配对势的角动量展开：

\[
V(\theta - \theta') = \sum_{l} V_l \, e^{i l (\theta - \theta')}
\]

Kohn-Luttinger理论指出，对于高角动量通道 \( l > 0 \)，配对势 \( V_l \) 可能为负（即吸引），从而促成非s波配对（如p波、d波等）。

此外，在多轨道或多谷系统中（如TBG或TMD），存在多个电子自由度（自旋、谷、轨道），配对可以发生在不同“味道”之间，带来更复杂的配对对称性（例如s±、d+id、p+ip等）。

在这类系统中还可能出现有限动量配对（如FFLO态），其配对函数表现为：

\[
\Delta(\mathbf{r}) = \Delta_0 e^{i \mathbf{Q} \cdot \mathbf{r}}
\]

其中 \( \mathbf{Q} \neq 0 \)，即配对粒子不再具有零总动量，可能由于能带不对称或外场作用引起。

传统与非常规配对机制的区分在于：吸引来源是声子（BCS）还是电子本身的波动效应（Kohn-Luttinger、反铁磁涨落等）；以及配对函数在动量空间的结构是否简单（s波）或具备复杂节点与相位变化（如d波、p波）。


## 2 电子-声子相互作用与声子重整化

电子-声子相互作用不仅构成了传统BCS配对机制的物理基础，也会反馈影响声子的频率和传播特性。本节将系统地总结与配对机制相关的声子耦合模型、声子重整化，以及更加量化的理论如 Holstein 模型与 Migdal-Eliashberg 理论。

### 2.1 变形势耦合模型（Deformation potential）

在二维晶格材料中，声子主要通过调节晶格间距或应变作用于电子系统。对纵向声子（LA）而言，电子耦合主要通过调节化学势来实现，这种机制称为变形势（deformation potential）耦合：

\[
\delta \mu(\mathbf{r}) = D \sum_{i=x,y} \partial_i u_i(\mathbf{r})
\]

其中：
- \( D \)：变形势耦合强度，典型值如石墨烯中 \( D \approx 20\,\text{eV} \)；
- \( u_i(\mathbf{r}) \)：原子的位移场；
- \( \partial_i u_i \)：局部体积应变。

这种耦合机制会在有效哈密顿量中引入电子密度与声子位移的耦合项：

\[
H_{e\text{-}ph} = D \int d^2r\, \rho(\mathbf{r})\, \nabla \cdot \mathbf{u}(\mathbf{r})
\]

在动量空间中，对应的哈密顿量为：

\[
H_{e\text{-}ph} = \sum_{\mathbf{k},\mathbf{q}} g_{\mathbf{q}}\, c^\dagger_{\mathbf{k+q}} c_{\mathbf{k}} (b_{\mathbf{q}} + b^\dagger_{-\mathbf{q}})
\quad\text{with}\quad g_{\mathbf{q}} \propto D\, \sqrt{\frac{\hbar}{2\rho_m \omega_{\mathbf{q}}}}
\]

其中 \( \rho_m \) 为面密度，\( \omega_{\mathbf{q}} \) 是声子频率。

### 2.2 声子频率的重整化机制

由于电子的反馈效应，声子的本征频率会被修正，具体体现在声速（或等效弹性模量）被软化：

\[
\tilde{c}^2 = c^2 \left(1 - \frac{D^2 |\chi(q, \omega)|}{\lambda + 2\mu} \right)
\]

其中：
- \( \chi(q, \omega) \)：电子密度-密度响应函数；
- \( \lambda, \mu \)：拉梅系数（弹性常数）；
- \( c \)：原始声速，\( \tilde{c} \)：修正后的声速。

在 \( q \to 0 \) 极限，长程库仑相互作用使得电子响应趋于零，软化效应减弱。若系统存在金属门或高介电常数环境，屏蔽增强，声子软化更加显著。

在双层或Moiré系统中，声子模式可分为偶模（even）和奇模（odd）：
- 偶模 → 引起整体电荷密度调制，可被长程库仑势反馈，促进配对；
- 奇模 → 层间反相震动，净电荷变化为零，无法引起长程吸引，作用较弱。

### 2.3 Holstein模型：简化电子-声子配对机制

Holstein模型是研究声子介导配对最基本的格点模型。其哈密顿量为：

\[
H = -t \sum_{\langle i,j \rangle,\sigma} c^\dagger_{i\sigma} c_{j\sigma}
+ \omega_0 \sum_i b_i^\dagger b_i
+ g \sum_{i,\sigma} (b_i^\dagger + b_i) n_{i\sigma}
\]

其中：第一项是电子跃迁；第二项为本征声子；第三项为局域电子-声子耦合，耦合常数为 \( g \)。

在弱耦合极限下，可通过积分出声子自由度得到电子-电子的有效吸引势：

\[
V_{\text{eff}}(\omega) = -\frac{2g^2 \omega_0}{\omega^2 - \omega_0^2}
\]

在低频下 \( \omega \to 0 \)，近似为：

\[
V_{\text{eff}} \approx -\frac{2g^2}{\omega_0}
\]

这意味着声子交换在低能区提供一个瞬时有效吸引，从而支持BCS配对。

### 2.4 Migdal-Eliashberg理论：量子化电子-声子超导理论

在BCS理论的基础上，Migdal-Eliashberg理论进一步考虑了频率依赖、自能修正与真实声子谱的影响，是描述强耦合超导体（如Pb、Nb₃Sn）的理论标准。

该理论的核心是两个自洽的格林函数方程（Eliashberg方程），定义了电子自能 \( \Sigma(\omega) \) 和配对函数 \( \Phi(\omega) \)：

\[
\Sigma(i\omega_n) = \pi T \sum_{m} \lambda(i\omega_n - i\omega_m) \frac{i\omega_m + \Sigma(i\omega_m)}{\sqrt{\omega_m^2 + \Delta^2(i\omega_m)}}
\]

\[
\Delta(i\omega_n) = \pi T \sum_m \left[ \lambda(i\omega_n - i\omega_m) - \mu \right] \frac{\Delta(i\omega_m)}{\sqrt{\omega_m^2 + \Delta^2(i\omega_m)}}
\]

其中：\( \lambda(\omega) \)：由声子谱 \( \alpha^2F(\omega) \) 决定的谱函数；\( \mu \)：屏蔽后的库仑伪势；\( \Delta(\omega) = \Phi(\omega) / [1 - \Sigma(\omega)/\omega] \)。

相比BCS理论，Eliashberg理论：更适用于强耦合（\( \lambda \sim 1 \)）；包含有限温度与动态效应；可用于预测非s波、非零频率对称性。

当声子谱集中于某个频率（如光学声子），该理论退化为Holstein近似，成为一个一致性的理论框架。


## 3 Kohn-Luttinger配对机制（KL理论核心）

Kohn-Luttinger机制展示了一个重要但非直观的结果：即便电子之间的基本相互作用是排斥性的，在高阶（主要是二阶）交换通道中，依然可能在某些角动量通道诱导出净吸引势，从而触发超导配对。

### 3.1 高阶交换引起的吸引

考虑一般费米液体中电子-电子的排斥相互作用 \( V(\mathbf{q}) > 0 \)，其一阶贡献不支持配对。然而，二阶修正（即包含粒子-空穴激发的屏蔽过程）可导致配对通道中的有效势重整。

KL理论对有效相互作用进行角动量分解：

\[
V(\theta - \theta') = \sum_{l} V_l\, e^{i l (\theta - \theta')}
\]

KL发现对于高角动量通道（\( l > 0 \)），由于费米面附近粒子-空穴过程的非解析性（logarithmic singularity），有：

\[
V_l^{\text{eff}} = V_l^{(1)} + V_l^{(2)} < 0
\]

即使 \( V_l^{(1)} > 0 \)，但 \( V_l^{(2)} \sim -V^2 \ln (E_F/T) \) 可能主导，并形成吸引核。这类配对通常具有非平庸角动量对称性（如p波、d波、f波等），其对称性由以下因素决定：费米面几何（如是否存在嵴线、热点）；有效相互作用在不同通道中的分布；能带的多重简并或轨道/谷自由度。这类配对机制在弱耦合下不需要声子参与，完全由电子间动能驱动。

### 3.2 KL方程的一般形式（无平移对称性）

当系统无平移对称性（如Moire超晶格或无序系统），不能简单地在动量空间中构造BCS态，而需在实空间中定义配对函数 \( \Delta(r_1, r_2) \)。KL配对满足如下积分方程：

\[
\Delta(r_1, r_2) = - \int d^2r_3 d^2r_4\, V_{\text{eff}}(r_1, r_2; r_3, r_4)\, \Delta(r_3, r_4)
\]

其中 \( V_{\text{eff}} \) 为包含RPA修正、Umklapp过程（格点动量非守恒）等在内的总有效相互作用。

若考虑离散格点系统（如扭转双层石墨烯、TMD双层），需显式考虑：局域的电子态结构；Umklapp过程带来的非对角动量耦合（即 \( \mathbf{G}_1 \neq \mathbf{G}_2 \)）；局域粒子-空穴过程对 \( V_{\text{eff}} \) 的调制。

这种机制等价于说：KL机制在实空间中通过局域粒子-空穴极化诱导非局域吸引势，进而形成配对。

### 3.3 能隙方程的矩阵表示与本征值问题

在具有准周期性（如Moire晶格）或Bloch结构的体系中，可以将配对函数 \( \Delta(r_1, r_2) \) 投影到Bloch态基底：

\[
\Delta_{n_1 n_2}(\mathbf{k}) = \int d^2r_1 d^2r_2\, \Phi^_{n_1 \mathbf{k}}(r_1)\, \Phi^_{n_2 -\mathbf{k}}(r_2)\, \Delta(r_1, r_2)
\]

配对方程被写作本征值问题：

\[
\Delta_{m_1 m_2}(\mathbf{k}) = \sum_{n_1 n_2, \mathbf{q}} \Gamma_{m_1 m_2, n_1 n_2}(\mathbf{k}, \mathbf{q})\, \Delta_{n_1 n_2}(\mathbf{q})
\]

其中 \( \Gamma \) 是配对核，其矩阵元包含静态RPA修正后的有效作用：

\[
\Gamma_{m_1 m_2, n_1 n_2}(\mathbf{k}, \mathbf{q}) = - \sum_{G_1,G_2} V_{\text{scr}}^{G_1,G_2}(\mathbf{k}-\mathbf{q})\, F_{m,n}^{G_1,G_2}(\mathbf{k},\mathbf{q})
\]

最终的临界温度判据是：

\[
\lambda_{\text{max}}(T_c) = 1
\]

即配对核的最大本征值（强度）随温度下降，当其达到1时系统进入配对不稳定，形成超导态。此方法不仅可判断是否存在配对，还能得到：最有利的配对通道（主导对称性）；能隙的动量分布；临界温度与有效相互作用的关系。


## 4 配对势的构造与筛选机制

为了研究电子之间是否可能形成配对态，需要构造一个包含所有相关物理过程的有效相互作用势 \( V_{\text{eff}} \)，其关键包括屏蔽效应、Umklapp修正、电子-声子贡献等。下面逐步构建这一框架。

### 4.1 RPA下的有效相互作用

考虑一个多体系统中电子间的静态相互作用，原始库伦势为 \( V_C(q) = \frac{2\pi e^2}{\epsilon |q|} \)。在实际材料中，该相互作用会被电子的集体响应所“屏蔽”，最常见的处理方法是RPA（Random Phase Approximation）：

\[
\varepsilon(\mathbf{q}, \omega) = 1 - V_C(q)\, \chi_0(\mathbf{q}, \omega)
\]

其中：\( \chi_0(\mathbf{q}, \omega) \)：无相互作用的密度响应函数；\( \varepsilon(q, \omega) \)：介电函数。

有效相互作用由下式给出：

\[
V_{\text{scr}}(\mathbf{q}, \omega) = \frac{V_C(q)}{\varepsilon(\mathbf{q}, \omega)} = \frac{V_C(q)}{1 - V_C(q)\chi_0(\mathbf{q}, \omega)}
\]

在考虑电子-声子贡献时，可引入耦合常数 \( g = \frac{D^2}{\lambda + 2\mu} \) 修正响应函数：

\[
\chi_{\text{eff}}(\mathbf{q}, \omega) = \chi_0(\mathbf{q}, \omega) - \frac{g \chi_0^2(\mathbf{q}, \omega)}{1 + g \chi_0(\mathbf{q}, \omega)}
\]

最终的有效屏蔽势为：

\[
[V_{\text{scr}}(\mathbf{q})]^{-1} = [V_C(\mathbf{q})]^{-1} - \chi_{\text{eff}}(\mathbf{q})
\]

该表达式不仅考虑了电子静态极化，还包含声子诱导的非局域相互作用，是分析配对势的核心输入。

### 4.2 Umklapp修正与简化

在晶格周期性系统中，动量守恒只在布里渊区内成立，因此存在Umklapp过程，即：

\[
\mathbf{k}_1 + \mathbf{k}_2 = \mathbf{k}_3 + \mathbf{k}_4 + \mathbf{G}
\]

这意味着有效相互作用 \( V_{\text{eff}}(\mathbf{k}, \mathbf{k}') \) 不再仅依赖于 \( \mathbf{k} - \mathbf{k}' \)，而具有非对角的 \( \mathbf{G}_1, \mathbf{G}_2 \) 成分。构造配对势时需考虑完整的矩阵结构：

\[
V_{\text{scr}}^{G_1, G_2}(\mathbf{q}) = \varepsilon^{-1}_{G_1, G_2}(\mathbf{q})\, V_C(\mathbf{q} + \mathbf{G}_2)
\]

若采用近似简化（flat metric），可以假设态间的 form factor 可因式分解：

\[
|M_G(\mathbf{k}, \mathbf{k+q})|^2 \approx f_G^2 \cdot |\langle \mathbf{k} | \mathbf{k+q} \rangle|^2
\]

此近似允许将响应函数写为：

\[
\chi_0^{G,G'}(\mathbf{q}) = f_G f_{G'}\, \tilde{\chi}(\mathbf{q})
\]

其中 \( f_G \) 是格点 form factor，反映波函数在不同 \( \mathbf{G} \) 成分上的分布。代入RPA公式后可得到：

\[
V_{\text{scr}}^{G,G'}(\mathbf{q}) = \delta_{G,G'} V_C(\mathbf{q} + \mathbf{G}) + \frac{V_C(\mathbf{q} + \mathbf{G}) V_C(\mathbf{q} + \mathbf{G}') f_G f_{G'} \tilde{\chi}}{1 - \sum_{G''} V_C(\mathbf{q} + \mathbf{G}'') f_{G''}^2 \tilde{\chi}}
\]

从该表达式中可观察到：Umklapp修正通过非对角项 \( G \neq G' \) 引入额外耦合；如果近似为 \( f_G \ll 1 \)，则这些项对总势贡献为负（吸引）；精确处理需使用数值计算，包括非因式 form factor \( \delta \chi^{G,G'} \) 修正项。这些修正决定了配对核的具体形态及其在各对称性通道中的强度，从而对最终是否发生超导配对起决定作用。

## 5 配对对称性与轨道结构分析

在理解超导配对态时，配对函数（能隙函数）\( \Delta(\mathbf{k}) \) 的对称性具有关键意义。它不仅决定了超导态的能谱结构，还反映了相互作用的性质及配对机制的来源。特别是在复杂多轨道或多谷材料中，配对结构的内部自由度更为丰富。

### 5.1 配对函数的对称性分类

配对函数在多自由度空间（自旋、轨道、层、谷）中可按对称性划分。若记配对算符为 \( \Delta_{\alpha\beta}(\mathbf{k}) = \langle c_{\mathbf{k}\alpha} c_{-\mathbf{k}\beta} \rangle \)，则需满足反对称性：

\[
\Delta_{\alpha\beta}(\mathbf{k}) = - \Delta_{\beta\alpha}(-\mathbf{k})
\]

根据这一约束，可以分类如下：

自旋对称性（Spin structure）Spin-singlet（反对称）：\( \Delta_{\uparrow\downarrow}(\mathbf{k}) = -\Delta_{\downarrow\uparrow}(\mathbf{k}) \)，需 \( \Delta(\mathbf{k}) = \Delta(-\mathbf{k}) \)（偶函数），如s波、d波。Spin-triplet（对称）：如 \( \Delta_{\uparrow\uparrow}(\mathbf{k}) \)，需 \( \Delta(\mathbf{k}) = -\Delta(-\mathbf{k}) \)（奇函数），如p波、f波。

- 轨道结构（Orbital/valley/layer）谷间配对（inter-valley）：\( \xi = \pm \)，可发生在 \( K \) 和 \( K' \) 之间；同一层/轨道内配对（intra-layer, intra-orbital）；交叉配对（interlayer, inter-sublattice）：适用于具有多轨道结构的系统，如Moire超晶格。

动量空间对称性（Angular momentum）
角动量 \( l = 0,1,2,\dots \)，对应配对函数的极角分量：
    \[
    \Delta(\mathbf{k}) \propto e^{i l \theta_\mathbf{k}}
    \]
    - \( l = 0 \)：s波（各向同性）；
    - \( l = 1 \)：p波（奇函数）；
    - \( l = 2 \)：d波（如 \( \cos k_x - \cos k_y \)）；
    - \( l = 3 \)：f波等。

结合自旋结构，总体配对波函数必须在两粒子交换下为反对称：

\[
\text{Total wavefunction} = \text{spin} \times \text{orbital} \times \text{momentum}
\]

这意味着：若动量为偶函数（如s波），则需要spin-singlet；若动量为奇函数（如p波），则需要spin-triplet。

### 5.2 因排斥势主导而引入的有符号变化

在Kohn-Luttinger类机制或其他由排斥主导的配对机制中，吸引配对态需在动量空间中“避开”排斥作用最强的区域。这导致配对函数必然带有有符号变化（sign change），表现为节点结构。

若有效配对势 \( V(\mathbf{k}, \mathbf{k}') > 0 \) 在某区域强烈排斥，则 \( \Delta(\mathbf{k}) \) 必须在该区域变化符号，从而使配对能量：
  \[
  E_{\text{pair}} \sim \int d\mathbf{k} d\mathbf{k}'\, V(\mathbf{k}, \mathbf{k}') \Delta^(\mathbf{k}) \Delta(\mathbf{k}') < 0
  \]

结果是形成节线（nodes）或符号交替结构，例如：d波配对中，\( \Delta(\mathbf{k}) = \cos k_x - \cos k_y \) 有四个节点；p波配对中，\( \Delta(\mathbf{k}) \propto \sin k_x \) 有两个节点。

此外，排斥主导下高角动量通道的配对更具优势，因为它们天然包含多个符号变化，可以适配粒子-空穴极化中的负反馈部分。

总结而言，配对函数的对称性不仅决定了超导态的基态性质（如节点、拓扑结构），也是识别其背后机制的关键。例如：s波 → 通常由声子介导；d/p波 → 多为排斥主导；混合对称性 → 多轨道/多自由度引入（如 s±、d+id）；符号变化（sign-reversal） → 是非常规超导的重要特征。


## 6 关键参数与材料依赖性

超导临界温度 \( T_c \) 受控于系统的一系列关键物理量，这些量反映了材料的能带结构、相互作用强度与电子响应能力。建立统一模型后，理解这些参数的变化如何影响 \( T_c \)，可指导理论向不同材料系统的迁移与推广。

### 6.1 影响 \( T_c \) 的关键物理量

1.费米能级处的态密度 \( D(\varepsilon_F) \)  
   态密度决定了可参与配对的电子数量，配对强度一般与之线性相关：
   \[
   \lambda = D(\varepsilon_F) \cdot |V_{\text{eff}}|
   \]
   若系统存在平带（flat band），\( D(\varepsilon_F) \to \infty \)，理论上可实现高温超导，但需警惕局域化导致的相干性丧失。

2.有效相互作用强度 \( V_{\text{eff}} \)  
   包含声子诱导项（负）与屏蔽后的排斥项（正），通常取决于材料的介电环境：
   \[
   V_{\text{eff}} = V_{\text{ph}} + \frac{V_C}{1 - V_C \chi_0}
   \]
   增大声子耦合（如提升 \( D \)）或降低介电常数 \( \varepsilon \) 均可增强 \( |V_{\text{eff}}| \)。

3.带宽与能带压缩
   狭窄带宽增强态密度，有利于配对形成。但过窄会削弱电子的相干传播，降低全局超导相干性，需权衡。

4.介电常数与外部屏蔽
   系统所处环境（如基底、栅极）会影响电子间相互作用的屏蔽：介电常数 \( \varepsilon_{\text{env}} \) 增大 → 减弱 \( V_C \)，缓解排斥；金属门可引入镜像电荷 → 实现“动态屏蔽”或“外部压制”。

5.晶格松弛与畸变
   微小的结构重构可显著影响能带形状与DOS（如TBG中松弛诱导的带展宽），需通过第一性原理或有效模型加以建模。

综合以上因素，BCS型配对公式下的临界温度估计为：

\[
T_c \sim \omega_D \cdot \exp\left(-\frac{1}{\lambda}\right)
\quad\text{with}\quad \lambda \propto D(\varepsilon_F) \cdot |V_{\text{eff}}|
\]

在KL或非BCS机制下，虽然没有解析表达式，但依然可以将最大本征值 \( \lambda_{\text{eig}} \to 1 \) 作为临界判据，依赖于 \( D \) 与 \( V_{\text{eff}} \) 的组合。

### 6.2 材料迁移策略

为将超导理论从某一模型体系（如TBG）迁移到其他材料体系，需遵循以下策略：

1. 提取核心输入参数：
   - 能带结构：获取低能电子态的色散关系 \( \varepsilon_n(\mathbf{k}) \)，用于计算DOS与极化函数；
   - 态密度 \( D(\varepsilon_F) \)：数值积分能带后得到；
   - 响应函数 \( \chi_0(\mathbf{q}) \)：通过Kubo公式或数值格林函数计算；
   - 相互作用势 \( V_C(q) \)：依据维度与介电常数决定其形式（2D: \( \sim 1/|q| \)，3D: \( \sim 1/q^2 \)）。

2. 替换TBG特有参数
   - 替代 Moiré 超晶格下的 form factor、Umklapp 权重等，改用目标体系的本征态投影；
   - 修改 \( D \)、\( V_{\text{ph}} \)、带宽等，反映新材料特性。

3. 重构配对核与本征值问题
   - 在目标体系构建新的配对核 \( \Gamma(\mathbf{k}, \mathbf{k}') \)，求解其最大本征值判断是否存在配对不稳定；
   - 若存在多个自由度（轨道、自旋、谷），配对核应包含所有态的投影信息。

4. 评估可行性与优化路径
   - 判断系统是否具备配对形成的基本条件（高DOS、有效吸引）；
   - 探索外部调控手段（如应变、电场、界面调节）以优化 \( T_c \)。


## 7 数值方法

在分析由声子或排斥诱导的超导配对机制时，配对势的非局域性、能带结构的复杂性、以及屏蔽修正等因素通常使问题难以解析求解，因此需要依赖数值方法构建并解出配对不稳定性的本征问题。以下为配对核构造与求解的通用数值流程与数学表达。

首先，定义配对函数在动量空间的投影：

\[
\Delta_{n_1 n_2}(\mathbf{k}) = \langle c_{n_1 \mathbf{k}} c_{n_2 -\mathbf{k}} \rangle
\]

对称性要求 \( \Delta_{n_1 n_2}(\mathbf{k}) = -\Delta_{n_2 n_1}(-\mathbf{k}) \)。配对方程在静态近似下可写为离散本征问题：

\[
\lambda\, \Delta_{m_1 m_2}(\mathbf{k}) = \sum_{\mathbf{q}, n_1, n_2} \Gamma_{m_1 m_2, n_1 n_2}(\mathbf{k}, \mathbf{q})\, \Delta_{n_1 n_2}(\mathbf{q})
\]

其中 \( \Gamma(\mathbf{k}, \mathbf{q}) \) 为配对核，其典型形式为：

\[
\Gamma_{m_1 m_2, n_1 n_2}(\mathbf{k}, \mathbf{q}) = -\sum_{G_1,G_2} V_{\text{scr}}^{G_1,G_2}(\mathbf{k} - \mathbf{q}) \, M^{G_1}_{m_1 n_1}(\mathbf{k}, \mathbf{q})\, [M^{G_2}_{m_2 n_2}(\mathbf{k}, \mathbf{q})]^
\]

其中 \( M^G_{mn}(\mathbf{k}, \mathbf{q}) \) 为波函数 form factor，即布洛赫态在不同 \( G \) 成分上的投影：

\[
M^G_{mn}(\mathbf{k}, \mathbf{q}) = \sum_{\alpha} u_{m \mathbf{k}}^\alpha\, [u_{n \mathbf{q}}^\alpha]^\, e^{i G \cdot r_\alpha}
\]

数值求解的一般步骤包括：

1. 动量格点截断与采样：将布里渊区离散为 \( N_k \) 个动量点，通常取 \( N_k \sim 100-1000 \) 以保证角动量分辨率；

2. 带结构截断：限制求和带指数在费米面附近，通常选取 \( n \in \{ n_F - N_b, ..., n_F + N_b \} \)；

3. 静态响应函数的计算：
   \[
   \chi_0^{G,G'}(\mathbf{q}) = \frac{1}{N_k} \sum_{\mathbf{k},n,m} \frac{f_{n \mathbf{k}} - f_{m \mathbf{k+q}}}{\varepsilon_{n \mathbf{k}} - \varepsilon_{m \mathbf{k+q}} + i\delta} \,
   M^G_{nm}(\mathbf{k}, \mathbf{k+q})\, [M^{G'}_{nm}(\mathbf{k}, \mathbf{k+q})]^
   \]
   用于构造介电函数 \( \varepsilon_{G,G'}(\mathbf{q}) = \delta_{G,G'} - V_C^{G}(\mathbf{q}) \chi_0^{G,G'}(\mathbf{q}) \)，再反求 \( V_{\text{scr}} \)。

4. 配对核矩阵构造：将 \( \Gamma(\mathbf{k}, \mathbf{q}) \) 在 \( (\mathbf{k},n_1,n_2) \) 空间展平为一个大矩阵 \( \Gamma_{a,b} \)，维度为 \( (N_k \cdot N_b^2) \times (N_k \cdot N_b^2) \)。

5. 本征值问题求解：对该矩阵进行本征分解，寻找最大本征值 \( \lambda_{\text{max}} \) 及其对应的配对函数 \( \Delta(\mathbf{k}) \)。配对临界条件为：
   \[
   \lambda_{\text{max}}(T_c) = 1
   \]

6. 对称性分析：通过分析本征解 \( \Delta(\mathbf{k}) \) 的相位与节点结构，识别配对通道（如s波、d波等），可借助角动量展开：
   \[
   \Delta(\mathbf{k}) = \sum_l \Delta_l\, e^{i l \theta_\mathbf{k}}
   \]


## 8 总结与理论迁移路线图

Kohn-Luttinger（KL）理论与RPA修正为理解排斥诱导的非常规超导提供了一套可推广的理论框架。在不依赖强电子-声子耦合的情况下，它解释了在排斥主导条件下，通过高阶粒子-空穴过程诱导有效吸引，从而促成角动量非平庸（如d波、f波等）的配对态。本节概述该框架的适用性与迁移策略，重点详解其在 flat-band kagome 等体系中的推广。

### 8.1 KL框架与RPA耦合适用范围

KL + RPA 理论基于的基本形式如下：

配对核方程（线性化 gap equation）：

\[
\lambda\, \Delta(\mathbf{k}) = \sum_{\mathbf{q}} \Gamma(\mathbf{k}, \mathbf{q})\, \Delta(\mathbf{q})
\]

其中配对核包含屏蔽后的有效相互作用：

\[
\Gamma(\mathbf{k}, \mathbf{q}) = - V_{\text{eff}}(\mathbf{k} - \mathbf{q}) \cdot \mathcal{F}(\mathbf{k}, \mathbf{q})
\]

有效相互作用通过RPA得到：

\[
V_{\text{eff}}(q) = \frac{V_C(q)}{1 - V_C(q) \chi_0(q)}
\]

此理论适用于：中等强度排斥主导；存在结构对称性或角动量分解通道；屏蔽结构（介电函数）可明确定义；可计算或拟合得出 \( \chi_0(\mathbf{q}) \)；需要数值求解 gap 本征值问题。


### 8.2 从 TBG 推广到kagome体系

物理特征：存在高度平坦的能带（极小带宽 \( W \)）；平带处态密度 \( D(\varepsilon_F) \to \infty \)，理论上有利于配对；三角格子下的kagome结构支持高对称性角动量展开（如 \( l = 2 \) 的 d波）；脆弱的相干性与可能的局域化效应是挑战。

建模步骤：

1. 能带结构输入  
   使用 tight-binding 模型提取 kagome 能带结构，其中 flat band 一般出现在顶层或中层：
   \[
   H_{\text{TB}} = -t \sum_{\langle ij \rangle} c_i^\dagger c_j + \cdots
   \]

   得到色散关系 \( \varepsilon_n(\mathbf{k}) \)，提取 \( D(\varepsilon_F) \)：

   \[
   D(\varepsilon_F) = \frac{1}{N_k} \sum_{\mathbf{k}} \delta(\varepsilon_{\text{flat}} - \varepsilon_F)
   \quad\Rightarrow\quad D \gg 1
   \]

2. 极化函数计算 \( \chi_0(q) \)  
   即便平带电子无群速度，也存在粒子-空穴过程：

   \[
   \chi_0(\mathbf{q}) = \sum_{\mathbf{k}} \frac{f(\varepsilon_{\mathbf{k}}) - f(\varepsilon_{\mathbf{k+q}})}{\varepsilon_{\mathbf{k}} - \varepsilon_{\mathbf{k+q}} + i\eta}
   \]

   可通过小能量展开近似评估。

3. 构造屏蔽势 \( V_{\text{eff}}(q) \)  
   在二维kagome中，裸库仑势为：

   \[
   V_C(q) = \frac{2\pi e^2}{\varepsilon |q|}
   \quad\Rightarrow\quad
   V_{\text{eff}}(q) = \frac{2\pi e^2}{\varepsilon |q| - 2\pi e^2 \chi_0(q)}
   \]

   在 \( D(\varepsilon_F) \) 极大时，屏蔽增强，可导致有效吸引区间。

4. 构建配对核并角动量展开  
   构造动量空间配对核 \( \Gamma(\mathbf{k}, \mathbf{q}) \)，投影到角动量通道：

   \[
   \Gamma_l = \int d\theta\, d\theta'\, e^{-i l \theta} \Gamma(\theta - \theta') e^{i l \theta'}
   \]

   分析最具吸引性的通道（如 \( l = 2 \)，即 d波）。

5. 求解本征值问题并提取 \( T_c \)  
   形成离散化后的线性本征值问题：

   \[
   \lambda \vec{\Delta} = \Gamma \cdot \vec{\Delta}
   \]

   若 \( \lambda_{\text{max}}(T_c) = 1 \)，则对应超导临界点。

6. 局域化效应修正  
   Flat band 易受局域化影响。可定义相干重正化因子：

   \[
   Z_{\text{coh}} = \langle \psi_{\mathbf{k}} | \mathcal{P}_{\text{nonloc}} | \psi_{\mathbf{k}} \rangle
   \quad\Rightarrow\quad
   \tilde{\lambda} = Z_{\text{coh}}^2 \cdot \lambda
   \]

   若 \( Z_{\text{coh}} \ll 1 \)，超导相干性被破坏，需通过调节结构或杂质优化。

总结：Flat-band kagome 为 KL + RPA 理论的重要适用候选；高DOS可有效放大微弱吸引；高角动量通道（如d波、f波）天然具备符号变化，可适配 KL 势；需结合局域态分析与结构调控提高配对稳定性。




