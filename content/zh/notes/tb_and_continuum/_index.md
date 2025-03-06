---
title: 凝聚态计算中的紧束缚模型与连续模型
date: 2024-06-17
---

在阅读石墨烯和TMD相关文献时，发现对紧束缚模型和连续模型的理解不够清晰，特别是转角体系的连续模型是如何构造的，故在此记录。

这里是一个关于 凝聚态计算中的紧束缚模型与连续模型 的笔记大纲，涵盖基本概念、数学框架、计算方法及应用等内容：  

---

## 1 紧束缚模型

### 1.1 简介

紧束缚模型是一种描述晶体中电子运动的有效哈密顿量模型，基于电子波函数主要局域于原子轨道的假设。它通过考虑电子在不同原子间的跃迁（跃迁积分）来描述能带结构。假设：电子的波函数主要由各个格点上的局域轨道（原子轨道或Wannier轨道）组成，哈密顿量仅包含有限个近邻原子之间的跃迁项（和高中化学学习的原子结构相关，例如$d$轨道等）。核心方程为：

\[
  H = \sum_{i} \epsilon_i c_i^\dagger c_i + \sum_{\langle i,j \rangle} t_{ij} c_i^\dagger c_j + h.c.
\]
其中，\( \epsilon_i \) 是原子轨道的能量，\( t_{ij} \) 是从 \( j \) 站点跃迁到 \( i \) 站点的跃迁矩阵元，\( c_i^\dagger, c_i \) 分别是电子的产生和湮灭算符。TB的优点在于：

- 适用于局域化电子（如过渡金属化合物中的 d 轨道电子）。  
- 结构清晰，物理直观，适用于分析电子在特定晶格结构中的行为。  
- 可以灵活调整哈密顿量的参数，以适配实验数据或第一性原理计算结果（DFT计算TB参数，TB进行大规模计算）。  

### 1.2 物理基础

紧束缚模型基于以下假设：  

1. 电子主要局域在原子轨道上：每个原子上的电子波函数 \( \phi_n(\mathbf{r} - \mathbf{R}_i) \) 局域于格点 \( \mathbf{R}_i \) 附近，其中 \( n \) 代表不同的轨道（如 s、p、d 轨道）。
2. 电子波函数在不同格点上的重叠较小：即波函数 \( \phi_n(\mathbf{r} - \mathbf{R}_i) \) 之间的重叠很小，因此电子的运动主要由近邻原子之间的跃迁决定。  
3. 哈密顿量主要由最近邻、次近邻等有限范围的跃迁项决定，远程跃迁可以忽略。  

基于这些假设，电子的总波函数可以写为局域轨道的线性组合：
\[
\Psi = \sum_{i, n} c_{i, n} \phi_n(\mathbf{r} - \mathbf{R}_i)
\]
其中，\( c_{i, n} \) 是待求的系数，决定了电子在不同原子轨道上的分布。

哈密顿量由电子在原子轨道上的能量项和它们在不同格点之间的跃迁项组成。对于单轨道的情况，一般形式为：  
\[
H = \sum_{i} \epsilon_i c_i^\dagger c_i + \sum_{\langle i,j \rangle} t_{ij} c_i^\dagger c_j + h.c.
\]
其中：\( \epsilon_i \) 是格点 \( i \) 处的局域轨道能量。\( t_{ij} \) 是电子从格点 \( j \) 跃迁到 \( i \) 的跃迁积分（hopping integral）。\( c_i^\dagger \) 和 \( c_i \) 分别是电子在格点 \( i \) 处的产生和湮灭算符。\( h.c. \) 代表取厄米共轭项，确保哈密顿量为厄米矩阵。  

### 1.3 计算方法

考虑一个一维晶格，每个格点上有一个局域轨道，电子只能在最近邻原子之间跃迁。哈密顿量为：
\[
H = \sum_{i} \epsilon c_i^\dagger c_i + \sum_{i} t (c_i^\dagger c_{i+1} + c_{i+1}^\dagger c_i)
\]
在动量空间下，可以利用傅里叶变换：
\[
c_i = \frac{1}{\sqrt{N}} \sum_k e^{i k R_i} c_k
\]
将哈密顿量变换为：
\[
H = \sum_k (\epsilon + 2t \cos k a) c_k^\dagger c_k
\]
能带色散关系：
\[
E(k) = \epsilon + 2t \cos k a
\]
其中 \( a \) 是晶格常数。该模型描述了一维电子系统中的能带结构，其能量色散关系呈现简谐函数形式。

在多轨道体系中，每个原子可能包含多个电子轨道（如 \( s, p, d \) 轨道）。此时，哈密顿量需要扩展为：
\[
H = \sum_{i, n} \epsilon_n c_{i,n}^\dagger c_{i,n} + \sum_{\langle i,j \rangle, m, n} t_{mn}^{ij} c_{i,m}^\dagger c_{j,n} + h.c.
\]
其中：\( m, n \) 代表不同的轨道类型（如 \( p_x \), \( p_y \), \( d_{xy} \) 等）。\( t_{mn}^{ij} \) 代表不同轨道之间的跃迁积分。  

此外，考虑到多邻近项（beyond nearest-neighbor hopping），哈密顿量可以推广为：
\[
H = \sum_i \epsilon_i c_i^\dagger c_i + \sum_{\langle i,j \rangle} t_{ij} c_i^\dagger c_j + \sum_{\langle\langle i,j \rangle\rangle} t'_{ij} c_i^\dagger c_j + h.c.
\]
其中：\( \langle i,j \rangle \) 表示最近邻跃迁。\( \langle\langle i,j \rangle\rangle \) 表示次近邻跃迁，跃迁积分 \( t' \) 可以用来描述远程耦合效应（如自旋轨道耦合）。  

在考虑自旋轨道耦合（spin-orbit coupling, SOC）时，需要将电子的自旋自由度纳入哈密顿量。例如，在Kane-Mele模型（石墨烯中的SOC修正）中，次近邻跃迁具有自旋轨道耦合项：
\[
H_{\text{SOC}} = i \lambda_{\text{SO}} \sum_{\langle\langle i,j \rangle\rangle} \nu_{ij} c_i^\dagger s^z c_j
\]
其中\( \lambda_{\text{SO}} \) 是自旋轨道耦合强度。\( \nu_{ij} = \pm 1 \) 取决于跃迁方向。\( s^z \) 是自旋 Pauli 矩阵。此项的引入会导致拓扑非平庸能带结构，如拓扑绝缘体的能隙打开。  

此外，在某些拓扑材料（如 Weyl 半金属）中，自旋轨道耦合会导致手性跃迁：
\[
H_{\text{Weyl}} = \sum_k \psi^\dagger (v_F \mathbf{k} \cdot \mathbf{\sigma}) \psi
\]
其中 \( v_F \) 是费米速度，\( \mathbf{\sigma} \) 是 Pauli 矩阵，导致Weyl点的形成。

## 2 连续模型（Continuum Model）

### 2.1 简介

连续模型是一种低能有效理论，通常用于描述电子在动量空间的行为，特别是接近某些特殊点（如布里渊区中心或狄拉克点）时的行为。它采用近似的连续微分方程来描述能带结构，而非离散格点模型。  主要假设：电子在低能范围内的行为可以用动量展开（k·p 方法）或有效场论的方法来描述，电子态可用平面波或其他平滑函数展开（类似于泰勒展开取前面几项的近似）。典型的连续模型方程（以Dirac型电子为例）：  
\[
  H = v_F (\sigma_x k_x + \sigma_y k_y)
\]
其中，\( v_F \) 是费米速度，\( \sigma_x, \sigma_y \) 是Pauli矩阵，\( k_x, k_y \) 是动量分量。这一方程描述了石墨烯等材料中的狄拉克费米子行为。 连续模型的优点在于：

- 适用于低能有效理论，特别是描述能带中的特殊点（如Dirac点或Weyl点）。  
- 计算复杂度较低，适用于研究物理机制和推导解析结果。  
- 可以与紧束缚模型连接，通过低能展开获得。  

### 2.2 物理基础

在固体物理中，电子的行为通常由哈密顿量描述：
\[
H = H_0 + H'
\]
其中，\( H_0 \) 是系统的主导项，而 \( H' \) 是相对较小的修正项。在能带结构中，我们通常选择费米面附近的低能电子态，以线性或二次展开的方式构造有效理论，从而得到简单而可解析的描述。

k·p方法是一种在布里渊区内某个特殊点（如 \(\Gamma\) 点、\(K\) 点）附近展开哈密顿量的微扰理论，基于能带理论中的动量算符近似：
\[
H(\mathbf{k}) = H_0 + \hbar \mathbf{k} \cdot \mathbf{p} + O(k^2)
\]
其中\( H_0 \) 是布里渊区特殊点处的哈密顿量。\( \mathbf{k} \cdot \mathbf{p} \) 近似描述了电子的动量依赖。高阶项 \( O(k^2) \) 进一步修正低能行为。k·p方法可以推广至考虑自旋轨道耦合（SOC）、应变效应等情况，例如在拓扑绝缘体与二维材料的研究中，它用于构造低能有效模型并计算拓扑不变量。

### 2.3 计算方法

在某些材料（如石墨烯、拓扑绝缘体）中，电子的低能有效行为类似于相对论性粒子，可以由Dirac方程描述：
\[
H = v_F (\sigma_x k_x + \sigma_y k_y)
\]
其中：\( v_F \) 是费米速度；\( \sigma_x, \sigma_y \) 是Pauli矩阵，描述伪自旋（如子晶格自由度）；\( k_x, k_y \) 是电子的波矢分量。

TMDs（过渡金属二硫族化合物，如 MoS\(_2\)）的低能有效模型可由k·p方法构造：
\[
H_{\text{TMD}} = at(\tau k_x \sigma_x + k_y \sigma_y) + \lambda \tau \sigma_z s_z
\]
其中：\( \tau = \pm1 \) 代表 \( K/K' \) 谷；\( \lambda \) 代表自旋轨道耦合（SOC）效应。

该模型可用于研究谷电子学和自旋霍尔效应。

## 3 总结

紧束缚模型（Tight-Binding Model, TBM）与连续模型（Continuum Model, CM）是凝聚态物理中研究电子结构的两种常见理论框架。紧束缚模型基于晶格上的离散哈密顿量，而连续模型通常用于低能近似，采用连续变量进行描述。连续模型通常可以从紧束缚模型的低能展开导出。

紧束缚模型中的电子态通常由离散格点上的局域轨道（如Wannier轨道）构成，而连续模型则基于动量空间的低能展开。为了从紧束缚模型导出连续模型，我们可以在动量空间进行低能近似，通常采用布里渊区中特定点（如 \( \Gamma \)、\( K \)附近的展开。  

考虑一维紧束缚模型，其哈密顿量可写为：
\[
H = \sum_n \left[ t c_n^\dagger c_{n+1} + t^* c_{n+1}^\dagger c_n \right]
\]
在动量表象下，该哈密顿量变为：
\[
H(k) = 2t \cos(k a)
\]
其中 \( a \) 是晶格常数，\( k \) 是第一布里渊区中的动量。

在小 \( k \) 近似（\( k \to 0 \) 处展开）下，我们可以将余弦项展开：
\[
H(k) \approx 2t - t k^2 a^2
\]
去掉常数项 \( 2t \)，得到一个类自由电子的有效哈密顿量：
\[
H_{\text{eff}} = - t a^2 k^2
\]
这表明在低能区，该系统表现为一个质量为 \( m^* = \frac{\hbar^2}{2t a^2} \) 的准粒子，自然过渡到连续模型的形式。

石墨烯的紧束缚哈密顿量为：
\[
H = -t \sum_{\langle i,j \rangle} c_i^\dagger c_j
\]
其中，\( t \) 是最近邻跃迁的跃迁能，\( \langle i,j \rangle \) 表示最近邻格点。

在动量空间，该哈密顿量可写为：
\[
H(\mathbf{k}) = -t \sum_{j=1}^{3} e^{i \mathbf{k} \cdot \mathbf{\delta}_j}
\]
其中 \( \mathbf{\delta}_j \) 是最近邻矢量。

在 \( K \) 点附近（布里渊区的角落），设 \( \mathbf{k} = \mathbf{K} + \mathbf{q} \)（其中 \( \mathbf{q} \) 是小量），可以展开：
\[
H(\mathbf{q}) \approx -t \left[ \sum_j e^{i \mathbf{K} \cdot \mathbf{\delta}_j} + i \sum_j (\mathbf{q} \cdot \mathbf{\delta}_j) e^{i \mathbf{K} \cdot \mathbf{\delta}_j} \right]
\]
通过线性化近似，我们得到Dirac形式：
\[
H_{\text{eff}} = \hbar v_F (\sigma_x q_x + \sigma_y q_y)
\]
其中 \( v_F = \frac{\sqrt{3}}{2} t a / \hbar \) 是费米速度。  

这一过程表明，在低能极限下，紧束缚模型的色散关系变为线性形式，对应于连续模型的Dirac方程。

## 4 终极实战：TBG 转角石墨烯  

转角双层石墨烯（TBG）是由两层单层石墨烯以一定角度（如魔角 \( \theta \approx 1.1^\circ \)）扭转而成的异质结构，具有超晶格周期性，展现出平带和强电子关联效应。

在研究TBG之前，先回顾单层石墨烯的紧束缚模型：
\[
H_{\text{SLG}} = -t \sum_{\langle i,j \rangle} c_i^\dagger c_j + \text{h.c.}
\]
其中： \( t \approx 2.7 \) eV 是碳原子之间的最近邻跃迁能；\( c_i^\dagger \) 和 \( c_j \) 分别是碳原子上的电子产生和湮灭算符；\( \langle i,j \rangle \) 表示最近邻碳原子对。

其动量空间表达式为：
\[
H_{\text{SLG}}(\mathbf{k}) = -t \sum_{j=1}^{3} e^{i \mathbf{k} \cdot \mathbf{\delta}_j}
\]
其中 \( \mathbf{\delta}_j \) 是最近邻向量。

在TBG中，每一层石墨烯的原子排列按照旋转角 \( \theta \) 发生扭转，并形成莫尔超晶格。双层系统的哈密顿量可以写成：

\[
H_{\text{TBG}} = H_{\text{L1}} + H_{\text{L2}} + H_{\text{int}}
\]

1. 第一层（L1）的哈密顿量：
   \[
   H_{\text{L1}} = -t \sum_{\langle i,j \rangle} c_{i}^{(1)\dagger} c_{j}^{(1)} + \text{h.c.}
   \]
2. 第二层（L2）的哈密顿量（旋转角 \( \theta \) 后）：
   \[
   H_{\text{L2}} = -t \sum_{\langle i,j \rangle} c_{i}^{(2)\dagger} c_{j}^{(2)} + \text{h.c.}
   \]
   其中，第二层的原子坐标经过一个旋转矩阵 \( R(\theta) \) 变换：
   \[
   \mathbf{r}_j^{(2)} = R(\theta) \mathbf{r}_j^{(1)}
   \]
3. 层间相互作用（最近邻层间耦合）：
   \[
   H_{\text{int}} = \sum_{i \in \text{L1}, j \in \text{L2}} T(\mathbf{r}_i^{(1)} - \mathbf{r}_j^{(2)}) c_i^{(1)\dagger} c_j^{(2)} + \text{h.c.}
   \]
   其中 \( T(\mathbf{r}) \) 是跨层跃迁矩阵元，其主要贡献是：
   \[
   T(\mathbf{r}) \approx t_{\perp} e^{-|\mathbf{r} - \mathbf{r}_0|/\lambda}
   \]
   其中 \( t_{\perp} \approx 0.3 \) eV 是垂直跃迁能，\( \lambda \) 是衰减长度。

看到这里已经开始手足无措了，由于 TBG 具有较大的莫尔周期，该紧束缚哈密顿量变得计算复杂，因此通常转而采用连续模型进行简化。

在小角度极限下对紧束缚模型进行展开，构造低能连续模型。这一方法最早由 Bistritzer 和 MacDonald 提出，被称为  Bistritzer-MacDonald 模型 。

首先，在单层石墨烯的低能极限，电子可以用 Dirac 费米子  近似：
\[
H_{\text{SLG}} = \hbar v_F (\sigma_x k_x + \sigma_y k_y)
\]
其中 \( v_F = \frac{\sqrt{3}}{2} t a / \hbar \) 是费米速度。

对双层石墨烯来说，每一层都可以用 Dirac 哈密顿量表示，但由于旋转角度的存在，第二层的动量需要进行旋转：
\[
\mathbf{k}^{(2)} = R(\theta) \mathbf{k}^{(1)}
\]

因此，TBG 的单粒子哈密顿量可写为：
\[
H_{\text{TBG}} =
\begin{bmatrix}
H_{\text{L1}}(\mathbf{k}) & T(\mathbf{r}) \\
T^\dagger(\mathbf{r}) & H_{\text{L2}}(\mathbf{k})
\end{bmatrix}
\]

在小角度极限下，TBG 的有效低能哈密顿量可以写为：
\[
H_{\text{BM}} =
\begin{bmatrix}
\hbar v_F (\sigma_x k_x + \sigma_y k_y) & T(\mathbf{r}) \\
T^\dagger(\mathbf{r}) & \hbar v_F (\sigma_x k_x' + \sigma_y k_y')
\end{bmatrix}
\]
其中 \( T(\mathbf{r}) \) 是层间耦合矩阵：
\[
T(\mathbf{r}) = 
\begin{bmatrix}
u & u' e^{-i G_1 \cdot \mathbf{r}} + u' e^{-i G_2 \cdot \mathbf{r}} + u' e^{-i G_3 \cdot \mathbf{r}} \\
u' e^{i G_1 \cdot \mathbf{r}} + u' e^{i G_2 \cdot \mathbf{r}} + u' e^{i G_3 \cdot \mathbf{r}} & u
\end{bmatrix}
\]
其中：\( G_1, G_2, G_3 \) 是莫尔布里渊区的基矢量；\( u \) ，\( u' \) 是层间耦合参数。

在魔角（约 \( 1.1^\circ \)） 附近，该哈密顿量导致平带的形成，从而引发强电子关联效应。
