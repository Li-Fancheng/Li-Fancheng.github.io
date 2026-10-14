# 趋化：从随机游走到集体模式形成
## 扩散、漂移、共享化学场与聚集失稳

> 西湖大学 2026 年秋季《定量生物学导论》课程学习笔记。

> **核心问题**：如果每个细胞只是随机运动并且略微偏好化学浓度较高的方向，那么为什么许多细胞放在一起时，可能自行形成空间上的聚集结构？
>
> 细胞分泌并响应同一个化学场，就能产生间接相互作用。这里从随机运动出发，分析这种反馈如何形成。

---

## 0. 模型概览

模型可以按以下机制逐步建立：

| 步骤 | 加入的机制 | 微观规则／连续模型 | 对应现象 |
|---|---|---|---|
| 1 | 无偏随机跳跃 | 左／右各 $1/2$ | 扩散 |
| 2 | 固定的方向偏好 | $p_+\neq p_-$ | 漂移 |
| 3 | 偏好取决于化学梯度 | $p_+(x)-p_-(x)\propto \partial_x c$ | 趋化 |
| 4 | 化学物质也会扩散和衰减 | $\partial_t c=D_c\partial_x^2c-\kappa c$ | 化学环境成为动态变量 |
| 5 | 细胞分泌吸引物 | 增加 $+\alpha\rho$ | 细胞–化学场的反馈 |
| 6 | 从均匀态加入微小涨落 | 联立两场方程 | 集体聚集 |

**外场驱动的积累与自发聚集需要区分：**

- 给定外部化学场后，细胞在某些位置**积累**，并不等同于群体从均匀环境中**自发形成图案**。
- 当细胞会改变它们正在响应的化学场时，才引入了群体间的间接耦合和反馈。

本文的图均由模型计算或数值模拟生成，用于展示机制。

---

## 1. 一维无偏随机游走与扩散

### 1.1 最简单的离散规则

考虑一条间距为 $a$ 的一维格点链。每隔 $\tau$ 时间，一个粒子做一次跳跃：

$$
X_{n+1}=X_n+a\eta_n,\qquad
\eta_n=\begin{cases}
+1,&\text{概率 }1/2,\\
-1,&\text{概率 }1/2.
\end{cases}
$$

这里 $X_n$ 表示第 $n$ 步的位置，$t=n\tau$。这个模型先描述一般粒子的随机运动。

![图 1：无偏随机游走的个体轨迹](figures/01_random_walk_trajectories.png)

*图 1. 多条独立随机游走的轨迹。每条轨迹完全随机，但群体统计规律非常简单。图像由 `code/random_walk.py` 生成。*

平均位移为零，并不表示粒子没有运动。位置分布的展宽需要用均方位移描述。

因为 $\langle\eta_n\rangle=0$，而不同时间的增量独立，

$$
\langle X_n-X_0\rangle=0,
\qquad
\boxed{\langle(X_n-X_0)^2\rangle=na^2.}
$$

利用 $t=n\tau$，定义

$$
\boxed{D=\frac{a^2}{2\tau}},
\qquad
\boxed{\langle(X(t)-X_0)^2\rangle=2Dt.}
$$

**物理解释**：单个粒子的方向不可预测，但分布的宽度按 $\sqrt{t}$ 增长。这就是扩散的一个最直接的统计定义。

![图 2：随机游走的均方位移与扩散分布](figures/02_random_walk_to_diffusion.png)

*图 2. 左：模拟得到的均方位移（mean-squared displacement，MSD）与 $2Dt$ 比较；右：许多粒子在固定时间的位置直方图与高斯解比较。*

### 1.2 从主方程到扩散方程

令 $P_i^n$ 为粒子在第 $n$ 步位于格点 $i$ 的概率。粒子能在下一步到达 $i$，只能是从左边往右跳，或从右边往左跳：

$$
\boxed{P_i^{n+1}=\frac12P_{i-1}^n+\frac12P_{i+1}^n.}
\tag{1}
$$

这称为**主方程（master equation）**，描述各状态概率如何随时间变化。它对应的仍是同一个随机游走过程。

将离散概率密度写成连续变量 $\rho(x,t)$：

$$
\rho(x,t+\tau)=\frac12\rho(x-a,t)+\frac12\rho(x+a,t).
$$

在粒子分布相对于格点间距足够平滑时，展开两边：

$$
\rho+\tau\partial_t\rho+\cdots
=\rho+\frac{a^2}{2}\partial_x^2\rho+\cdots.
$$

取 $a,\tau\to0$，保持 $D=a^2/(2\tau)$ 有限，得到

$$
\boxed{\partial_t\rho=D\partial_x^2\rho.}
\tag{2}
$$

还可以把它写成局部守恒的形式：

$$
\partial_t\rho=-\partial_x J,\qquad
\boxed{J_{\rm diff}=-D\partial_x\rho.}
\tag{3}
$$

这里 $\rho$ 可理解为一个粒子的概率密度，也可理解为许多独立粒子的平均数密度（只差归一化）。**随机游走、扩散方程和 菲克定律（Fick's law），是同一个简单过程的三种描述。**

> 扩散通量指向低密度区域，是随机运动的统计结果。它与趋化漂移的竞争决定了密度起伏是衰减还是增长。

---

## 2. 有偏随机游走与漂移

### 2.1 不对称的跳跃概率

保持步长 $a$、步时 $\tau$ 不变；唯一的改变是向右的概率略大：

$$
p_+=\frac12+\varepsilon,\qquad
p_-=\frac12-\varepsilon,
\qquad |\varepsilon|\leq\frac12.
$$

此时每一步的平均位移不再是零：

$$
\langle X_{n+1}-X_n\rangle=a(p_+-p_-)=2a\varepsilon.
$$

因此

$$
\boxed{\langle X(t)-X_0\rangle=vt},
\qquad
\boxed{v=\frac{2a\varepsilon}{\tau}}.
\tag{4}
$$

![图 3：无偏与有偏随机游走](figures/03_biased_random_walk.png)

*图 3. 唯一的变化是从 $p_+=0.50$ 改为 $p_+=0.58$。群体中心产生线性漂移，分布仍保持宽度，随机性没有消失。*

偏置也会改变方差。精确结果为

$$
\operatorname{Var}(X_n-X_0)=na^2\bigl[1-(p_+-p_-)^2\bigr].
$$

所以强偏置会改变扩散部分的系数。在**弱偏置（weak bias）** 极限下（$\varepsilon\sim a$，并使 $a^2/\tau$ 有限），其主导扩散系数仍为 $D=a^2/(2\tau)$，同时 $v$ 可以保持有限。

宏观方程变为

$$
\boxed{\partial_t\rho=D\partial_x^2\rho-v\partial_x\rho}
\qquad\text{（}v\text{ 为常数）}.
\tag{5}
$$

相比无偏过程，方程多了一个**漂移项（drift）**，描述分布中心的定向移动。

### 2.2 当偏置随位置变化时

假设 $\varepsilon=\varepsilon(x)$。这时概率方程需要认真写成

$$
P_i^{n+1}=p_+(x_{i-1})P_{i-1}^n+p_-(x_{i+1})P_{i+1}^n.
\tag{6}
$$

注意概率必须在**起跳位置**计算，而不是在到达位置计算。展开式 (6)：

$$
\partial_t\rho
=\frac{a^2}{2\tau}\partial_x^2\rho
-\frac{2a}{\tau}\partial_x(\varepsilon(x)\rho).
$$

令 $v(x)=2a\varepsilon(x)/\tau$，得到

$$
\boxed{\partial_t\rho=D\partial_x^2\rho-\partial_x[v(x)\rho].}
\tag{7}
$$

这一步很值得强调：当 $v$ 随位置改变时，输运项是 **$-\partial_x(v\rho)$，而不是简单的 $-v\partial_x\rho$**。原因是粒子总数必须守恒。

> 化学浓度梯度可以使局部方向偏置 $\varepsilon(x)$ 随位置变化。

---

## 3. 趋化：由化学浓度梯度决定方向偏置

### 3.1 化学场与趋化响应

设 $c(x)$ 为外界给定的**化学吸引物浓度（chemoattractant concentration）**，即能引起细胞向其高浓度方向移动的物质的浓度。采用线性有效响应：

$$
\boxed{v(x)=\chi\partial_x c(x)},\qquad \chi\gt 0.
\tag{8}
$$

其中 $\chi$ 称为**趋化敏感系数（chemotactic sensitivity）**。$\chi\gt 0$ 表示向较高浓度漂移。与上节的格点模型对应的局部跳跃概率为

$$
p_\pm(x)=\frac12\pm\frac{\chi\tau}{2a}\partial_x c(x).
\tag{9}
$$

该线性概率规则只适用于偏置较弱、两个概率均处于 $[0,1]$ 的情形。它是对运动偏置的有效描述；细菌内部的受体与信号网络被归入响应系数。

代入式 (7)：

$$
\boxed{\partial_t\rho
=D\partial_x^2\rho
-\partial_x\bigl(\chi\rho\partial_x c\bigr).}
\tag{10}
$$

或者更直观地，用通量表达：

$$
\boxed{J=\underbrace{-D\partial_x\rho}_{\text{随机扩散}}
+\underbrace{\chi\rho\partial_xc}_{\text{定向趋化}}.}
\tag{11}
$$

扩散使密度起伏逐渐平滑，趋化漂移则把细胞输运到化学浓度更高的区域。

### 3.2 外加吸引物浓度峰

假设 $c(x)$ 在空间中心有一个峰，例如

$$
c(x)=c_{\max}\exp\left[-\frac{(x-x_0)^2}{2w^2}\right]
$$

（这里为了直观省略了有限边界的处理；代码使用周期边界下的对应光滑剖面）。若初始细胞密度为常数，演化后会在 $x_0$ 附近增加。

![图 4：固定化学场中的趋化积累](figures/04_chemotaxis_in_fixed_field.png)

*图 4. 上：预先给定、始终不随细胞变化的 $c(x)$；下：细胞密度随时间逐渐向峰附近积累。黑色虚线是零通量稳态的理论形状，**不表示**图中 $t=25$ 已达到稳态。*

在零通量稳态，

$$
J=-D\partial_x\rho+\chi\rho\partial_x c=0.
$$

因此

$$
\partial_x(\ln\rho)=\frac\chi D\partial_xc
\quad\Longrightarrow\quad
\boxed{\rho_{\rm st}(x)\propto\exp\!\left(\frac\chi D c(x)\right).}
\tag{12}
$$

比例常数由粒子总数确定。式 (12) 不仅告诉我们粒子会聚集到哪里，也说明了 $\chi/D$ 决定漂移相对于扩散的强弱。

### 3.3 外场积累与集体自组织

这种积累由外场驱动：每个粒子都独立响应预先给定的化学分布，即使没有其他粒子，也会向浓度峰移动。密度峰的位置由外界设定。

自发的集体模式形成（collective pattern formation）还需要群体反馈，也就是**细胞能够改变化学场**。

---

## 4. 化学场的动力学与细胞反馈

### 4.1 化学物质本身会发生什么？

最自然的两条规则是：

1. 化学物质在空间中扩散，扩散系数 $D_c$；
2. 化学物质在局部降解，降解率 $\kappa$。

所以首先写出

$$
\boxed{\partial_t c=D_c\partial_x^2c-\kappa c.}
\tag{13}
$$

此时化学场会演化，但仍不依赖细胞密度。若每个细胞以相同速率分泌吸引物，场的动力学就需要加入细胞源项。

如果单位长度上的细胞数是 $\rho(x,t)$，那么局部生产率应与 $\rho$ 成正比：

$$
\boxed{\partial_t c=D_c\partial_x^2 c+\alpha\rho-\kappa c.}
\tag{14}
$$

这里 $\alpha$ 表示每个细胞的信号分泌强度；数值模拟采用一致的无量纲单位。

### 4.2 细胞密度与化学场的耦合

把式 (10) 和式 (14) 放在一起：

$$
\boxed{\begin{aligned}
\partial_t\rho
&=D\partial_x^2\rho-\chi\partial_x(\rho\partial_x c),\\[3pt]
\partial_t c
&=D_c\partial_x^2 c+\alpha\rho-\kappa c.
\end{aligned}}
\tag{15}
$$

这是一类最简的 **Keller–Segel 趋化模型**。密度方程中的扩散和漂移来自随机游走，化学场方程中的各项分别对应扩散、分泌和降解。

![图 5：细胞与化学场形成的反馈回路](figures/05_field_mediated_feedback.png)

*图 5. 多一点细胞 $\to$ 多一点吸引物 $\to$ 更多细胞移入该处 $\to$ 密度增加。可编辑矢量版本见 [`figures/05_field_mediated_feedback.svg`](figures/05_field_mediated_feedback.svg)。*

共享化学场产生了间接的**细胞间相互作用（cell–cell interaction）**：细胞并不需要彼此直接接触，而是通过它们共同分泌、扩散和感知的化学物质互相影响。

### 4.3 两个简单的守恒关系

为避免边界通量造成干扰，后续数值模拟采用周期边界条件。积分第一个方程得

$$
\boxed{\frac{d}{dt}\int_0^L\rho(x,t)\,dx=0.}
\tag{16}
$$

细胞总数守恒。积分第二个方程得

$$
\boxed{\frac{d}{dt}\int_0^L c(x,t)\,dx
=\alpha\int_0^L\rho(x,t)\,dx-\kappa\int_0^L c(x,t)\,dx.}
\tag{17}
$$

总化学物质不守恒，因为细胞持续分泌且物质会降解。这两个积分关系可用于检查数值计算是否满足相应的守恒与收支关系。

---

## 5. 从均匀态到局部聚集

### 5.1 扩散与趋化正反馈的竞争

设初始细胞密度接近均匀：

$$
\rho(x,0)=\rho_0+\text{微小涨落}.
$$

在一个局部偶然有稍多细胞的地方：

- 扩散通量倾向于消除这个密度起伏；
- 细胞分泌的吸引物较多，浓度梯度把更多细胞向那里吸引；
- 后一过程可能进一步放大原来的起伏。

因此自发聚集的物理机制可以概括为

$$
\boxed{\text{扩散平滑}\quad\text{vs.}\quad
\text{趋化正反馈}.}
$$

其他参数不变时，增大 $\chi$ 会增强趋化正反馈；反馈能否胜过扩散，需要结合稳定性条件判断。

### 5.2 使用同一套初始条件，改变 $\chi$

下面的结果都由 `code/collective_1d.py` 积分式 (15) 得到，采用一维周期区间 $0\leq x\lt L$：

| 参数 | 数值 | 作用 |
|---|---:|---|
| $L$ | $20$ | 系统长度 |
| $D$ | $0.3$ | 细胞扩散 |
| $D_c$ | $1.0$ | 化学场扩散 |
| $\alpha$ | $1.0$ | 分泌强度 |
| $\kappa$ | $1.0$ | 降解率 |
| $\rho_0$ | $1.0$ | 均匀背景细胞密度 |
| $\chi$ | $0.15$ 或 $0.90$ | **唯一改变的参数** |

![图 6：弱反馈与强反馈下的一维时空密度图](figures/06_collective_1d_space_time.png)

*图 6. 横轴是空间，纵轴是时间，颜色表示 $\rho/\rho_0$。左图：微小起伏逐渐衰减；右图：起伏被反馈放大，形成局部聚集。两个色条采用各自的刻度，比较时要读数值，不能只比较颜色。*

![图 7：细胞密度与化学物质分布随时间的关系](figures/07_collective_1d_snapshots.png)

*图 7. 上：密度起伏的增长；下：在选定时刻的细胞密度与吸引物浓度峰相关。这里展示的是有限时间演化过程，未声称找到了唯一稳定的最终图案。*

这里需要区分三点：

1. **初态的涨落非常小**。真正大的密度峰是动力学放大出来的。
2. 如果初值**严格**均匀、方程也完全无噪声，均匀态仍然是精确解；必须有某个小扰动才会触发失稳。
3. **聚集的出现**与**最终是否得到稳定且间距规则的图案**是不同的问题。这个简化模型足以展示反馈失稳，但并没有保证稳定斑点 的形成。

### 5.3 从个体运动到群体反馈

$$
\underbrace{\text{单次随机跳跃}}_{\text{随机轨迹}}
\longrightarrow
\underbrace{\text{密度与通量}}_{\text{扩散／漂移}}
\longrightarrow
\underbrace{\text{细胞产生的化学场}}_{\text{动态环境}}
\longrightarrow
\underbrace{\text{群体反馈}}_{\text{聚集失稳}}.
$$

这条对应关系把个体跳跃、密度输运和共享场反馈连接起来，也说明了方程中每一项的物理来源。

---

## 6. 线性稳定性与聚集阈值

式 (15) 存在一个均匀解

$$
\rho=\rho_0,\qquad c=c_0=\frac{\alpha\rho_0}{\kappa}.
$$

对其施加微扰 $\rho=\rho_0+\delta\rho$、$c=c_0+\delta c$，保留一阶项：

$$
\begin{aligned}
\partial_t\delta\rho
&=D\partial_x^2\delta\rho-\chi\rho_0\partial_x^2\delta c,\\
\partial_t\delta c
&=D_c\partial_x^2\delta c+\alpha\delta\rho-\kappa\delta c.
\end{aligned}
$$

令 $\delta\rho,\delta c\propto e^{\sigma t+iqx}$。**如果化学场相对于细胞密度变化得够快**，可近似认为该傅里叶模态（Fourier mode）的化学场已达到局部准稳态：

$$
\delta c_q\simeq \frac{\alpha}{\kappa+D_cq^2}\delta\rho_q.
$$

代回密度方程，得到近似增长率

$$
\boxed{
\sigma(q)\simeq q^2\left[
\frac{\chi\alpha\rho_0}{\kappa+D_cq^2}-D
\right].}
\tag{18}
$$

特别是长波极限，控制竞争的无量纲参数是

$$
\boxed{\Gamma=\frac{\chi\alpha\rho_0}{D\kappa}},
\qquad \ell_c=\sqrt{\frac{D_c}{\kappa}}.
\tag{19}
$$

于是

$$
\sigma(q)\simeq Dq^2\left[\frac{\Gamma}{1+(q\ell_c)^2}-1\right].
$$

- 对无限大系统，$\Gamma\lt 1$ 时所有非零空间模式衰减；$\Gamma\gt 1$ 时存在增长的长波模式。
- 对长度为 $L$ 的周期系统，允许的非零波数是 $q_m=2\pi m/L$。只有当至少某一个允许的模式满足 $\sigma(q_m)\gt 0$ 才会形成线性失稳。对本模型，失稳阈值为 $\Gamma\gt 1+(2\pi\ell_c/L)^2$。
- 与“化学场快速响应”近似不同，**线性失稳阈值本身**也可由完整的两场线性方程精确验证，因为该 $2\times2$ 动力学矩阵的迹为负，行列式变号决定失稳。

![图 8：两组参数下的小扰动增长率](figures/08_optional_growth_rate.png)

*图 8. 用式 (18) 计算的准稳态近似增长率，黑色水平线是 $\sigma=0$；虚线竖线标示有限周期系统允许的最小非零波数。*

这一分析给出线性失稳的起始判据（onset criterion），但**不能单独决定最终图案的形状与长期稳定性**。

---

## 7. 二维模型与短时聚集

真正的空间图案通常在二维中更直观。我们只需将 $\partial_x$ 换为 $\nabla$：

$$
\begin{aligned}
\partial_t\rho &=D\nabla^2\rho-\chi\nabla\!\cdot(\rho\nabla c),\\
\partial_t c&=D_c\nabla^2c+\alpha\rho-\kappa c.
\end{aligned}
$$

![图 9：二维模拟中的短时聚集结构](figures/09_optional_2d_patterns.png)

*图 9. `code/collective_2d.py` 的二维周期模拟：从接近均匀的状态出现局部高密度区域。这是**短时动力学示意**，不是稳定 spots 的预测，更不是细菌菌落的定量拟合。最简二维模型在更长时间可能出现越来越尖的峰，数值网格和缺失的饱和机制会影响峰值。*

---

## 8. 模型假设与适用范围

模型保留随机运动、梯度响应和吸引物分泌三个机制，适用范围由以下简化决定：

- 位置依赖的跳跃概率是运动偏置的有效模型。真实大肠杆菌通过信号的时间变化调节直游与转向（run-and-tumble）；这一微观信号过程需要额外变量描述，参见 Berg 与 Brown 的实验 [2]。
- 高密度区域中的拥挤、趋化响应饱和及物质耗尽，可能影响长期演化和峰值。这些效应在当前方程中尚未包含。
- 化学波、信号兴奋性与移动带（traveling bands）涉及其他机制，需要另建模型。
- 当前模型的细胞间耦合来自共享化学场。分泌、扩散和感知共同构成反馈，使原本独立的运动变成集体动力学。

这些条件限定了模型能解释的现象，也指出了定量研究真实细菌系统时需要补充的机制。

---

## 9. 补充计算与参数比较

### 9.1 有偏随机游走的均值与方差

取 $p_+=0.6$、$p_-=0.4$，则

$\langle X_n-X_0\rangle=0.2na$，方差为 $na^2(1-0.2^2)=0.96na^2$。相对无偏过程，平均位移改变，方差系数也略有变化。

### 9.2 位置依赖漂移的守恒形式

周期边界条件下，$\frac{d}{dt}\int_0^L\rho\,dx=-[J(L)-J(0)]=0$。因此漂移项应写为 $-\partial_x(v\rho)$。当 $v(x)$ 不为常数时，只写 $-v\partial_x\rho$ 会漏掉 $-\rho\partial_xv$。

### 9.3 固定化学场中的稳态分布

零通量条件给出 $\rho_{\rm st}\propto e^{\chi c/D}$。当 $\chi\to0$ 时，趋化漂移消失，稳态趋于均匀分布。

### 9.4 切断分泌或趋化响应

在式（15）的模型中，$\alpha=0$ 或 $\chi=0$ 都会切断由细胞分泌吸引物形成的正反馈。不过，$\alpha=0$ 时若仍有外加的非均匀化学场 $c$，细胞仍可在外场驱动下积累。

### 9.5 有限周期系统的阈值

取 $L=20,D=0.3,D_c=1,\alpha=\kappa=\rho_0=1$，则 $\chi_c=(D/\alpha\rho_0)(\kappa+D_c(2\pi/L)^2)\simeq0.3296$。图 6 中的 $\chi=0.15$ 位于稳定侧，$\chi=0.90$ 位于失稳侧。

### 9.6 参数变化与数值结果

`code/collective_1d.py` 中可以分别改变 $\chi$、$D$ 或 $\kappa$，比较密度起伏的增长。增大 $\chi$ 加强趋化响应，增大 $D$ 加强扩散平滑，增大 $\kappa$ 则缩短化学信号的寿命。第 6 节的判据可用于判断这些参数变化是否会跨过失稳阈值。

---

## 10. 数值计算与图像复现

复现命令见 [`README.md`](README.md)。图像由下列脚本生成：

| 文件 | 对应内容 | 输出图片 |
|---|---|---|
| `code/random_walk.py` | 第 1–2 节 | 图 1–3 |
| `code/fixed_field.py` | 第 3 节 | 图 4 |
| `code/feedback_diagram.py` | 第 4 节 | 图 5（PNG + SVG） |
| `code/collective_1d.py` | 第 5、6 节 | 图 6–8 |
| `code/collective_2d.py` | 第 7 节 | 图 9 |
| `code/make_all_figures.py` | 一键生成所有图 | 全部 |

代码采用 NumPy + Matplotlib。连续方程的计算采用**守恒通量**和**漂移项迎风离散**，并检查细胞总数守恒；这些处理不是另加物理机制，而是为了把已经写下来的方程正确地离散到计算机上。

---

## 11. 参考文献与延伸阅读

1. **E. F. Keller & L. A. Segel (1971)**, “Model for chemotaxis,” *Journal of Theoretical Biology* **30**, 225–234. [DOI: 10.1016/0022-5193(71)90050-6](https://doi.org/10.1016/0022-5193(71)90050-6). —— 经典趋化通量模型；与本文从随机运动到通量的推导相关。
2. **H. C. Berg & D. A. Brown (1972)**, “Chemotaxis in *Escherichia coli* analysed by Three-dimensional Tracking,” *Nature* **239**, 500–504. [DOI: 10.1038/239500a0](https://doi.org/10.1038/239500a0). —— 对真实细菌改变游动方向这一现象的经典实验。
3. **E. O. Budrene & H. C. Berg (1991)**, “Complex patterns formed by motile cells of *Escherichia coli*,” *Nature* **349**, 630–633. [DOI: 10.1038/349630a0](https://doi.org/10.1038/349630a0). —— 细菌通过分泌吸引物形成 斑点与条纹 的实验现象。
4. **E. F. Keller & L. A. Segel (1971)**, “Traveling bands of chemotactic bacteria: A theoretical analysis,” *Journal of Theoretical Biology* **30**, 235–248. [DOI: 10.1016/0022-5193(71)90051-8](https://doi.org/10.1016/0022-5193(71)90051-8). —— 讨论消耗营养物所驱动的移动带，可与本文的分泌反馈模型比较。

**机制之间的联系**：从随机游走得到扩散，从方向偏置得到漂移，从化学梯度得到趋化，而**从细胞和化学场的双向反馈**才得到集体自组织的可能性。
