# 基因调控网络的理论：反馈、动力学、噪声与集体行为
## 互抑制开关、振荡、随机切换与网络响应


> 西湖大学 2026 年秋季《定量生物学导论》课程学习笔记。

---

## 概述：反馈、时间结构与随机性

基因调控网络（gene regulatory network，GRN）由基因、转录因子以及相互之间的激活和抑制关系组成。许多细胞的基因组基本相同，却能保持不同的表达状态，甚至在外部刺激消失后仍维持这种差异。理解这一现象，可以从四类问题入手：

1. **记忆**：怎样的局部调控机制允许一个系统在相同外部环境中保持两种表达状态？
2. **时间结构**：为什么某些调控网络趋向稳态，而另一些能够持续振荡？
3. **随机性**：如果基因表达来自离散、随机的反应，所谓“稳定状态”还有什么意义？
4. **集体行为**：当网络包含许多基因时，能否在不知道所有分子细节的情况下，对扰动的传播作出一般性判断？

这四个问题分别联系着非线性动力系统中的多稳态、Hopf（霍普夫）分岔、随机过程中的大偏差，以及复杂系统中的扰动传播。

**核心区分：网络拓扑决定可能的反馈结构；调控函数和时间尺度决定反馈是否真正转化为某种动力学行为；噪声进一步决定这些行为在有限时间内是否可靠。**


---

## 第一部分　从分子调控到动力学方程

### 1. GRN 中的“状态”是什么？

一个简化的描述是：基因 $i$ 的表达状态由相应的 mRNA（信使 RNA）数量 $m_i$ 和蛋白质数量或浓度 $p_i$ 表示；其中某些蛋白质是转录因子，能调控其他基因的转录速率。外部信号记为 $u$。

最简单的连续时间方程为

$$
\begin{aligned}
\frac{dm_i}{dt}&=r_i(\mathbf p,u)-\delta_{m,i}m_i,\\
\frac{dp_i}{dt}&=k_i m_i-\delta_{p,i}p_i.
\end{aligned}
\tag{1}
$$

这里 $r_i$ 是转录产生率，$k_i$ 是翻译速率，$\delta_{m,i}$ 和 $\delta_{p,i}$ 是有效降解或稀释速率。调控网络体现在函数 $r_i(\mathbf p,u)$ 中。

这个模型包含两个基本机制：

- 表达是**产生与去除的动态平衡**，而非基因简单地“开”或“关”；
- 调控关系给出产生率如何依赖其他变量，因而可能产生反馈。

一个细胞状态在此是向量 $\mathbf x=(x_1,\ldots,x_N)$ 的某个区域，而不是必须逐个基因指定为 0 或 1。后面引入的布尔描述是进一步粗粒化的结果。

### 2. Hill（希尔）函数与启动子占据率

先研究一个转录因子 $X$ 与启动子 $P$ 的结合。将有效结合过程写为

$$
P+nX\rightleftharpoons PX_n.
$$

如果结合相对于表达动力学足够快，可以近似看作平衡，定义有效解离常数 $K$，使得

$$
\frac{[PX_n]}{[P]}=\left(\frac{x}{K}\right)^n,
\qquad x=[X].
$$

假设启动子只有游离和结合这两类状态，那么结合占据率是

$$
\theta(x)
=\frac{[PX_n]}{[P]+[PX_n]}
=\frac{x^n}{K^n+x^n}.
\tag{2}
$$

当结合促进转录时，可以令

$$
r_{+}(x)=r_{\min}+(r_{\max}-r_{\min})\frac{x^n}{K^n+x^n};
\tag{3}
$$

当结合抑制转录时，可以令

$$
r_{-}(x)=r_{\min}+(r_{\max}-r_{\min})\frac{K^n}{K^n+x^n}.
\tag{4}
$$

这些函数具有阈值附近的非线性响应和高浓度时的饱和。

**重要限定。** 在更复杂的启动子模型中，$n$ 常是拟合得到的**有效 Hill 指数**，不必恰好等于真实的结合分子数；多个结合位点并不自动意味着恰好得到式 (2)。这里使用 Hill 函数，是为了保留理解反馈动力学所需的最少非线性。

**正反馈提供放大机制，但双稳态还要求适当的反馈强度与非线性。** 这些条件决定产生曲线能否与降解曲线形成多个交点。

### 3. 从 mRNA—蛋白质模型到单变量方程

考虑只含一种蛋白质的模型：

$$
\dot m=r(p)-\delta_m m,
\qquad
\dot p=k m-\delta_p p.
\tag{5}
$$

如果 $\delta_m\gg\delta_p$，且输入 $r(p)$ 的变化相对于 mRNA 弛豫也足够慢，则在蛋白质的时间尺度上，mRNA 近似满足

$$
m(t)\simeq\frac{r(p(t))}{\delta_m}.
$$

于是

$$
\dot p\simeq\frac{k}{\delta_m}r(p)-\delta_p p.
\tag{6}
$$

取浓度尺度 $K$，时间尺度 $\delta_p^{-1}$，令

$$
x=\frac pK,
\qquad \tau=\delta_p t,
$$

则方程成为

$$
\frac{dx}{d\tau}=F(x)-x,
\qquad
F(x)=\frac{k\,r(Kx)}{\delta_m\delta_p K}.
\tag{7}
$$

推广到多基因网络，可用

$$
\boxed{\frac{dx_i}{d\tau}=F_i(x_1,\ldots,x_N;u)-\gamma_i x_i.}
\tag{8}
$$

这里 $F_i$ 表示经其他分子调控后的有效表达，$\gamma_i$ 是相对于所选时间单位的降解率。

**为什么不能永远消去 mRNA？** 因为 mRNA 与蛋白质的两级动力学本身就产生响应滞后。对于振荡问题，忽略中间变量可能改变失稳条件，甚至消除振荡机制。一个合理的粗粒化，必须保留当前科学问题依赖的慢变量、增益和相位滞后。

### 4. 反馈回路的符号

在一个给定状态 $\mathbf x$ 附近，用偏导数定义有符号调控边：

$$
J_{ij}^{(\mathrm{reg})}=\frac{\partial F_i}{\partial x_j}.
$$

若 $J_{ij}^{(\mathrm{reg})}\gt 0$，则 $j\to i$ 是局部激活；若小于零，则是局部抑制。沿一个闭合回路，把所有边的符号相乘：

- 结果为正，称为**正回路**（例如一个自激活，或两个基因相互抑制）；
- 结果为负，称为**负回路**（例如三个基因依次相互抑制）。

René Thomas 提出的理论原则，在适当的连续或离散调控动力学假设下，表明正回路是出现多个稳态的**必要而非充分条件**，负回路是出现持续振荡的**必要而非充分条件** [1,2]。

这里的“必要条件”有明确模型类别和假设，不能无条件地推广到任意有时滞、受外部驱动或含特殊非局部相互作用的系统。

---

## 第二部分　正反馈：从局部非线性到细胞记忆

### 5. 自激活基因的最小模型

一个自激活基因具有表达促进自身表达的作用。考虑

$$
\boxed{\dot x=b+\alpha\frac{x^n}{1+x^n}-x,}
\tag{9}
$$

其中时间已经无量纲化；$b\ge0$ 表示基础表达，$\alpha\gt 0$ 表示最大诱导产生率，$n$ 控制响应陡峭程度。

记

$$
f(x)=b+\alpha\frac{x^n}{1+x^n}-x.
$$

固定点由 $f(x^*)=0$ 给出。对于一维连续动力系统，

$$
f'(x^*)\lt 0\ \Rightarrow\ \text{局部渐近稳定},
\qquad
f'(x^*)\gt 0\ \Rightarrow\ \text{不稳定}.
\tag{10}
$$

**物理图像：** 在固定点处，生产率等于损失率；如果浓度略微增高后损失率反过来大于生产率，系统就回到固定点，这个固定点才稳定。

#### 5.1 一个能够完全解析的例子

先设 $b=0$、$n=2$：

$$
\dot x=\alpha\frac{x^2}{1+x^2}-x.
\tag{11}
$$

固定点满足

$$
x\left(\frac{\alpha x}{1+x^2}-1\right)=0.
$$

因此，一个固定点是 $x_0=0$，另外两个满足

$$
x^2-\alpha x+1=0,
$$

从而

$$
\boxed{x_{\pm}=\frac{\alpha\pm\sqrt{\alpha^2-4}}2.}
\tag{12}
$$

只有 $\alpha\ge2$ 时，它们才是非负实固定点。当 $\alpha\gt 2$ 时，

$$
0\lt x_-\lt 1\lt x_+,
\qquad x_-x_+=1.
$$

接下来判定稳定性。对式 (11) 求导：

$$
f'(x)=\frac{2\alpha x}{(1+x^2)^2}-1.
$$

在非零固定点，由 $\alpha=(1+x^2)/x$ 可得

$$
f'(x^*)=\frac{2}{1+(x^*)^2}-1
=\frac{1-(x^*)^2}{1+(x^*)^2}.
\tag{13}
$$

于是

$$
\begin{array}{c|c|c}
\text{固定点}&f'(x^*)&\text{稳定性}\\\hline
x_0=0&-1&\text{稳定}\\
x_-\lt 1&\gt 0&\text{不稳定}\\
x_+\gt 1&\lt 0&\text{稳定}
\end{array}
$$

特别地，当 $\alpha=5/2$，

$$
\boxed{x_0=0,\quad x_-=\tfrac12,\quad x_+=2.}
\tag{14}
$$

相线结构为

```text
x = 0              x = 1/2                    x = 2
  ●  <------------   ○   ----------------->    ●
   稳定             不稳定                    稳定
```

$1/2$ 是确定性动力学中的**阈值**：初值在其两侧，将分别趋向低表达和高表达吸引子。在相同环境、相同动力学参数下，长期结果仍取决于历史初值，这就是细胞记忆的最小机制。

#### 5.2 为什么“反馈强度”存在阈值？

更一般地，仍取 $b=0$，但允许 $n\gt 1$。正的固定点需要满足

$$
\alpha=\frac{1+x^n}{x^{n-1}}=x+x^{1-n}.
\tag{15}
$$

把右边视为 $x\gt 0$ 的函数，求其极小值：

$$
\frac{d\alpha}{dx}=1-(n-1)x^{-n}=0
\quad\Longrightarrow\quad x_c^n=n-1.
$$

代回去，得到两个正固定点出现的临界值

$$
\boxed{\alpha_c=n(n-1)^{-(n-1)/n},\qquad n\gt 1.}
\tag{16}
$$

$n=2$ 时，$\alpha_c=2$，与式 (12) 一致。

当 $n=1$ 且 $b=0$ 时，$f(x)=\alpha x/(1+x)-x$ 只有 $x=0$ 和至多一个正固定点；一般并不存在上述意义的两种稳定表达状态。**反馈增益和非线性缺一不可。**

### 6. 双稳态、迟滞与记忆

式 (11) 已经有双稳态，但由于 $b=0$，$x=0$ 对所有 $\alpha$ 都是固定点。单独扫描 $\alpha$，不能据此得到完整的双向迟滞环。迟滞（hysteresis）指**控制参数正向与反向扫描时，状态在不同阈值发生突变**。

双稳态附近的一种常用正规形（normal form，即保留局部关键非线性的简化方程）为

$$
\boxed{\dot z=h+a z-z^3,\qquad a\gt 0.}
\tag{17}
$$

这里 $z$ 表示以某个参考表达状态为中心的序参量，因此可以有正负值，$h$ 则表示偏向某一表达状态的外部信号。

固定点满足

$$
h+a z-z^3=0.
$$

鞍结分岔（saddle-node bifurcation）发生在

$$
f(z)=0,\qquad f'(z)=a-3z^2=0.
$$

后一个条件给出

$$
z_{\mathrm{fold}}=\pm\sqrt{\frac a3}.
$$

再代回固定点条件：

$$
h=z^3-az=-\frac{2a}{3}z.
$$

因而两个跳变阈值为

$$
\boxed{h_{\pm}=\pm\frac{2a^{3/2}}{3\sqrt3}.}
\tag{18}
$$

当 $|h|\lt 2a^{3/2}/(3\sqrt3)$，有三个固定点：外侧两点稳定，中间点不稳定。缓慢增大 $h$ 时，系统沿低表达分支运动，直至该分支消失才跳到高表达状态；反向减小 $h$，则在另一个阈值才跳回。

**三个不同概念：**

- *双稳态*：在固定参数下共存两个稳定状态；
- *记忆*：撤去有限幅度的外部刺激后，系统仍保留被切换后的状态；
- *迟滞*：缓慢往返扫描参数，前向与后向切换点不同。

三者密切相关，但不应把它们当作同义词。尤其当扫描速度不够慢时，即使一个系统没有静态双稳态，也可能表现出由动力学滞后造成的路径依赖。

#### 6.1 保留基础表达

为了避免 $x=0$ 成为严格零表达边界，我们返回式 (9) 并令 $b\gt 0$。此时固定点条件及折叠点条件分别为

$$
b+\alpha\frac{x^n}{1+x^n}-x=0,
\qquad
\alpha\frac{n x^{n-1}}{(1+x^n)^2}-1=0.
$$

消去 $\alpha$，用折叠处的 $x$ 参数化分岔边界，得到

$$
\boxed{
\alpha_{\rm fold}(x)=\frac{(1+x^n)^2}{n x^{n-1}},
\qquad
b_{\rm fold}(x)=\frac{x}{n}(n-1-x^n).
}
\tag{19}
$$

式 (19) 明确表明：在合适的 $(b,\alpha)$ 范围内，两个鞍结分岔界定了存在三个固定点的区域。基础表达不必摧毁记忆，它通常只是改变双稳态所在的参数范围。

对于后面研究的随机切换，采用 $b\gt 0$ 尤其关键，因为在式 (11) 的纯自激活模型中，分子数严格为零时不会产生新的分子，零表达态是随机动力学的**吸收态**，无法自发逃逸。

---

## 第三部分　两个基因：互抑制开关与对称性破缺

### 7. 相互抑制为什么也是正反馈？

令 $x$ 与 $y$ 为两种互相抑制对方表达的转录因子：

$$
\boxed{
\begin{aligned}
\dot x&=\frac{\alpha}{1+y^n}-x,\\
\dot y&=\frac{\alpha}{1+x^n}-y.
\end{aligned}}
\tag{20}
$$

增加 $x$ 会抑制 $y$，而 $y$ 的减少又解除对 $x$ 的抑制。因此，两条负边构成**整体正反馈**。

这个例子在生物学上与互斥的细胞命运选择联系紧密，在动力系统中则对应多稳态和对称性破缺。Gardner、Cantor 和 Collins 于 2000 年构建的合成遗传开关提供了经典实验例子 [3]。

### 8. $n=2$ 时互抑制开关的解析解

令 $n=2$。稳态满足

$$
x(1+y^2)=\alpha,
\qquad
y(1+x^2)=\alpha.
\tag{21}
$$

两式相减：

$$
x-y+xy^2-yx^2=(x-y)(1-xy)=0.
\tag{22}
$$

因此稳态只能属于两类：

**第一类：对称稳态 $x=y=s$。**

$$
s(1+s^2)=\alpha.
\tag{23}
$$

对任意 $\alpha\gt 0$，该方程有唯一正解。

**第二类：非对称稳态 $xy=1$。**

将 $xy=1$ 代入 (21)，得到 $x+y=\alpha$，因此

$$
\boxed{
(x,y)=\left(\frac{\alpha+\sqrt{\alpha^2-4}}2,
\frac{\alpha-\sqrt{\alpha^2-4}}2\right)
}
\tag{24}
$$

及交换 $x,y$ 后的另一个解。非对称稳态只有在 $\alpha\ge2$ 时存在。

这两个分支代表“$X$ 高、$Y$ 低”和“$X$ 低、$Y$ 高”的两种表达状态。

#### 8.1 用雅可比矩阵判断稳定性

系统 (20) 在 $n=2$ 时的雅可比矩阵（Jacobian matrix）为

$$
J(x,y)=
\begin{pmatrix}
-1&-\dfrac{2\alpha y}{(1+y^2)^2}\\[5pt]
-\dfrac{2\alpha x}{(1+x^2)^2}&-1
\end{pmatrix}.
\tag{25}
$$

对称稳态 $x=y=s$ 处，定义

$$
g=\frac{2\alpha s}{(1+s^2)^2}
=\frac{2s^2}{1+s^2},
$$

第二个等号利用了 (23)。此时

$$
J_s=\begin{pmatrix}-1&-g\\-g&-1\end{pmatrix},
$$

其两个特征向量分别是 $(1,1)$ 和 $(1,-1)$，特征值为

$$
\boxed{\lambda_{\rm sym}=-1-g,
\qquad\lambda_{\rm antisym}=-1+g.}
\tag{26}
$$

由 (23)，$\alpha=2$ 对应 $s=1$；此时 $g=1$。因此：

- 当 $\alpha\lt 2$，对称稳态对两种模式都稳定；
- 当 $\alpha\gt 2$，反对称扰动增长，对称稳态变成不稳定的鞍点。

对非对称稳态，利用 $xy=1$ 可将非对角元改写成

$$
J_{12}=-\frac{2}{1+y^2},
\qquad
J_{21}=-\frac{2}{1+x^2}.
$$

由于

$$
(1+x^2)(1+y^2)=x^2+y^2+2=(x+y)^2=\alpha^2,
$$

所以两个特征值是

$$
\boxed{\lambda_{\pm}=-1\pm\frac2\alpha.}
\tag{27}
$$

当 $\alpha\gt 2$ 时，两者都为负，故**两个非对称稳态都稳定**。

我们得到完整的动力学结论：

$$
\boxed{
\alpha\lt 2:\ \text{唯一对称稳定态};
\qquad
\alpha\gt 2:\ \text{两个稳定非对称态 + 一个对称鞍点}.
}
\tag{28}
$$

$\alpha=2$ 为理想对称模型的临界值，发生对称性破缺。由于模型在交换 $x\leftrightarrow y$ 下不变，两个稳定态成对出现；在真实网络中轻微的参数不对称会偏好其中一侧，但未必立即消除双稳态。

#### 8.2 散度与周期轨道的限制

对这个最简二维模型，向量场的散度为

$$
\frac{\partial \dot x}{\partial x}
+\frac{\partial \dot y}{\partial y}=-2\lt 0.
\tag{29}
$$

由 Bendixson–Dulac 判据，在第一象限的单连通区域内不存在周期轨道。**这一个互相抑制回路很适合做开关，但在式 (20) 的动力学假设下不能成为持续振荡器。**

所以，虽然两个基因都包含负调控边，但这个回路的整体符号和可实现的动力学，与下面的三基因负反馈回路根本不同。

---

## 第四部分　负反馈：振荡何时真正出现？

### 9. 为什么负反馈通常稳定，却有时振荡？

单纯负反馈是一种偏差纠正机制：表达过高时压低生产率，表达过低时解除抑制。因此“有负反馈”并不足以推出振荡。

持续振荡要求系统纠正偏差时发生适当的超调。在线性稳定性语言中，这需要足够大的**反馈增益**与足够强的**相位滞后**，以使一对复特征值穿过虚轴。

### 10. 三基因环形抑制振荡器

三基因环形抑制振荡器（repressilator）由三种依次抑制对方表达的蛋白质构成：

$$
X_1\dashv X_2,\qquad
X_2\dashv X_3,\qquad
X_3\dashv X_1.
$$

三个负边相乘，整体构成负反馈回路。其对称模型写为

$$
\boxed{
\dot x_i=\frac{\alpha}{1+x_{i-1}^{n}}-x_i,
\qquad i=1,2,3\pmod 3.
}
\tag{30}
$$

这里 $x_0\equiv x_3$。Elowitz 和 Leibler 于 2000 年构建了这一著名的合成振荡器 [4]。

#### 10.1 对称固定点

设 $x_1=x_2=x_3=s$，则

$$
s=\frac{\alpha}{1+s^n},
\qquad\Longleftrightarrow\qquad
\alpha=s(1+s^n).
\tag{31}
$$

因为右侧对 $s\gt 0$ 单调递增，所以任意 $\alpha\gt 0$ 都有唯一的对称正固定点。

#### 10.2 在固定点附近线性化

记 $x_i=s+\eta_i$。对抑制函数求导并定义正的局部增益

$$
g=-\left.\frac{d}{dx}\frac{\alpha}{1+x^n}\right|_{x=s}
=\frac{n\alpha s^{n-1}}{(1+s^n)^2}.
$$

利用 (31) 可消去 $\alpha$：

$$
\boxed{g=\frac{n s^n}{1+s^n}.}
\tag{32}
$$

线性方程于是是

$$
\dot\eta_i=-\eta_i-g\eta_{i-1}.
\tag{33}
$$

采用环上的傅里叶模态（Fourier mode）

$$
\eta_i(\tau)=A\,e^{\lambda\tau}
\exp\!\left(\frac{2\pi\mathrm i k i}{3}\right),
\qquad k=0,1,2.
$$

因为 $\eta_{i-1}=e^{-2\pi\mathrm i k/3}\eta_i$，所以

$$
\lambda_k=-1-g\exp\!\left(-\frac{2\pi\mathrm i k}{3}\right).
\tag{34}
$$

将三个 $k$ 代入，得到

$$
\boxed{
\begin{aligned}
\lambda_0&=-1-g,\\
\lambda_{1,2}&=-1+\frac g2\ \pm\ \mathrm i\frac{\sqrt3 g}{2}.
\end{aligned}}
\tag{35}
$$

一个模态始终衰减，另两个则是一对带有相位差的复共轭模态，其实部随 $g$ 增大而上升。

#### 10.3 Hopf 阈值与 Hill 指数的作用

复特征值实部为零时，

$$
-1+\frac g2=0
\quad\Longrightarrow\quad
\boxed{g_c=2.}
\tag{36}
$$

由于 $g=ns^n/(1+s^n)\lt n$，要让 $g$ 穿过 2，必须有

$$
\boxed{n\gt 2.}
\tag{37}
$$

进一步令 $g=2$，得到

$$
\frac{n s^n}{1+s^n}=2
\quad\Longrightarrow\quad
s_H^n=\frac{2}{n-2}.
$$

对应的临界表达强度为

$$
\boxed{
\alpha_H=\frac{n}{n-2}
\left(\frac{2}{n-2}\right)^{1/n},
\qquad n\gt 2.
}
\tag{38}
$$

例如 $n=3$ 时，$\alpha_H=3\,2^{1/3}\simeq3.78$。在失稳阈值处，复特征值的虚部绝对值是

$$
\omega_H=\frac{\sqrt3 g_c}{2}=\sqrt3.
$$

因此阈值附近的无量纲振荡周期为

$$
T_H=\frac{2\pi}{\sqrt3}.
\tag{39}
$$

如果最初时间单位为蛋白质降解时间 $\delta_p^{-1}$，则物理周期为 $2\pi/(\sqrt3\,\delta_p)$。这只是**临界点附近的线性结果**，并不是在所有 $\alpha$ 下都不变的振荡周期。

#### 10.4 固定点失稳与极限环稳定性

式 (35)–(38) 证明了：当 $\alpha$ 增大并穿过阈值时，对称固定点的一对共轭特征值穿过虚轴；在满足 Hopf 分岔通常的非退化条件时，附近会出现周期解。

但**线性分析本身并不确定新生周期轨道稳定还是不稳定**。需要研究非线性正规形（第一 Lyapunov 系数），或者利用适用于特定系统的进一步分析，才能确定分岔是超临界还是亚临界。把“存在负回路”“存在 Hopf 阈值”“存在稳定、可观测的极限环”明确分开，是理论推理应有的精度。

#### 10.5 mRNA 动力学对振荡条件的影响

保留 mRNA 和蛋白质，写成

$$
\dot m_i=r(p_{i-1})-\delta_m m_i,
\qquad
\dot p_i=k m_i-\delta_p p_i.
\tag{40}
$$

在对称固定点附近，设 $q=-r'(p_*)\gt 0$。对线性扰动代入 $e^{\lambda t}$，得

$$
(\lambda+\delta_m)\,\widehat m_i=-q\,\widehat p_{i-1},
$$

$$
(\lambda+\delta_p)\,\widehat p_i=k\,\widehat m_i.
$$

消去 $\widehat m_i$：

$$
\boxed{
(\lambda+\delta_m)(\lambda+\delta_p)\widehat p_i
=-kq\,\widehat p_{i-1}.
}
\tag{41}
$$

与只保留蛋白质的一阶响应相比，现在每个调控步骤多了一个动态因子 $1/(\lambda+\delta_m)$。它给不同频率的信号带来额外的幅度衰减与相位滞后，因此临界增益与频率都会改变。

结论不是“模型越复杂越好”，而是：**当研究对象恰好是由时间滞后驱动的振荡时，某个原本看似可消去的变量可能不能消去。**

#### 10.6 奇数与偶数抑制环

把式 (33) 推广为 $N$ 个完全对称的环形抑制变量：

$$
\dot\eta_i=-\eta_i-g\eta_{i-1},\qquad i\pmod N.
$$

傅里叶模态给出

$$
\lambda_k=-1-g e^{-2\pi\mathrm ik/N}.
\tag{42}
$$

若 $N$ 为偶数，则存在 $k=N/2$，其特征值为 $\lambda=-1+g$，最先发生的可以是**实特征值失稳**：这是偶数个抑制构成整体正反馈的反映。

若 $N\ge3$ 为奇数，最容易失稳的一对模态位于相角最接近 $\pi$ 的位置，故

$$
\operatorname{Re}\lambda_{\max}=-1+g\cos\!\left(\frac\pi N\right),
$$

相应的线性振荡失稳阈值是

$$
\boxed{g_c=\sec\!\left(\frac\pi N\right).}
\tag{43}
$$

$N=3$ 时再次得到 $g_c=2$。网络尺寸与回路奇偶性，可以改变同一局部抑制作用的集体动力学后果。

---

## 第五部分　随机基因表达：吸引子为什么不是永久状态？

### 11. 最简单的随机反应：恒定生成与独立降解

以蛋白质或 mRNA 分子数 $N$ 为变量，考虑两个反应：

$$
\varnothing\xrightarrow{\beta}X,
\qquad
X\xrightarrow{\gamma}\varnothing.
\tag{44}
$$

第一种反应以恒定速率 $\beta$ 增加一个分子；第二种反应以总速率 $\gamma n$，从分子数 $n$ 减少到 $n-1$。令 $P_n(t)=\Pr(N(t)=n)$。

描述分子数概率演化的主方程（master equation）为

$$
\boxed{
\frac{dP_n}{dt}
=\beta P_{n-1}+\gamma(n+1)P_{n+1}
-(\beta+\gamma n)P_n,
}
\tag{45}
$$

并约定 $P_{-1}=0$。

其中前两项代表流入状态 $n$ 的概率流，后两项是离开状态 $n$ 的概率流。

#### 11.1 详细平衡与稳态分布

对于这个一维出生–死亡过程（birth–death process），平稳时相邻状态之间的净概率流为零，即满足详细平衡（detailed balance）：

$$
\beta P_n^*=\gamma(n+1)P_{n+1}^*.
$$

于是

$$
P_{n+1}^*=\frac{\beta}{\gamma(n+1)}P_n^*.
$$

反复迭代：

$$
P_n^*=\frac1{n!}\left(\frac\beta\gamma\right)^n P_0^*.
$$

利用归一化 $\sum_{n=0}^{\infty}P_n^*=1$，得到

$$
P_0^*=e^{-\beta/\gamma}.
$$

因此

$$
\boxed{
P_n^*=e^{-\mu}\frac{\mu^n}{n!},
\qquad
\mu=\frac\beta\gamma.
}
\tag{46}
$$

这是泊松分布（Poisson distribution），因此

$$
\boxed{\langle N\rangle=\mu,
\qquad\operatorname{Var}(N)=\mu,
\qquad\frac{\sqrt{\operatorname{Var}(N)}}{\langle N\rangle}=\mu^{-1/2}.}
\tag{47}
$$

即使完全没有细胞之间的参数差异，单个细胞的表达也会有随机涨落。当平均分子数很小时，相对涨落尤其显著。

#### 11.2 从矩方程再核对一次

利用反应带来的跃迁，直接得到

$$
\frac{d\langle N\rangle}{dt}=\beta-\gamma\langle N\rangle,
$$

$$
\frac{d\operatorname{Var}(N)}{dt}
=\beta+\gamma\langle N\rangle
-2\gamma\operatorname{Var}(N).
\tag{48}
$$

稳态时第一式给 $\mu=\beta/\gamma$，第二式给 $\operatorname{Var}(N)=\mu$，与式 (46) 相符。

**生物学含义：** 确定性方程描述平均变化，并不意味着单细胞沿着完全确定的轨道运动。Elowitz 等的单细胞实验是讨论基因表达噪声的经典起点 [5]。

### 12. 启动子切换与超泊松噪声

恒定生成率往往过于理想。考虑启动子在 OFF 与 ON 两种状态之间随机切换：

$$
\mathrm{OFF}\xrightleftharpoons[k_{\rm off}]{k_{\rm on}}\mathrm{ON}.
$$

只有在 ON 状态，以速率 $r$ 产生 mRNA；每个 mRNA 以速率 $\gamma$ 降解。用 $S(t)\in\{0,1\}$ 表示启动子是否处于 ON 状态。令

$$
\kappa=k_{\rm on}+k_{\rm off},
\qquad
\pi=\langle S\rangle_*=\frac{k_{\rm on}}\kappa.
$$

平均 mRNA 数 $M$ 满足

$$
\frac{d\langle M\rangle}{dt}=r\langle S\rangle-\gamma\langle M\rangle,
$$

所以

$$
\boxed{\langle M\rangle_*=\frac{r\pi}{\gamma}.}
\tag{49}
$$

为了计算方差，需要混合矩 $\langle SM\rangle$。由生成、降解和启动子切换分别计算其变化率：

$$
\frac{d\langle SM\rangle}{dt}
=r\langle S\rangle+k_{\rm on}\langle M\rangle
-(\gamma+\kappa)\langle SM\rangle.
\tag{50}
$$

在稳态代入 $\langle S\rangle=\pi$ 和式 (49)，可得

$$
\operatorname{Cov}(S,M)
=\langle SM\rangle_*-\pi\langle M\rangle_*
=\frac{r\pi(1-\pi)}{\gamma+\kappa}.
\tag{51}
$$

再利用二阶矩方程，

$$
\frac{d\operatorname{Var}(M)}{dt}
=r\langle S\rangle+\gamma\langle M\rangle
+2r\operatorname{Cov}(S,M)-2\gamma\operatorname{Var}(M),
$$

在稳态得到

$$
\boxed{
\operatorname{Var}(M)
=\langle M\rangle
+\frac{r^2\pi(1-\pi)}{\gamma(\gamma+\kappa)}.
}
\tag{52}
$$

方差与均值之比称为 Fano 因子（Fano factor），这里为

$$
\boxed{
\frac{\operatorname{Var}(M)}{\langle M\rangle}
=1+\frac{r(1-\pi)}{\gamma+\kappa}\gt 1
}
\tag{53}
$$

（当 $0\lt \pi\lt 1$ 且 $r\gt 0$ 时严格大于 1）。

mRNA 的涨落既来自离散生成、降解事件产生的散粒噪声（shot noise），也来自对产生率共同起调制作用的启动子状态。启动子切换越慢，表达速率越容易长时间维持高或低状态，方差就越大。这是时间相关性放大表达噪声的一个可解析例子。

---

## 第六部分　噪声驱动的状态切换：从势景观到大偏差

### 13. 为什么确定性稳定态会变成随机亚稳态？

回到双稳态模型。确定性动力学

$$
\dot x=f(x)
$$

拥有低表达稳定点 $x_L$、不稳定阈值 $x_U$ 与高表达稳定点 $x_H$：

$$
x_L\lt x_U\lt x_H.
$$

有了噪声以后，即使系统处于 $x_L$ 附近，也可能偶然越过 $x_U$，随后确定性漂移便将其带向 $x_H$。这不是吸引子“消失”，而是吸引子从严格永久的终态变成**具有有限寿命的亚稳态（metastable state）**。

#### 13.1 朗之万描述与 Fokker–Planck 方程

先采用加性白噪声的朗之万（Langevin）描述：

$$
\boxed{dx=f(x)\,dt+\sqrt{2D}\,dW_t,}
\tag{54}
$$

其中 $W_t$ 为维纳过程（Wiener process，即标准布朗运动），$D$ 为扩散强度。对应的 Fokker–Planck 方程（概率密度的演化方程）为

$$
\frac{\partial\rho}{\partial t}
=-\frac{\partial}{\partial x}[f(x)\rho]
+D\frac{\partial^2\rho}{\partial x^2}.
\tag{55}
$$

在一维中可定义

$$
U(x)=-\int^x f(y)\,dy,
\qquad f(x)=-U'(x).
\tag{56}
$$

稳态概率流为

$$
j=f\rho-D\rho'.
$$

在零概率流边界条件下令 $j=0$，得到

$$
\frac{\rho'}\rho=\frac fD=-\frac{U'}D
\quad\Longrightarrow\quad
\boxed{\rho_*(x)\propto e^{-U(x)/D}.}
\tag{57}
$$

稳定固定点对应 $U$ 的局部极小值；两稳态之间的不稳定固定点对应势垒顶部。对于式 (9) 的 $n=2$ 情形，

$$
f(x)=b+\alpha\frac{x^2}{1+x^2}-x,
$$

直接积分得到

$$
\boxed{
U(x)=\frac{x^2}{2}-bx-\alpha\bigl(x-\arctan x\bigr)+\text{constant}.
}
\tag{58}
$$

例如，选择 $b=0.02$、$\alpha=2.5$，就有三个正固定点，约为 $0.0211$、$0.4659$ 和 $2.0329$。外侧两个是稳定态；中间的点形成势垒。

#### 13.2 Kramers 逃逸与记忆寿命

设 $x_a$ 是某一势阱底部，$x_b$ 是分隔两个盆地的势垒顶部，

$$
\Delta U=U(x_b)-U(x_a)\gt 0.
$$

在一维、加性弱噪声、充分光滑的非退化势阱及适当边界条件下，逃逸时间的 Kramers 近似为

$$
\boxed{
\tau_{a\to\mathrm{other}}
\simeq\frac{2\pi}{\sqrt{U''(x_a)|U''(x_b)|}}
\exp\!\left(\frac{\Delta U}{D}\right).
}
\tag{59}
$$

指数因子比前因子更值得关注：噪声强度的微小变化，或反馈使势垒略微提高，都能造成状态寿命的巨大差异。

因此，**细胞记忆不要求状态永久不可逃逸；它只要求平均切换时间远大于相关的生物学观察时间。**

注意 $b=0$ 的自激活模型不能直接拿来讨论从零分子数态向高表达态的真实随机逃逸，因为它在 $N=0$ 处没有产生事件。保持 $b\gt 0$ 可避免这个问题。

### 14. 从化学主方程到稀有事件势垒

化学反应的产生与损失率通常随状态变化，噪声因而也是**状态依赖**的。相应的扩散近似和逃逸势垒可以从反应规则推导。

设系统大小为 $\Omega$，分子数 $N\simeq\Omega x$。在状态 $N$，产生一个分子的速率为 $\Omega B(x)$，损失一个分子的速率为 $\Omega R(x)$。写作

$$
N\xrightarrow{\ \Omega B(N/\Omega)\ }N+1,
\qquad
N\xrightarrow{\ \Omega R(N/\Omega)\ }N-1.
$$

这里 $B$ 和 $R$ 都是非负的宏观速率函数。平均场漂移是

$$
\boxed{\dot x=B(x)-R(x).}
\tag{60}
$$

#### 14.1 系统大小展开给出扩散近似

一步跃迁使 $x$ 改变 $\pm1/\Omega$。对主方程进行 Kramers–Moyal 展开至二阶，得到

$$
\frac{\partial\rho}{\partial t}
=-\partial_x\!\left[(B-R)\rho\right]
+\frac{1}{2\Omega}\partial_x^2\!\left[(B+R)\rho\right]
+O(\Omega^{-2}).
\tag{61}
$$

对应的 Itô 随机微分方程为

$$
\boxed{
dx=[B(x)-R(x)]\,dt
+\sqrt{\frac{B(x)+R(x)}{\Omega}}\,dW_t.
}
\tag{62}
$$

它说明两件事：

1. 大系统的典型表达涨落相对变小，噪声幅度是 $O(\Omega^{-1/2})$；
2. 随机性由产生与损失事件共同贡献，因此噪声通常依赖 $x$，并不等于随意加入一个常数扩散项。

#### 14.2 罕见切换与 WKB 近似

扩散近似通常能够正确描述稳定点附近的典型小涨落，但跨过吸引盆地边界是指数罕见事件，所需的“尾部”概率有时不能靠只保留二阶展开精确描述。

WKB 近似（Wentzel–Kramers–Brillouin approximation）在这里用于求解小噪声下的主导指数行为。设大 $\Omega$ 下的稳态概率具有以下形式：

$$
P_N\asymp\exp[-\Omega S(x)],
\qquad x=N/\Omega.
\tag{63}
$$

这里的 $S(x)$ 是大偏差作用量，也对应稳态分布的率函数（rate function）。化学主方程在平稳状态下为

$$
0=\Omega B(x-\Omega^{-1})P_{N-1}
+\Omega R(x+\Omega^{-1})P_{N+1}
-\Omega[B(x)+R(x)]P_N.
$$

由 (63)，在主阶有

$$
\frac{P_{N-1}}{P_N}\simeq e^{S'(x)},
\qquad
\frac{P_{N+1}}{P_N}\simeq e^{-S'(x)}.
$$

忽略相对 $O(\Omega^{-1})$ 的变化后，得到零能量 Hamilton–Jacobi 方程

$$
0=B(x)(e^q-1)+R(x)(e^{-q}-1),
\qquad q=S'(x).
\tag{64}
$$

令指数变量代换后，这成为一个二次代数方程：令 $z=e^q$，则

$$
Bz^2-(B+R)z+R=0
\quad\Longrightarrow\quad
(z-1)(Bz-R)=0.
$$

一个解是 $q=0$，对应确定性弛豫分支；另一个为

$$
\boxed{S'(x)=\ln\frac{R(x)}{B(x)}.}
\tag{65}
$$

这个非零动量分支描述从稳定态向阈值点的罕见“逆流”涨落轨道。若 $x_a$ 是稳定点、$x_b$ 是通往另一吸引盆地的中间不稳定点，则主导的逃逸势垒是

$$
\boxed{
\Delta S=\int_{x_a}^{x_b}\ln\!\frac{R(x)}{B(x)}\,dx,
\qquad
\tau_{a\to b}\asymp e^{\Omega\Delta S}.
}
\tag{66}
$$

这里 $\asymp$ 表示只给出系统大小 $\Omega$ 下的主导指数阶，未给出前因子。对于一维出生–死亡过程，这也可以从相邻态平稳概率的比值直接得到。

以双稳态的自激活模型为例，可以选

$$
B(x)=b+\alpha\frac{x^2}{1+x^2},
\qquad R(x)=x,
\tag{67}
$$

其中 $b\gt 0$。式 (66) 就是该具体化学跳跃模型的逃逸指数，而不能不加区分地替换为“确定性势函数 $U$ 的高度 / 某个常数噪声强度”。

**最重要的物理区别：** 一维 Langevin 常噪声中的 $U$、一般出生–死亡模型的大偏差作用量 $S$，以及多维非平衡系统中常说的“细胞命运景观”，是有关联但不相同的概念。

### 15. 为什么不能把所有 GRN 都画成一个普通能量景观？

在多基因系统中，漂移是向量场 $\mathbf f(\mathbf x)$。一般没有标量函数 $U(\mathbf x)$ 使

$$
\mathbf f=-\nabla U.
$$

存在稳态概率环流时，平稳分布不能仅用详细平衡来解释。多维随机动力学可以定义大偏差准势（quasipotential），用来描述罕见跃迁的代价，但真实轨迹通常不沿准势的负梯度运动。

Waddington 景观（Waddington landscape，即用地形类比细胞状态与分化路径）有助于直观理解状态选择，却不能直接代替具体 GRN 的动力学方程。

---

## 第七部分　许多基因：从网络细节到集体响应

### 16. 布尔网络与粗粒化状态

此前的微分方程适合分析少量基因的精细动力学。但当 $N$ 很大时，我们可能首先关心更粗的问题：**一个很小的表达扰动在网络中会消失，还是会传播到越来越多的基因？**

布尔模型令每个基因只有两种状态 $\sigma_i(t)\in\{0,1\}$，并以同步离散时间更新：

$$
\boxed{
\sigma_i(t+1)=F_i\bigl(\sigma_{j_1}(t),\ldots,\sigma_{j_K}(t)\bigr).
}
\tag{68}
$$

固定的布尔规则和连线确定有限状态空间上的吸引子。在这类**同步确定性**模型中，轨道最终必须进入固定点或周期轨道，但是否可以把它们对应到真实细胞状态，不能仅凭模型构造本身保证。

### 17. Derrida–Pomeau 退火近似

考虑一个简单随机布尔网络系综：

- 每个节点随机选择 $K$ 个输入；
- 每个节点真值表中的每个输出独立地以概率 $p$ 取 1，以概率 $1-p$ 取 0；
- 用退火近似（annealed approximation）忽略长期演化中固定连线、固定真值表造成的相关性（可理解为每步对相关随机结构重新抽样）。

取两个彼此接近的初始状态 $\boldsymbol\sigma$ 和 $\boldsymbol\sigma'$，用归一化的汉明距离（Hamming distance，即两个状态中不同分量所占的比例）

$$
d_t=\frac1N\sum_{i=1}^N
\mathbf 1\bigl[\sigma_i(t)\ne\sigma_i'(t)\bigr]
\tag{69}
$$

描述两种状态的差异。

**第一步：** 对一个节点，它的 $K$ 个输入在两个网络副本中完全相同的概率近似为

$$
(1-d_t)^K.
$$

因此至少一个输入不同的概率为

$$
1-(1-d_t)^K.
$$

**第二步：** 如果输入组合不同，它们对应随机真值表中的两个不同表项。两个输出一个为 1、另一个为 0 的概率是

$$
p(1-p)+(1-p)p=2p(1-p).
$$

**第三步：** 将两种独立概率相乘，得到下一步的期望差异比例

$$
\boxed{
d_{t+1}=2p(1-p)\left[1-(1-d_t)^K\right].
}
\tag{70}
$$

这就是退火近似下的扰动传播映射 [6]。

对很小的扰动 $d_t\ll1$，有

$$
1-(1-d_t)^K=K d_t+O(d_t^2),
$$

所以

$$
\boxed{
d_{t+1}=\Lambda d_t+O(d_t^2),
\qquad \Lambda=2Kp(1-p).
}
\tag{71}
$$

于是

$$
\boxed{
\begin{cases}
\Lambda\lt 1:&\text{小扰动平均收缩（ordered）},\\
\Lambda\gt 1:&\text{小扰动平均传播},\\
\Lambda=1:&\text{该近似下的临界分界}.
\end{cases}}
\tag{72}
$$

尤其当 $p=1/2$，

$$
\Lambda=K/2,
\qquad \boxed{K_c=2.}
\tag{73}
$$

这个结果的价值在于：不必给每一个基因指定生化参数，仍能就一个明确的集体动力学性质给出可解析的相边界。

### 18. 平均入度与传播阈值的适用范围

对于一个具体 布尔规则 $F_i$，可以定义其对第 $j$ 个输入的敏感性：

$$
s_{ij}=\Pr_{\boldsymbol\sigma}
\bigl[F_i(\boldsymbol\sigma)\ne F_i(\boldsymbol\sigma^{(j)})\bigr],
\tag{74}
$$

其中 $\boldsymbol\sigma^{(j)}$ 表示只翻转第 $j$ 个输入的组合，概率是对指定的输入分布取平均。节点的总敏感性为

$$
s_i=\sum_{j\in\mathrm{inputs}(i)}s_{ij}.
\tag{75}
$$

对独立随机真值表，每个 $s_{ij}$ 的系综平均为 $2p(1-p)$，故 $\langle s_i\rangle=2Kp(1-p)$；式 (71) 正是这一具体选择的结果。

真实逻辑函数可能具有渠化（canalization）：某个输入取特定值时，输出便被确定，其他输入不再影响它。节点的度分布也可能不均匀，更新可能异步进行。因此，**具有相同 $K$ 的网络，不一定有相同的扰动传播性质**。传播阈值还受逻辑敏感性、连线结构与动力学相关性影响。

必须区分三个层次：

1. **严格模型结论**：在指定随机系综与退火近似下，有 $K_c=2$；
2. **理论启发**：基因网络可能在稳健性和响应灵敏度之间面临权衡；
3. **需要额外证据的生物学主张**：“真实 GRN 通常自组织到临界点”。仅凭前两个结论，不能得到第三个结论。

---

## 第八部分　补充例子：无需反馈也能实现时间过滤

### 19. 一致前馈回路与持续信号检测

无反馈的前馈模体（feed-forward motif）也能对信号持续时间作出选择性响应。模体指网络中反复出现的小型连接结构；这里考虑的是一致前馈回路，即两条输入到输出的路径具有相同的调控符号。

考虑输入 $u$ 一方面直接激活输出 $Z$，另一方面激活中间变量 $Y$，而 $Y$ 也激活 $Z$。设输出采用近似 与门（AND）规则，即需要 $u$ 和 $Y$ 都足够高，$Z$ 才显著表达。

为简化，令输入由 0 突然变成 1，并令

$$
\dot y=\beta_y u-\gamma_y y,
\qquad y(0)=0.
\tag{76}
$$

在 $u=1$ 的持续时间内，解为

$$
y(t)=y_\infty(1-e^{-\gamma_y t}),
\qquad y_\infty=\frac{\beta_y}{\gamma_y}.
\tag{77}
$$

若输出需要 $y\gt y_c$，其中 $0\lt y_c\lt y_\infty$，则启动时间为

$$
\boxed{
t_{\rm on}=\frac1{\gamma_y}
\ln\!\left(\frac{y_\infty}{y_\infty-y_c}\right).
}
\tag{78}
$$

如果输入脉冲长度 $T\lt t_{\rm on}$，中间变量来不及上升到阈值，输出就几乎不启动；反之，持续输入可以启动输出。由于输出同时依赖输入的直接通路，一旦 $u$ 关闭，输出产生率可以迅速停止。

这种机制可以滤掉短脉冲，称为**持续信号检测（persistence detection）** [7]。它依赖具体的组合逻辑和中间变量的时间尺度，不能仅由前馈回路的连线推出。

---

## 第九部分　网络结构与动力学机制

### 20. 几类最小机制的比较

| 网络/过程 | 最小理论描述 | 可能出现的行为 | 必须额外检查的条件 |
|---|---|---|---|
| 自激活 | 一维非线性 ODE | 双稳态、记忆、迟滞 | 非线性增益是否足以产生三个固定点；控制参数是否真的跨越双向分岔 |
| 双基因互抑 | 二维 ODE | 对称性破缺、互斥表达状态 | 非对称稳定态是否存在；是否存在其他慢过程 |
| 三基因互抑环 | 三维 ODE | 振荡失稳、周期动力学 | 环路增益、相位滞后、Hopf 性质和极限环稳定性 |
| 随机出生–死亡 | 主方程 | 分子数涨落 | 平均分子数、反应事件与相关时间 |
| 随机双稳态 | Langevin / WKB | 噪声诱导切换、亚稳态记忆 | 非零逃逸路径、势垒或大偏差作用量 |
| 随机布尔网络 | 离散随机网络 | 扰动收缩或扩散 | 更新规则、逻辑敏感性、连线相关性、近似前提 |
| 一致前馈 | 有限时间的前馈滤波 | 输入持续时间识别 | 组合逻辑与中间变量时间尺度 |

### 21. 结论及其条件

1. **一个基因调控图不是一个动力学模型。** 调控函数、降解机制与时间尺度才把图变成预测。
2. **正反馈是多稳态的结构来源，但不是多稳态的充分条件。** 强度和非线性必须达到阈值。
3. **负反馈可以稳定，也可以振荡。** 决定其命运的是环路增益和相位滞后，而不仅是负边的数目。
4. **噪声既带来局部表达不确定性，也能驱动吸引子间极其罕见的跳变。** 典型涨落与罕见切换由不同的近似控制。
5. **大规模网络可以有统计规律，但特定随机系综的临界点不等于真实 GRN 的普适定律。**

### 22. 几个容易混淆的判断

**问题 A：** 观察到某个细胞表达状态持续数天，能否断言其调控网络有双稳态？

不能。缓慢弛豫、长寿命上游信号、外部环境的持久差异以及真正的多稳态，都可能产生长时间状态保持。需要进一步区分状态保持机制和状态切换性质。

**问题 B：** 某个三基因负反馈系统的 Hill 指数约为 1.5，是否可以断言它不能振荡？

不能。只能说**式 (30) 的三变量、瞬时抑制、对称一阶模型**在 $n\le2$ 时不会通过我们推导的机制发生 Hopf 失稳。如果系统具有显著转录—翻译滞后、额外状态或其他动态相互作用，结论需要重新分析。

**问题 C：** 为什么不把现实中所有已知调控关系都加进模型？

因为解释力不随细节数量单调增加。模型的价值来自它能否识别一个确定机制，并给出可区分的动力学预言。添加生物学细节只有在它改变了所研究问题的关键机制、临界条件或可观测预测时才有明确收益。

---

## 补充推导

### 1. 单基因自激活的临界点

对于 $\dot x=\alpha x^n/(1+x^n)-x$，取 $b=0,n\gt 1$，正固定点出现的最小 $\alpha$ 为

$$
\alpha_c=n(n-1)^{-(n-1)/n}.
$$

**推导：** 在正固定点上 $\alpha=x+x^{1-n}$。对其求导，极小值条件为 $x^n=n-1$，代回即得。进一步求二阶导数可知对应的是极小值。若 $\alpha\lt \alpha_c$，没有正固定点；若 $\alpha\gt \alpha_c$，有两个正固定点，且通过线性化可证明小的正固定点不稳定、大的稳定。

### 2. 互抑制开关的临界模态

在式 (20) 中取一般 $n\gt 1$，设对称稳态为 $s$。它对反对称扰动的稳定性由下面的特征值决定。

**推导：** 因为 $\alpha=s(1+s^n)$，对称雅可比矩阵 的非对角元为 $-g$，其中 $g=ns^n/(1+s^n)$。反对称特征值为 $\lambda=-1+g$，因此临界条件是 $g=1$，即 $s^n=1/(n-1)$。所以

$$
\boxed{
\alpha_c^{\mathrm{toggle}}=
\frac{n}{n-1}\left(\frac1{n-1}\right)^{1/n}
=n(n-1)^{-(1+1/n)}.
}
$$

$n=2$ 时回到 $\alpha_c=2$。注意这是对称分支的线性失稳条件，并不自动给出所有可能的不对称模型分支结构。

### 3. 三基因负反馈的相位模态

Repressilator 在 Hopf 阈值附近的一对临界特征向量可以选为

$$
(1,e^{2\pi\mathrm i/3},e^{4\pi\mathrm i/3})
$$

及其复共轭。

**推导：** 环形移位算符作用于上述向量，只乘上 $e^{-2\pi\mathrm i/3}$。故它是线性化矩阵的特征向量，并具有式 (35) 的特征值。向量的三个分量相位依次相差 $2\pi/3$，所以临界振荡不是三个基因同步升降，而是沿调控环依次变化的相位波。

### 4. 从反应跳跃到化学噪声

考虑跃迁 $N\to N+1$ 的速率 $\Omega B(N/\Omega)$，以及 $N\to N-1$ 的速率 $\Omega R(N/\Omega)$。对 $x=N/\Omega$，漂移和扩散系数由每次跃迁的增量与发生率共同决定。

**推导：** 单次跃迁为 $\pm1/\Omega$。单位时间内的一阶变化量为 $(\Omega B-\Omega R)/\Omega=B-R$；二阶增量为 $(\Omega B+\Omega R)/\Omega^2=(B+R)/\Omega$。因此得到式 (61)、(62)。但是大偏差需要完整的跃迁生成函数 $B(e^q-1)+R(e^{-q}-1)$，只保留 $q$ 的二阶项不能一般性地还原 $q=\ln(R/B)$。

### 5. 布尔网络的临界条件

对于有偏随机真值表 $p\ne1/2$ 的网络，临界入度 $K_c$ 由扰动传播率决定。

**推导：** 由 $\Lambda=2Kp(1-p)$，得到

$$
\boxed{K_c=\frac1{2p(1-p)}.}
$$

每条输入边的平均敏感性为 $2p(1-p)$，依赖的是随机真值表输出相互独立的假设。具体调控函数可能具有渠化，网络也可能存在度异质性、相关性和异步更新，因此传播阈值不能仅由 $K$ 与 $p$ 确定。

---

## 参考文献与延伸阅读

下列文献对应本文讨论的具体机制。

1. **Thomas 反馈回路原则**：R. Thomas, “The role of feedback circuits: Positive feedback circuits are a necessary condition for positive real eigenvalues of the Jacobian matrix,” *Berichte der Bunsengesellschaft für physikalische Chemie* **98**, 1148–1151 (1994). [DOI: 10.1002/bbpc.19940980916](https://doi.org/10.1002/bbpc.19940980916)。另见 C. Soulé, “Mathematical approaches to differentiation and gene regulation,” *Comptes Rendus Biologies* **329**, 13–20 (2006), [DOI: 10.1016/j.crvi.2005.10.002](https://doi.org/10.1016/j.crvi.2005.10.002)。关于必要条件的严格表述应参考具体形式化假设。
2. **回路、状态与吸引子的理论综述**：D. Thieffry, “Dynamical roles of biological regulatory circuits,” *Briefings in Bioinformatics* **8**, 220–225 (2007). [期刊页面](https://academic.oup.com/bib/article/8/4/220/221430)。
3. **双稳态合成开关**：T. S. Gardner, C. R. Cantor and J. J. Collins, “Construction of a genetic toggle switch in *Escherichia coli*,” *Nature* **403**, 339–342 (2000). [DOI: 10.1038/35002131](https://doi.org/10.1038/35002131)。
4. **三基因振荡器**：M. B. Elowitz and S. Leibler, “A synthetic oscillatory network of transcriptional regulators,” *Nature* **403**, 335–338 (2000). [DOI: 10.1038/35002125](https://doi.org/10.1038/35002125)。
5. **单细胞表达噪声**：M. B. Elowitz, A. J. Levine, E. D. Siggia and P. S. Swain, “Stochastic gene expression in a single cell,” *Science* **297**, 1183–1186 (2002). [DOI: 10.1126/science.1070919](https://doi.org/10.1126/science.1070919)。
6. **随机布尔网络的退火近似**：B. Derrida and Y. Pomeau, “Random networks of automata: A simple annealed approximation,” *Europhysics Letters* **1**, 45–49 (1986). [DOI: 10.1209/0295-5075/1/2/001](https://doi.org/10.1209/0295-5075/1/2/001)。
7. **网络模体与信号处理**：U. Alon, “Network motifs: theory and experimental approaches,” *Nature Reviews Genetics* **8**, 450–461 (2007). [DOI: 10.1038/nrg2102](https://doi.org/10.1038/nrg2102)。另见 U. Alon, *An Introduction to Systems Biology: Design Principles of Biological Circuits*，讨论生物回路的设计原理。

---
