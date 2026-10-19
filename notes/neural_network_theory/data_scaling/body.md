> **V · Scaling Laws · Note 18**
>
> ← [Note 17](/notes/neural_network_theory/time_scaling/) · [系列目录](/notes/neural_network_theory/) · [Note 19](/notes/neural_network_theory/model_scaling/) →

**可解性层次：** 高斯序列代理的可解标度。

## 1. 数据幂律需要说明数据在观测什么

给模型更多样本，测试误差可能近似 $L(D)-L_\infty\propto D^{-\beta}$。指数不能仅从模型参数量决定；还要知道目标尾部、观测噪声、采样几何以及正则化怎样随 $D$ 变化。

本篇采用 Note 07 的高斯序列代理，完整推导一个数据指数。它给出透明的偏差–方差机制；真实随机输入样本的核回归还需控制经验算子和本征函数统计。

## 2. 可解模型与目标正则性

设 $\eta_i\asymp i^{-p}$、$a_i^2\asymp i^{-q}$，其中 $p\gt 1,q\gt 1$。前者使核迹有限，后者保证目标 $L^2$ 能量有限。观察

$$
z_i=a_i+\frac\sigma{\sqrt D}\xi_i,
\qquad\hat a_i=\frac{\eta_i}{\eta_i+\lambda}z_i.
$$

风险为

$$
E(D,\lambda)=\sum_i\left(\frac\lambda{\eta_i+\lambda}\right)^2a_i^2
+\frac{\sigma^2}{D}\sum_i\left(\frac{\eta_i}{\eta_i+\lambda}\right)^2.
$$

记有效截止模态 $k_\lambda\asymp\lambda^{-1/p}$。大于它的模态被强烈抑制，小于它的模态主要受观测噪声限制。这里 $i$ 是排序后的模态编号；若原始系统有多维频率及简并度，必须先把简并度计入编号谱，才能使用下面的指数。

## 3. 逐项估计偏差，不跳过 ridge 饱和

高模态尾部满足

$$
B_{\rm tail}\asymp\sum_{i\gt k_\lambda}i^{-q}
\asymp k_\lambda^{1-q}.
$$

低模态残余偏差则为

$$
B_{\rm low}\asymp\lambda^2\sum_{i\leq k_\lambda}i^{2p-q}.
$$

若 $1\lt q\lt 2p+1$，求和由上端主导，两部分同阶 $k_\lambda^{1-q}$。若 $q=2p+1$，出现 $\lambda^2\log k_\lambda$；若 $q\gt 2p+1$，低模态求和收敛，偏差变成 $\lambda^2\asymp k_\lambda^{-2p}$，这是普通 ridge 的平滑度饱和。

因此不能对任意光滑目标都一直使用 $B\sim k^{1-q}$。换用更高阶谱滤波器或硬截止，饱和性质也会改变。

### 连续谱积分还能保留领先常数

取 $\eta_i=c i^{-p}$、$a_i^2\sim A i^{-q}$，令 $k=(c/\lambda)^{1/p}$。在不饱和区域，对完整滤波作积分，

$$
B\sim A k^{1-q}\int_0^\infty
\frac{u^{2p-q}}{(1+u^p)^2}du,
\qquad
V\sim\frac{\sigma^2k}{D}\int_0^\infty
\frac{du}{(1+u^p)^2}.
$$

第一个积分在零端要求 $q\lt 2p+1$，在无穷远要求 $q\gt 1$，重新给出偏差估计的有效区间。饱和的本质，是误差主导位置从高编号尾部转移到有限个低模态，而不是积分方法失效。

第二个积分为 $\Gamma(1/p)\Gamma(2-1/p)/p$，与 Note 07 的噪声有效维数相同。因此有效截止模态的图像可以精确到领先常数，而不只匹配指数。

## 4. 方差与最优数据指数

方差的低模态贡献约为 $\sigma^2 k_\lambda/D$。高模态部分为

$$
V_{\rm tail}\asymp\frac{\sigma^2}{D\lambda^2}
\sum_{i\gt k_\lambda}i^{-2p}
\asymp\frac{\sigma^2}{D}k_\lambda,
$$

所以两者同阶。在不饱和区间 $1\lt q\lt 2p+1$，

$$
E\asymp C_bk^{-(q-1)}+C_v\frac{\sigma^2 k}{D}.
$$

求导或配平得到

$$
k_*\asymp(D/\sigma^2)^{1/q},\qquad
\lambda_*\asymp(D/\sigma^2)^{-p/q},
$$

$$
\boxed{E_*(D)\asymp(\sigma^2/D)^{(q-1)/q},
\qquad\beta=\frac{q-1}{q}.}
$$

在饱和区间，偏差指数用 $2p$ 替换 $q-1$，给出 $\beta=2p/(2p+1)$，临界处带对数修正。这个可解代理里 $p$ 控制所需 $\lambda$，在不饱和区间最终的最优 $D$ 指数主要由目标能量尾部决定；这依赖当前用 $a_i$ 定义目标的方式，换成 RKHS 标准化系数时表达会重写。

### 与数值图一致的完整常数：$p=q=2$

取 $\eta_i=i^{-2}$、$a_i^2=A i^{-2}$，无限模态单位能量归一化给出 $A=6/\pi^2$。两个积分都可算出：$\int_0^\infty u^2/(1+u^2)^2du=\int_0^\infty(1+u^2)^{-2}du=\pi/4$。于是

$$
E(k)\sim\frac{A\pi}{4k}+\frac{\pi\sigma^2k}{4D},
\quad k_*\sim\sqrt{\frac{AD}{\sigma^2}},
\quad\lambda_*\sim\frac{\sigma^2}{AD},
$$

$$
\boxed{E_*(D)\sim\frac{\sqrt6}{2}\frac\sigma{\sqrt D}.}
$$

图中 $\sigma=0.3$，领先系数约为 0.3674。在远离低模态离散效应和高模态截断的区间，$\sqrt D\,E(D)$ 应接近这个数。只匹配指数通常不够约束理论，常数与交叉尺度提供了更严格的比较。

### 函数正则性与维数怎样进入指数

若一个 $d_x$ 维问题按频率排列后的目标尾部满足 $\sum_{i\gt k}a_i^2\asymp k^{-2r/d_x}$，当前模型中便有 $q-1=2r/d_x$，于是

$$
\beta=\frac{2r}{2r+d_x}.
$$

这里 $r$ 是由尾部衰减定义的有效正则性。这与 Sobolev 非参数回归的速率结构相同：更平滑的目标高阶能量衰减更快，更多输入维数使同一分辨率包含更多模态。

这不意味着可以从现实网络的数据指数直接反推出输入的真实维数。数据可能位于低维流形，核模态可能不是普通 Fourier 波，目标也可能只依赖少量特殊方向。恰在幂律边界上的系数序列，还可能需要对数修正才能属于相应的严格 Sobolev 空间。指数对应的是一组尾部与观测假设，并不能独立鉴定数据几何。

### Source 条件解释 ridge 的饱和

另一种目标描述为 $a_i=\eta_i^s b_i$，$\sum_i b_i^2\lt \infty$。偏差是 $\sum_i[\lambda\eta_i^s/(\eta_i+\lambda)]^2b_i^2$。对于 $0\lt s\leq1$，令 $\eta_i=\lambda u$，乘子可由常数乘 $\lambda^s$ 控制，故 $B\lesssim\lambda^{2s}$。当 $s\gt 1$，有限低模态的残余收缩使普通 ridge 的典型领先偏差停在 $\lambda^2$。

这与按 $q$ 分类的饱和，是同一个滤波机制的两种写法。谱硬截断、迭代正则化或高阶滤波器可以改变目标光滑度与偏差的关系，因此数据指数既反映目标，也反映估计器。

## 5. 噪声、早停和无噪声极限

若评价有噪声测试标签的 MSE，应加上 $\sigma^2$ 的不可约底噪；拟合幂律时要对超额风险做拟合。若学习的是确定性函数而只有训练观测含噪，则底噪的含义另行指定。

早停滤波 $1-e^{-\eta_it}$ 也产生 $k_t\asymp t^{1/p}$；但有限数据中的噪声滤波和经验谱必须重新计算，不能只把 $D$ 换成 $t$ 就宣称得到同一学习曲线。

尤其在 $\sigma=0$ 的高斯序列模型中，所有 $z_i=a_i$ 已被直接无噪声观测，$D$ 不再控制信息量，令 $\lambda\to0$ 就能恢复目标。它不能预测真实有限点采样的无噪声数据指数。若另定义“只观察前 $D$ 个模态”的任务，尾部给出 $D^{-(q-1)}$，但那已经是另一个观测模型。

### 新样本、重复观测与相关数据

在序列实验中，$D$ 被定义为每个谱系数噪声方差缩小到 $\sigma^2/D$ 的观测预算。真实数据中，增加独立输入还会改善方向覆盖与经验算子估计；重复相同输入主要平均该点的观测噪声，未必增加新的可辨识模态。

若样本具有时间相关，均值估计的方差包含自协方差和，而不是简单的 $1/D$。不同模态还可能有不同相关时间，因此未必存在一个有效样本数替换所有 $D$。去重、数据质量与采样分布变化，都可能同时改变噪声和目标，不能只按条目数比较同一学习曲线。

处理 token 数能反映计算量，却不自动等于独立信息量。数据大量重复时，$BD^{-\beta}$ 的统计含义必须重新定义；Note 20 的成本优化只有在这种操作定义固定后才有可比较的指数。

## 6. 数值结果

![高斯序列模型的数据风险和最优有效模态数。](figures/experiment.svg)

*图 1 · 自制可解序列实验。* 使用 100000 个模态、$p=q=2$、归一化目标能量、$\sigma=0.3$。每个 $D$ 在对数网格上最小化解析期望风险。左图比较所得风险与 $D^{-1/2}$；右图显示最优 $\lambda$ 对应的有效模态数。此次生成的中间区间斜率为 $-0.4998$。有限模态截断和超参数网格限制了可用拟合窗口。

复现命令为 `python experiments.py --note 18`。在 $p=2$ 时，$q=4$ 对应不饱和指数 $3/4$；$q=7$ 则进入 ridge 饱和区间，指数为 $4/5$。这些是风险式的解析预测，图中实际计算采用 $p=q=2$。

这里按已知教师选取最优 $\lambda$，是一条 oracle 基线。真实学习程序需要验证集或其他可执行规则，选择误差与计算预算不在当前解析风险里。

## 7. 开放问题

从谱代理走到真实数据，要研究本征函数是否集中、采样是否覆盖关键区域、输入分布是否随数据规模变化，以及特征学习是否使 $p,q$ 本身改变。这些因素都可能改变指数或产生交叉区间。

**相关文献。** [Bordelon、Canatar 与 Pehlevan](https://proceedings.mlr.press/v119/bordelon20a.html)讨论核学习曲线的谱依赖；[Bahri 等，Explaining Neural Scaling Laws](https://arxiv.org/abs/2102.06701)比较不同幂律机制。本文推导的指数属于明确的序列观测和滤波规则，不能直接贴到语言模型损失上。

---

← [Note 17](/notes/neural_network_theory/time_scaling/) · [系列目录](/notes/neural_network_theory/) · [Note 19](/notes/neural_network_theory/model_scaling/) →
