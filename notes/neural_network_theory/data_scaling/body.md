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

## 5. 噪声、早停和无噪声极限

若评价有噪声测试标签的 MSE，应加上 $\sigma^2$ 的不可约底噪；拟合幂律时要对超额风险做拟合。若学习的是确定性函数而只有训练观测含噪，则底噪的含义另行指定。

早停滤波 $1-e^{-\eta_it}$ 也产生 $k_t\asymp t^{1/p}$；但有限数据中的噪声滤波和经验谱必须重新计算，不能只把 $D$ 换成 $t$ 就宣称得到同一学习曲线。

尤其在 $\sigma=0$ 的高斯序列模型中，所有 $z_i=a_i$ 已被直接无噪声观测，$D$ 不再控制信息量，令 $\lambda\to0$ 就能恢复目标。它不能预测真实有限点采样的无噪声数据指数。若另定义“只观察前 $D$ 个模态”的任务，尾部给出 $D^{-(q-1)}$，但那已经是另一个观测模型。

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
