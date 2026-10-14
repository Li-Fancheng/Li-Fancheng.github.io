# From Single-Cell Division Variability to Population Growth {#from-single-cell-noise-to-population-growth}
## Size Inheritance, Cell-Cycle Regulation, and Lineage Correlations {#cell-cycle-regulation-temporal-correlations-and-stochastic-dynamics}

> Personal notes for Introduction to Quantitative Biology, Westlake University, Fall 2026.

> **Central question.** Can fluctuations in single-cell division times translate into differences in long-term population growth? How does cell-size regulation affect this translation?

<!-- Legacy heading anchors below preserve links to relocated sections. -->

## Overview {#main-lecture}

<div id="0-the-physical-problem-and-the-roadmap"></div>

<div id="notation"></div>

<div id="lecture-roadmap-and-learning-objectives"></div>

Even cells with the same genetic background in the same environment do not divide at the same time: some generations are short and others long. Meanwhile, cell sizes fluctuate while their distribution can remain stable. Single-cell tracking lets us measure both phenomena and poses a prediction problem: **is the distribution of single-cell generation times sufficient to infer how fast a population multiplies?**

The most direct estimate treats the mean generation time measured along a lineage as the population doubling time. Yet lineages contribute different numbers of descendants, and earlier-dividing ancestors usually leave more offspring. [Hashimoto et al. (2016)](https://doi.org/10.1073/pnas.1519412113) observed population doubling times shorter than mean generation times; the experimental figure is in [Appendix A.2](#a2-experimental-doubling-time-comparison). This observation calls for an explicit statistical model connecting single-cell averages to population growth.

A second question is essential: **can an early division become a timing lead that persists over many generations?** Early division may reflect faster biomass accumulation, or division at a smaller size. In the latter case, daughters inherit a smaller birth size, and size feedback can affect when they divide next. The duration of one generation therefore involves growth, division regulation, and size inheritance.

We first use experimental size stability to establish the relation between growth and division, then use independent division times as a baseline for understanding how population sampling favors faster histories. The adder then shows how size deviations propagate to later generations, distinguishing local timing fluctuations from differences that accumulate over many generations. Finally, we summarize how size regulation constrains whole lineage histories. The main analysis holds every cell's size-growth rate fixed to isolate division timing; variability in growth rates is an extension discussed at the end.

**Size homeostasis** describes correction of size deviations; **stationarity** means that a statistical distribution does not change with generation. They are related but distinct. We write cell size as $v$ and use volume for the conservation argument.

We distinguish the **population growth rate** $\Lambda=\lim_{t\to\infty}t^{-1}\ln N(t)$, when this limit exists, from the common **single-cell size growth rate** $\lambda$ in $\dot v=\lambda v$. Branching calculations first determine the exponent of $M(t)=\mathbb E[N(t)]$; identifying it with the exponent of a typical population requires the branching regularity conditions discussed in [Appendix C](#appendix-c). The volume argument instead applies directly to a realized population.

| Symbol | Meaning |
|---|---|
| $v_n$, $v_n^{\mathrm d}$ | Birth and division sizes in generation $n$ |
| $\tau_n$, $T_n=\sum_{j=0}^{n-1}\tau_j$ | One generation time and accumulated time |
| $\mu=\mathbb E[\tau_n]$ | Mean chronological generation time |
| $\Delta$ | Mean added size in the adder |
| $C_y(k)$, $S_y(\omega)$ | Stationary covariance and two-sided PSD, with $\omega$ in radians per generation |

A **chronological lineage** follows one daughter chosen with probability $1/2$ at each division, without weighting by her eventual descendants. A **retrospective lineage** traces the ancestors of a cell sampled from a later population snapshot. These procedures generally produce different statistics.

---

<div id="1-why-an-independent-division-timer-does-not-regulate-size"></div>

## 1. Experimental starting point: variable division times, stable sizes {#1-cell-size-homeostasis-and-the-failure-of-a-timer}

Single-cell tracking records division times, birth and division sizes, and growth trajectories between divisions. Initially small or large cells recover toward normal sizes, and generation time depends on birth size. The observations of [Taheri-Araghi et al. (2015)](https://doi.org/10.1016/j.cub.2014.12.009) motivate comparing division rules that respond differently to inherited birth size. Explaining generation-time variability requires accounting for this response.

<div id="11-exponential-growth-between-divisions"></div>

### 1.1 Exponential growth and symmetric division {#problem-1-what-does-it-take-to-remove-size-memory}

For some bacterial growth conditions, a useful approximation is exponential growth at the same constant rate $\lambda>0$ for every cell:

$$
\frac{dv}{dt}=\lambda v,\qquad v(t)=v_n e^{\lambda t}.
$$

If division occurs after $\tau_n$ and each daughter receives half the mother's volume, then

$$
v_n^{\mathrm d}=v_n e^{\lambda\tau_n},\qquad
v_{n+1}=\frac12v_n e^{\lambda\tau_n}.\tag{1.1}
$$

This relation separates the division rule, which chooses $\tau_n$, from the growth and partition kinematics. At the same birth size and the same $\lambda$, earlier division means a smaller division size and smaller daughters. A change in division timing therefore changes the next generation's initial state.

<div id="12-the-timer-becomes-a-random-walk-in-log-size"></div>

### 1.2 Independent division times omit feedback from inherited size {#12-a-timer-is-a-random-walk-in-log-size}

A size-independent timer draws i.i.d. positive times with mean $\mu$ and variance $\sigma_\tau^2$. Set $y_n=\ln v_n$. Taking logarithms of (1.1) and iterating gives

$$
y_{n+1}-y_n=\lambda\tau_n-\ln2,\qquad
y_n-y_0=\sum_{j=0}^{n-1}(\lambda\tau_j-\ln2).\tag{1.2}
$$

Independence now makes both consequences transparent:

$$
\mathbb E[y_n-y_0]=n(\lambda\mu-\ln2),\qquad
\operatorname{Var}(y_n-y_0)=n\lambda^2\sigma_\tau^2.\tag{1.3}
$$

Tuning $\mu=\ln2/\lambda$ removes the drift but leaves unbounded variance. **Zero drift is not restoring feedback.** Even the noise-free balanced timer gives $v_{n+1}=v_n$, preserving an initial size perturbation forever. A restoring map instead shrinks the perturbation by a factor $|a|<1$ each generation. Eliminating noise prevents wandering; feedback corrects an existing deviation.

![Timer trajectories wander while restoring feedback bounds size fluctuations.](figures/01_timer_and_homeostasis.svg)

*Fig 1. Why feedback matters.* Orange trajectories obey the zero-drift log-size timer (1.2); blue trajectories obey an illustrative restoring **log-size** map with $a=1/2$ and the same innovations. The blue map is not the exact additive-size adder. Panel B compares analytic variances from a fixed initial size. All trajectories are synthetic.

### 1.3 Timer, adder, and sizer

Three idealized rules illustrate what feedback means:

| Rule | Division criterion | Noise-free birth-size response near balance |
|---|---|---|
| Timer | Wait a fixed duration, independent of birth size | $v_{n+1}=v_n$: neutral, slope $1$ |
| Adder | Add a fixed size $\Delta$ | $v_{n+1}=(v_n+\Delta)/2$: restoring, slope $1/2$ |
| Sizer | Reach a fixed division size $v_*$ | $v_{n+1}=v_*/2$: no inherited size deviation, slope $0$ |

The slopes describe loss of size memory, rather than unique molecular mechanisms. A noisy adder or sizer can retain a stationary size distribution; a noisy independent timer with the kinematics above cannot. A unified local linear model is developed in [Appendix F](#appendix-f).

![Single-cell microscopy and birth-size dependence testing timer and sizer predictions.](figures/taheri_araghi_2015_f1_highres.png)

*Fig 2. Experimental size control.* Original Figure 1 from [Taheri-Araghi et al. (2015)](https://doi.org/10.1016/j.cub.2014.12.009). In panel D, generation time decreases with birth size while division size increases. These trends challenge a size-independent timer and a fixed-threshold sizer, respectively. Panels A–C give microscopy and growth-condition context.

These observations and models show that generation times cannot be specified freely without regard to birth size. We now temporarily set size feedback aside and use independent division times as a prediction baseline: why does the single-cell mean underestimate the population's multiplication rate in this model?

---

<div id="2-eulerlotka-theory-why-division-time-noise-can-accelerate-a-population"></div>

## 2. Predicting population growth from single-cell statistics: an independent baseline {#2-division-time-noise-and-the-population-growth-paradox}

For a chronological mean generation time $\mu$, the simplest estimate is $\Lambda=\ln2/\mu$. An independent branching model explains what this estimate omits: histories with different division times contribute different numbers of descendants to the later population. The model first answers a statistical prediction question; applying it to size-regulating cells also requires checking the physiological constraints on growth and division.

<div id="21-from-the-renewal-equation-to-the-population-growth-rate"></div>

### 2.1 Euler–Lotka and sampling bias {#problem-2-correlation-without-changing-a-histogram}

Consider a binary branching model with no death: every newborn draws an independent positive lifetime from a density $p(\tau)$ and produces two daughters following the same rule. There is no shared random environment. Let $M(t)$ be the expected number of living descendants of one newborn. Conditioning on its first division gives

$$
M(t)=\mathbb P(\tau>t)+2\int_0^t p(\tau)M(t-\tau)\,d\tau.\tag{2.1}
$$

In an exponentially growing solution, a division delayed by $\tau$ leaves the descendants less time to multiply, supplying a factor $e^{-\Lambda\tau}$. Balancing the factor of two from reproduction gives the **Euler–Lotka equation**:

$$
\boxed{2\mathbb E_p[e^{-\Lambda\tau}]=1.}\tag{2.2}
$$

Its positive root is unique because the left side decreases continuously from $2$ to $0$. [Appendix A](#appendix-a) derives it by Laplace transformation and discusses renewal amplitudes and synchronized divisions.

The same discount explains population-biased sampling. In this independent model, complete ancestral lifetimes in a retrospective sample have density

$$
q(\tau)=2e^{-\Lambda\tau}p(\tau),\qquad \int_0^\infty q(\tau)\,d\tau=1.\tag{2.3}
$$

Short-lived ancestors leave more descendants and receive more weight. This $q$ is a distribution of **completed ancestral lifetimes**, not the age or eventual lifetime distribution of all cells currently alive. In a correlated model, the weighting acts on whole histories and need not reduce to (2.3).

![Renewal histories and growth weighting toward shorter lifetimes.](figures/06_renewal_and_sampling.svg)

*Fig 3. Earlier divisions leave more descendants.* Panel A shows the survivor and two-subtree terms in (2.1). Panel B uses a gamma lifetime density with shape and rate $4$ (mean $1$); its Euler–Lotka rate is $4(2^{1/4}-1)$. Growth weighting shifts the ancestral lifetime distribution toward short times. Dotted lines mark the two means.

<div id="22-jensens-inequality-and-the-population-growth-advantage"></div>

### 2.2 Jensen's inequality and the reference population

A deterministic lifetime $\mu$ gives growth rate $\ln2/\mu$. For variable lifetimes with the **same chronological mean** $\mu<\infty$, convexity yields

$$
\frac12=\mathbb E_p[e^{-\Lambda\tau}]
\ge e^{-\Lambda\mathbb E_p[\tau]}=e^{-\Lambda\mu}.
$$

Taking logarithms,

$$
\boxed{\Lambda\ge\frac{\ln2}{\mu}.}\tag{2.4}
$$

The inequality is strict for a nondegenerate independent lifetime distribution. Earlier divisions contribute disproportionately to later populations; late divisions do not compensate for that benefit at fixed mean time. For weak noise, with $\Lambda\sigma_\tau\ll1$ and a regular fluctuation-amplitude expansion whose higher cumulants are negligible at second order,

$$
\Lambda\simeq\frac{\ln2}{\mu}
+\frac{(\ln2)^2}{2\mu^3}\sigma_\tau^2.\tag{2.5}
$$

The step-by-step expansion is in [Appendix B](#appendix-b); an information-theoretic interpretation is in [Appendix H](#appendix-h). The experimental doubling-time comparison of [Hashimoto et al. (2016)](https://doi.org/10.1073/pnas.1519412113) is retained as Fig 7 in Appendix A. The reference in (2.4) is an independent-lifetime population with mean $\mu$, **not automatically a population whose cells grow in size at rate $\lambda$**.

<div id="24-the-apparent-paradox-volume-conservation"></div>

### 2.3 Cell-number growth and biomass growth: what is held fixed? {#23-the-conservation-constraint}

Equation (2.4) holds mean generation time $\mu$ fixed, without fixing cell sizes or biomass accumulation rates. To distinguish earlier division from faster growth, return to cells growing at one common constant size-growth rate $\lambda$. Define

$$
V_{\mathrm{tot}}(t)=\sum_{i\text{ alive at }t}v_i(t).
$$

Between divisions, $\dot V_{\mathrm{tot}}=\lambda V_{\mathrm{tot}}$. At division, daughters' volumes sum to the mother's volume, so there is no jump. Thus, with no death or biomass loss and no explosion of division events,

$$
\boxed{V_{\mathrm{tot}}(t)=V_{\mathrm{tot}}(0)e^{\lambda t}.}\tag{2.6}
$$

Divisions increase cell number but redistribute existing volume. Write the mean size in a population snapshot as $\overline v_{\mathrm{snap}}(t)=V_{\mathrm{tot}}(t)/N(t)$. Then exactly

$$
\frac{\ln N(t)}t=\lambda+
\frac{\ln V_{\mathrm{tot}}(0)-\ln\overline v_{\mathrm{snap}}(t)}t.
$$

Consequently,

$$
\boxed{\Lambda=\lambda\quad\text{if}\quad
\lim_{t\to\infty}\frac{\ln\overline v_{\mathrm{snap}}(t)}t=0.}\tag{2.7}
$$

Convergence of the snapshot mean to a finite positive value is sufficient; bounded periodic modulation also suffices. A stationary **chronological birth-size** distribution alone does not establish this snapshot condition, because the population weights cells differently.

At fixed $\lambda$, a cell-number exponent $\Lambda>\lambda$ would require the mean snapshot size to decrease exponentially. The independent timer's growth increment relative to the same-mean-generation-time baseline therefore cannot be interpreted directly as a proliferation benefit available to size-homeostatic cells. It also does not establish a growth cost of size regulation.

As a model consistency check, stationarity along a chronological lineage with $\mathbb E|\ln v|<\infty$ gives $\mu=\ln2/\lambda$ by taking expectations of (1.2). Substituting this mean into a noisy independent model would predict $\Lambda>\lambda$, incompatible with (2.7) under its snapshot condition. Section 1 explains why: the independent timer itself cannot maintain size homeostasis. Predicting a real population from single-cell measurements requires describing how division times respond to inherited sizes and how that response affects whole lineage histories.

---

<div id="6-exactly-solvable-example-ii-the-adder-as-a-size-homeostatic-stochastic-map"></div>

## 3. Can an early division create a lasting lead? Size inheritance and the adder {#3-size-feedback-the-adder-and-exact-cancellation}

We now follow the consequences of a timing change. At the same birth size and a fixed growth rate, an earlier division leaves smaller daughters. If the next generation divides after adding a fixed size, its smaller birth size tends to lengthen its generation time. The adder provides a concrete example of this feedback. We first calculate how size deviations recover, then ask whether they can make accumulated time depart from the usual pace over many generations.

<div id="61-definition-of-the-adder"></div>

### 3.1 The stochastic map and recovery of mean size {#problem-3-can-noisy-lifetimes-have-zero-long-run-noise}

An adder divides at

$$
v_n^{\mathrm d}=v_n+\Delta+\eta_n,\qquad
v_{n+1}=\frac12(v_n+\Delta+\eta_n).\tag{3.1}
$$

Here $\eta_n$ are i.i.d. added-size fluctuations of mean zero and finite variance $\sigma_\eta^2$, independent of inherited size. We require $\Delta+\eta_n>0$ so added size and lifetimes are positive; an untruncated Gaussian noise is only a small-fluctuation idealization. Experiments motivate approximately birth-size-independent addition, not by themselves the i.i.d. noise assumption.

![Birth-size-independent added size in bacterial experiments.](figures/taheri_araghi_2015_f2_highres.png)

*Fig 4. Experimental input to the adder.* Original Figure 2 from [Taheri-Araghi et al. (2015)](https://doi.org/10.1016/j.cub.2014.12.009). Mean added size is approximately independent of birth size within each condition, and rescaled conditional distributions collapse. Panels A–C cover *E. coli*, *B. subtilis*, and size mutants.

Taking expectations in (3.1) shows

$$
\mathbb E[v_n]-\Delta=2^{-n}(\mathbb E[v_0]-\Delta).\tag{3.2}
$$

Thus the stationary mean is $\Delta$ and a mean-size perturbation is halved each generation. For $x_n=v_n-\Delta$, the map is $x_{n+1}=x_n/2+\eta_n/2$. In stationarity, independence of the new noise gives $\sigma_x^2=\sigma_x^2/4+\sigma_\eta^2/4$, hence

$$
\operatorname{Var}(v)=\frac{\sigma_\eta^2}{3},\qquad
\rho_v(k)=2^{-k},\quad k\ge0.\tag{3.3}
$$

The autocorrelation follows by iterating the restoring map: future innovations have zero covariance with $x_n$, leaving $\operatorname{Cov}(x_n,x_{n+k})=2^{-k}\sigma_x^2$. Noise fixes the width; feedback fixes memory. [Appendix E](#appendix-e) contains the full calculations and experimental relaxation figure.

<div id="64-connecting-the-adder-to-generation-times"></div>

### 3.2 How generation time records size changes {#32-from-growth-to-generation-time}

Exponential growth and symmetric division require

$$
v_ne^{\lambda\tau_n}=v_n^{\mathrm d}=2v_{n+1}.
$$

Taking logarithms gives $\lambda\tau_n=\ln2+\ln v_{n+1}-\ln v_n$, or

$$
\boxed{\tau_n=\frac{\ln2}{\lambda}
+\frac{\ln v_{n+1}-\ln v_n}{\lambda}.}\tag{3.4}
$$

A generation time is a **difference** of adjacent log sizes. This is an exact identity for positive sizes with common constant $\lambda$ and symmetric division. No adder assumption, stationarity, or weak-noise expansion was used.

<div id="65-the-telescoping-identity-the-central-result-of-the-lecture"></div>

### 3.3 From individual cycles to accumulated time: telescoping {#33-telescoping-the-theoretical-core}

Summing (3.4) from $0$ to $n-1$,

$$
\begin{aligned}
T_n&=\sum_{j=0}^{n-1}\tau_j\\
&=\frac{n\ln2}{\lambda}
+\frac1\lambda\big[(\ln v_1-\ln v_0)
 +(\ln v_2-\ln v_1)+\cdots+(\ln v_n-\ln v_{n-1})\big].
\end{aligned}
$$

Each interior size appears once with each sign, leaving

$$
\boxed{T_n=\frac{n\ln2}{\lambda}
+\frac{\ln v_n-\ln v_0}{\lambda}.}\tag{3.5}
$$

![Interior log-size contributions cancel, leaving only endpoint sizes.](figures/09_telescoping_identity.svg)

*Fig 5. Timing fluctuations become a boundary term.* Blue and purple log-size terms cancel between adjacent generations; only the initial orange and final green terms remain. This cancellation is exact for every $n$.

**The fluctuations of individual generation times need not accumulate over many generations.** For independent times, every generation supplies a new contribution to the timing uncertainty. Here the same inherited size contributes with opposite signs in neighboring generations. Feedback maintains the endpoint sizes within a stationary distribution, so the boundary term cannot develop an extensive variance.

To state this precisely, let $y_n=\ln v_n$ be stationary with finite variance. Equation (3.5) gives

$$
\operatorname{Var}(T_n)
=\frac{2\operatorname{Var}(y)-2\operatorname{Cov}(y_n,y_0)}{\lambda^2}
\le\frac{4\operatorname{Var}(y)}{\lambda^2}.\tag{3.6}
$$

Thus $D_\tau\equiv\lim_{n\to\infty}\operatorname{Var}(T_n)/n=0$; mixing is not needed for this bound. If the endpoint covariance tends to zero, the variance tends to $2\operatorname{Var}(\ln v)/\lambda^2$. Stationarity with finite first log moment also gives $\mu=\ln2/\lambda$ exactly. Ergodicity supplies the corresponding long-record mean. Finite birth-size variance alone does not guarantee finite log-size variance near $v=0$.

<div id="66-exact-consequence-for-population-growth"></div>

### 3.4 The exact growth consequence and its conditions

Population growth weights unusually fast **whole histories**, so a variance bound alone is not enough to control their exponential weights. From the exact identity,

$$
\mathbb E[e^{-sT_n}]
=e^{-ns\ln2/\lambda}
\mathbb E\left[\left(\frac{v_n}{v_0}\right)^{-s/\lambda}\right].\tag{3.7}
$$

If the endpoint expectation has logarithm $o(n)$ for $s$ near the selected growth rate, the exponential rate per generation is simply $-s\ln2/\lambda$. Reproduction supplies $2^n$ paths; its exponential increase balances the time discount when $\ln2-s\ln2/\lambda=0$, giving **$\Lambda=\lambda$** under the branching conditions formalized in Section 4 and [Appendix C](#appendix-c).

A transparent sufficient endpoint condition is stationary positive sizes with finite positive and inverse moments $\mathbb E[v^{2s/\lambda}]$ and $\mathbb E[v^{-2s/\lambda}]$ near $s=\lambda$: Cauchy–Schwarz bounds the expectation uniformly, while Jensen bounds it below by $1$. Sizes bounded above and away from zero are an especially simple case. These are sufficient conditions, not automatic consequences of the adder's finite size variance; [Appendix G](#appendix-g) makes the hierarchy explicit.

This lineage explanation complements the realized-population volume argument in Section 2.3. Both give exact growth constraints with their respective additional conditions. A sizer or another homeostatic rule can obey the same telescoping cancellation. Adder-specific results are the slope $1/2$, variance $\sigma_\eta^2/3$, and birth-size ACF $2^{-k}$, not the identity (3.5).

<div id="67-a-counterintuitive-prediction-size-is-positively-correlated-time-is-negatively-correlated"></div>

### 3.5 Positive size memory, negative timing memory

For small relative size fluctuations, $\ln(\Delta+x_n)\simeq\ln\Delta+x_n/\Delta$, so

$$
\delta\tau_n\equiv\tau_n-\frac{\ln2}{\lambda}
\simeq\frac{x_{n+1}-x_n}{\lambda\Delta},\qquad
\rho_\tau(k)\simeq-2^{-(k+1)},\quad k\ge1.\tag{3.8}
$$

Birth sizes remain positively correlated, yet generation times are negatively correlated in this linearized adder. A long cycle tends to produce a larger daughter, which needs less time to add the next increment. This interpretation refers to the stated model; the negative ACF formula is an approximation, whereas telescoping is exact.

![Positive size correlations and negative timing correlations lead to bounded accumulated variance.](figures/04_noise_cancellation.svg)

*Fig 6. Noisy divisions, zero long-run timing noise.* Panel A compares exact adder birth-size correlations with linearized timing correlations (3.8). With matched one-generation timing variance, panel B shows $\operatorname{Var}(T_n)/C_\tau(0)=n$ for independent times and $2(1-2^{-n})$ for the linearized stationary adder. Full covariance calculations are in Appendix E.

---

<div id="4-why-integrated-temporal-correlations-control-growth"></div>

## 4. Which lineage statistics predict population growth? {#4-temporal-correlations-and-long-time-fluctuations}

The adder shows that a one-generation histogram does not tell us whether a fast generation is followed by another fast generation or compensated by a slower one. In general, growth calculations require the statistics of whole histories. At second order in a weak-noise expansion, that information reduces to the integrated covariance across generations. We develop both levels to establish when a correlation function can predict growth.

<div id="3-correlated-division-times-a-generalized-eulerlotka-relation"></div>

### 4.1 Correlated Euler–Lotka {#problem-4-acf-versus-psd}

Let $T_n$ be the time to the generation-$n$ birth along a chronological daughter path. Binary reproduction creates $2^n$ paths. The expected sum of their discounted birth weights is $2^n\mathbb E_{\mathrm{chron}}[e^{-sT_n}]$, even when paths share ancestry. Define

$$
\psi(s)=\lim_{n\to\infty}\frac1n\ln\mathbb E_{\mathrm{chron}}[e^{-sT_n}].
$$

For a locally specified, nonexplosive branching process with finite expected birth counts on bounded time intervals, a well-defined limit $\psi$, and a growth exponent selected by the discounted-birth convergence boundary,

$$
\boxed{\ln2+\psi(\Lambda)=0.}\tag{4.1}
$$

Stationary/ergodic lineage statistics alone are not sufficient: the limit and its connection to a population exponent must also hold. A shared random environment or nonergodic mixture can make expected and typical growth differ. [Appendix C](#appendix-c) gives the many-to-one derivation and precise convergence argument; [Appendix D](#appendix-d) supplies a complete finite-state Markov example.

Independence gives $\mathbb E[e^{-sT_n}]=(\mathbb E[e^{-s\tau}])^n$, recovering (2.2). Correlation changes this whole-history weight. Jensen applied to $T_n$ still yields $\psi(s)\ge-s\mu$ and hence $\Lambda\ge\ln2/\mu$ when (4.1) holds. Unlike independent noisy lifetimes, correlated noisy lifetimes can attain equality: Section 3 showed how.

<div id="41-variance-of-the-accumulated-generation-time"></div>

<div id="43-an-important-qualification"></div>

### 4.2 The integrated covariance

For stationary times with finite variance, define $C_\tau(k)=\operatorname{Cov}(\tau_j,\tau_{j+k})$. Summing all covariance pairs gives

$$
\operatorname{Var}(T_n)=nC_\tau(0)
+2\sum_{k=1}^{n-1}(n-k)C_\tau(k).\tag{4.2}
$$

There are $n$ diagonal terms, and $n-k$ pairs separated by $k$, counted in both orders. If $\sum_{k\ge1}|C_\tau(k)|<\infty$, divide by $n$ and pass to the limit:

$$
\boxed{D_\tau=C_\tau(0)+2\sum_{k=1}^\infty C_\tau(k).}\tag{4.3}
$$

$D_\tau$ is the long-run variance per generation, in units of time squared. It is nonnegative since it is a limit of variances divided by $n$. The one-generation variance measures local noise; the covariance sum measures whether future cycles reinforce or compensate for it. Independent times have $D_\tau=\sigma_\tau^2$, persistence generally increases the sum, and anticorrelations can cancel it. Absolute summability justifies the covariance representation; the direct variance-limit definition can still make sense when the sum does not converge.

With a regular weak-noise expansion of asymptotic cumulants, the independent result (2.5) becomes

$$
\boxed{\Lambda\simeq\frac{\ln2}{\mu}
+\frac{(\ln2)^2}{2\mu^3}D_\tau.}\tag{4.4}
$$

The replacement is physical: the second cumulant of a long history is $\operatorname{Var}(T_n)$, not $n\operatorname{Var}(\tau)$. [Appendix B](#appendix-b) gives the expansion once, including the required interchange of limits and control of higher cumulant rates. Small single-generation variance alone is insufficient when persistence becomes large.

**$D_\tau=0$ removes the second-order weak-noise correction; by itself it does not guarantee exact $\Lambda=\lambda$.** Exact equality in Section 3 uses the stronger endpoint structure and exponential-moment condition, or the independent volume argument. Neither follows from retaining only a quadratic cumulant term.

<div id="7-from-temporal-correlations-to-power-spectra"></div>

### 4.3 Zero frequency is the same mechanism

For centered generation-time fluctuations, use the two-sided PSD convention

$$
S_\tau(\omega)=\sum_{k=-\infty}^\infty C_\tau(k)e^{-i\omega k},
\qquad -\pi\le\omega\le\pi.
$$

With absolutely summable covariances,

$$
\boxed{S_\tau(0)=D_\tau.}\tag{4.5}
$$

The zero-frequency component measures fluctuations that persist across many generations. For $y_n=\ln v_n$, the exact timing identity is $\delta\tau_n=(y_{n+1}-y_n)/\lambda$. Differencing multiplies a Fourier mode by $e^{i\omega}-1$, so for a stationary log-size process with a regular spectral density,

$$
S_\tau(\omega)=\frac{4\sin^2(\omega/2)}{\lambda^2}S_y(\omega).\tag{4.6}
$$

If $S_y$ is finite at zero, differencing removes zero-frequency timing noise. Slowly varying size deviations appear with almost equal values in adjacent cycles and subtract out. Telescoping, integrated covariance cancellation, and suppression at zero frequency are three expressions of this same mechanism. The full AR(1) spectra are in [Appendix F](#appendix-f); finite-frequency oscillatory modes address a separate question in the [Further Topic](#further-topic).

---

<div id="9-synthesis-three-exact-ideas-and-one-controlled-approximation"></div>

<div id="what-remains-biologically-open"></div>

## 5. Returning to biology: which single-cell differences affect long-term growth? {#5-synthesis-and-outlook}

### 5.1 Division variability and long-term proliferation differences

A mean generation time measured along a lineage is generally insufficient to predict the population doubling time. In the independent model, earlier-dividing ancestors leave more descendants, and Euler–Lotka supplies the statistical weighting. In cells with a fixed size-growth rate, however, earlier division also changes daughters' birth sizes. Their subsequent response determines whether the timing differences accumulate.

The adder demonstrates a specific mechanism: size feedback corrects inherited deviations and produces negative generation-time correlations in the small-fluctuation limit. More generally, at constant $\lambda$ with symmetric division, accumulated timing deviations depend only on the initial and final log sizes. Stationary log sizes with finite variance keep the accumulated variance bounded, giving $D_\tau=0$. Individual cycles still vary, but their timing fluctuations have no persistent variance accumulation.

The exact population exponent additionally requires the snapshot mean-size condition, or subexponential endpoint weights and a valid correlated branching-growth relation. Under the respective conditions, $\Lambda=\lambda$. This is a joint constraint from fixed growth rate and size stability; the analysis does not establish a tradeoff between size regulation and proliferation. For more general correlated lifetime models, $D_\tau$ or $S_\tau(0)$ controls the growth correction only within a valid second-order weak-noise expansion. Its vanishing alone cannot establish an exact result at arbitrary noise strength.

### 5.2 From one-generation distributions to lineage histories

A one-generation distribution describes how often short and long cycles occur, without specifying their order along a lineage. The two-state model in [Appendix D](#appendix-d) makes the distinction explicit: strict alternation fixes the total duration of every two generations, whereas persistence produces runs of fast or slow generations. Identical one-generation distributions can therefore yield different time weights for whole histories and different population growth rates. Connecting single cells to populations requires retaining the links between generations.

Size regulation further restricts which timing histories can occur. Under the growth and partition conditions used here, a lineage that keeps accumulating a timing lead or delay must have birth sizes that keep decreasing or increasing, respectively. Homeostasis constrains that long-term accumulation, beyond the width of the size distribution. The adder makes the feedback mechanism explicitly calculable; telescoping extends to other size-control rules under the same growth and partition conditions. The central distinction is that a one-generation distribution describes local fluctuations, while size inheritance and feedback determine whether they can become lasting differences between lineage histories.

### 5.3 From division timing to variable growth capacity

Holding $\lambda$ fixed isolates how division regulation rearranges timing for a given biomass accumulation law. The next biological question is: **what changes when cells accumulate biomass at different rates and pass those differences to their daughters?** When growth rates vary across generations, there is no common proportionality factor between generation time and log-size change, so accumulated time no longer telescopes in the manner used here. Growth-rate differences and their intergenerational correlations then enter the growth problem, beyond the conclusions of the fixed-$\lambda$ model. [Lin and Amir (2020)](https://doi.org/10.1103/PhysRevE.101.012401) studies these effects. How growth-rate variability and random partition asymmetry change the relations used here is discussed in [Appendix G](#appendix-g).


---

<div id="optional-technical-appendices"></div>

## Technical Appendices

These appendices collect complete derivations, solvable models, and conditions for the results. Equation labels use appendix letters; figure numbering continues from the main text.

<div id="appendix-d-the-asymptotic-amplitude-of-the-renewal-solution"></div>

### Appendix A. Renewal theory, population sampling, and asymptotic amplitude {#appendix-a}

This appendix supplies the renewal calculation behind Section 2. The exponent belongs first to the expected population; standard nonexplosive, supercritical age-dependent branching assumptions are needed to identify the typical exponent.

#### A.1 Laplace-transform derivation

To obtain the population growth rate from (2.1), take Laplace transforms.

Write $S(t)=\mathbb P(\tau>t)$ and define

$$
\widehat M(s)=\int_0^\infty e^{-st}M(t)\,dt,
\qquad
\widehat p(s)=\int_0^\infty e^{-s\tau}p(\tau)\,d\tau
             =\mathbb E[e^{-s\tau}],
$$

with $s$ initially chosen large enough for the integrals to converge. Since $S(t)=\int_t^\infty p(\tau)\,d\tau$, its transform is $\widehat S(s)=[1-\widehat p(s)]/s$. The convolution in (2.1) becomes a product, giving

$$
\widehat M(s)=\widehat S(s)+2\widehat p(s)\widehat M(s)
\quad\Rightarrow\quad
\widehat M(s)=\frac{1-\widehat p(s)}{s[1-2\widehat p(s)]}.
$$

The denominator vanishes at a positive rate $s=\Lambda$ satisfying the **Euler–Lotka equation**,

$$
\boxed{2\,\mathbb E[e^{-\Lambda\tau}]=1.}
$$

The function $2\widehat p(s)$ decreases continuously and strictly from $2$ to $0$ on the positive real axis, so this root is unique. It is the leading pole of $\widehat M(s)$; since $1/(s-\Lambda)$ is the Laplace transform of $e^{\Lambda t}$, the renewal asymptotics for the lifetime density assumed here give

$$
M(t)\sim Ae^{\Lambda t},
\qquad
\Lambda=\lim_{t\to\infty}\frac{\log M(t)}{t}.
$$

Here $A=\lim_{t\to\infty}e^{-\Lambda t}M(t)>0$ is the asymptotic amplitude for one initial newborn cell, and $\sim$ means $M(t)/(Ae^{\Lambda t})\to1$. Its amplitude and the qualification for synchronized divisions are calculated below.

#### A.2 Experimental doubling-time comparison {#a2-experimental-doubling-time-comparison}

![Measured generation-time variability, population doubling time, and growth-rate gain in E. coli.](figures/hashimoto_2016_fig03.jpg)

*Fig 7. Experimental mean generation time and population doubling time.* Original Figure 3 from [Hashimoto et al., *PNAS* (2016)](https://doi.org/10.1073/pnas.1519412113). Panel B compares mean generation time with population doubling time; ratios above one indicate faster population doubling. Panel C relates the gain to generation-time variability; panel E compares Euler–Lotka predictions with measurements. Panels A and D show the lifetime distribution and model, respectively. Original figure reproduced without modification; source: [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4812751/#fig03).

**Notation and scope.** In this paper's gain definition, $\lambda=\ln2/\langle\tau\rangle_g$; our $\lambda$ denotes the exponential growth rate of cell size. These definitions must be distinguished when comparing this figure with (2.7). The figure supports the experimental doubling-time comparison; it does not establish that a size-homeostatic population with one common constant size-growth rate violates mass conservation.

#### A.3 The renewal amplitude and lattice exception

Section 2.1 obtains the population growth exponent from the positive root $\Lambda$ of $1-2\widehat p(s)=0$. No pole lies further to the right, since

$$
|2\widehat p(s)|\le2\widehat p(\operatorname{Re}s)<1
\qquad\text{for }\operatorname{Re}s>\Lambda.
$$

For a lifetime density, the inequality $|2\widehat p(\Lambda+i\omega)|<1$ also holds for $\omega\ne0$, excluding other poles on the same vertical line. Moreover,

$$
-2\widehat p'(\Lambda)
=2\mathbb E_p[\tau e^{-\Lambda\tau}]
=\mathbb E_q[\tau]>0,
$$

so the pole is simple. Expanding the transform from Section 2.1 about this pole gives

$$
\widehat M(s)=\frac{A}{s-\Lambda}+O(1),
\qquad
A=\frac{1-\widehat p(\Lambda)}{\Lambda[-2\widehat p'(\Lambda)]}.
\tag{A.1}
$$

Using $\widehat p(\Lambda)=1/2$, the residue becomes

$$
\boxed{
A=\frac{1}{4\Lambda\,\mathbb E_p[\tau e^{-\Lambda\tau}]}
 =\frac{1}{2\Lambda\,\mathbb E_q[\tau]}.
}\tag{A.2}
$$

The nonlattice renewal asymptotics applicable to a lifetime density then give $e^{-\Lambda t}M(t)\to A$. This amplitude is for one initial newborn cell; it need not equal $M(0)=1$, because the initial transient is not described by the asymptotic exponential.

Fixed lifetimes are an exception to constant-amplitude asymptotics. If every cell divides after exactly $\tau_0$, all division events are synchronized and

$$
M(t)=2^{\lfloor t/\tau_0\rfloor},
\qquad
\Lambda=\frac{\ln2}{\tau_0}.
\tag{A.3}
$$

The logarithmic growth rate still tends to $\Lambda$, but $e^{-\Lambda t}M(t)$ oscillates periodically instead of converging to a constant. Thus the noise-free growth exponent used in Section 2.2 remains valid, while the stronger limit defining $A$ requires the nonlattice condition.


### Appendix B. Independent and correlated weak-noise calculations {#appendix-b}

The two expansions below share one idea: expand the logarithm of the chronological-history weight. Equation (2.5) is the independent result; (4.4) is its long-history counterpart.

<div id="23-a-weak-noise-estimate"></div>

#### B.1 Independent lifetime expansion

The inequality establishes faster growth relative to a deterministic baseline with the same mean generation time, but does not give the size of the increment. To estimate the excess, write the mean and variance of the chronological lifetime distribution explicitly as

$$
\mu=\bar\tau=\int_0^\infty \tau p(\tau)\,d\tau,
\qquad
\sigma_\tau^2=\int_0^\infty (\tau-\mu)^2p(\tau)\,d\tau.
$$

We expand around the mean lifetime, so the small quantity in the exponential is $\Lambda(\tau-\mu)$. Assume $\Lambda\sigma_\tau\ll1$ and a regular expansion in the fluctuation amplitude, with third and higher centered moments negligible at second order. Factoring out $e^{-\Lambda\mu}$ and expanding the remaining exponential gives

$$
\begin{aligned}
\mathbb E_p[e^{-\Lambda\tau}]
&=\int_0^\infty p(\tau)e^{-\Lambda\tau}\,d\tau\\
&=e^{-\Lambda\mu}\int_0^\infty p(\tau)e^{-\Lambda(\tau-\mu)}\,d\tau\\
&\simeq e^{-\Lambda\mu}\int_0^\infty p(\tau)
  \left[1-\Lambda(\tau-\mu)+\frac{\Lambda^2}{2}(\tau-\mu)^2\right]d\tau\\
&=e^{-\Lambda\mu}\left(1+\frac{\Lambda^2}{2}\sigma_\tau^2\right).
\end{aligned}
$$

The constant term integrates to $1$, the linear term vanishes because $\int_0^\infty(\tau-\mu)p(\tau)\,d\tau=0$, and the quadratic term gives the variance. Taking the logarithm and using $\ln(1+x)=x+O(x^2)$ with $x=\Lambda^2\sigma_\tau^2/2$ yields

$$
\ln\mathbb E_p[e^{-\Lambda\tau}]
\simeq-\Lambda\mu+\ln\left(1+\frac{\Lambda^2}{2}\sigma_\tau^2\right)
\simeq-\Lambda\mu+\frac{\Lambda^2}{2}\sigma_\tau^2.
$$

Here $x$ is already second order in the fluctuation amplitude, so the additional $x^2$ term from the logarithm is fourth order. Substituting into the logarithm of Euler–Lotka gives the second-order cumulant expansion

$$
0=\ln2+\ln\mathbb E_p[e^{-\Lambda\tau}]
\simeq\ln2-\Lambda\mu+\frac{\Lambda^2}{2}\sigma_\tau^2.
\tag{B.1}
$$

To solve for the growth rate, write $\Lambda=\Lambda_0+\delta\Lambda$, where $\Lambda_0=\ln2/\mu$ is the noise-free rate. There is no first-order noise term in (B.1), so $\delta\Lambda$ also starts at second order. Substituting and retaining terms through this order gives

$$
\begin{aligned}
0
&\simeq\ln2-\mu(\Lambda_0+\delta\Lambda)
       +\frac{(\Lambda_0+\delta\Lambda)^2}{2}\sigma_\tau^2\\
&\simeq-\mu\,\delta\Lambda+\frac{\Lambda_0^2}{2}\sigma_\tau^2.
\end{aligned}
$$

The product $\delta\Lambda\,\sigma_\tau^2$ is fourth order and is discarded, which allows $\Lambda^2$ to be replaced by $\Lambda_0^2$ in the variance term. Solving gives $\delta\Lambda\simeq\Lambda_0^2\sigma_\tau^2/(2\mu)$, hence

$$
\Lambda\simeq
\frac{\ln2}{\mu}
+\frac{(\ln2)^2}{2\mu^3}\sigma_\tau^2.
$$

In the independent model, the leading correction is controlled by the lifetime variance. The correlated calculation below replaces the single-generation variance by the long-run variance of accumulated time.

#### B.2 Covariance bookkeeping

Assume the $\tau_n$ form a stationary process. Write

$$
C_\tau(k)
=\operatorname{Cov}(\tau_n,\tau_{n+k}),
\qquad C_\tau(-k)=C_\tau(k).
$$

Since $T_n=\sum_{j=0}^{n-1}\tau_j$,

$$
\operatorname{Var}(T_n)
=\sum_{i,j=0}^{n-1}
\operatorname{Cov}(\tau_i,\tau_j).
$$

Viewing this sum as a covariance matrix makes the bookkeeping explicit: diagonal entries contribute the individual variances, and off-diagonal entries occur in symmetric pairs.

![A covariance matrix highlights mirrored off-diagonal contributions, while accumulated-time variance grows at different rates for alternating, independent, and persistent correlations.](figures/07_accumulated_time_covariances.svg)

*Fig 8. Counting the fluctuations that accumulate.* Each square in panel A contributes to $\operatorname{Var}(T_n)$. The circled squares are the same covariance counted as $(i,j)$ and $(j,i)$, explaining the factor $2$ below. Panel B evaluates (4.2) for $C_\tau(k)=\sigma_\tau^2a^{|k|}$ with $a=-1/2,0,1/2$; panel A shows the $a=-1/2$ case. All three processes have the same one-generation variance, but their long-time slopes $D_\tau/\sigma_\tau^2=(1+a)/(1-a)$ differ. These are analytic stationary model covariances.

There are $n$ terms of separation zero, and $2(n-k)$ ordered pairs separated by $k$. Hence the exact finite-$n$ formula is

$$
\operatorname{Var}(T_n)
=nC_\tau(0)
+2\sum_{k=1}^{n-1}(n-k)C_\tau(k).
$$

If correlations are absolutely summable, division by $n$ and passage to the limit give

$$
\boxed{
D_\tau
\equiv\lim_{n\to\infty}\frac{\operatorname{Var}(T_n)}{n}
=C_\tau(0)+2\sum_{k=1}^{\infty}C_\tau(k).
}
$$

$D_\tau$ has units of time squared per generation. It measures *long-run fluctuations of accumulated division time*, not just fluctuations of a single generation. It is nonnegative because it is a limit of variances divided by $n$.

For independent division times, all nonzero-lag covariances vanish and $D_\tau=\sigma_\tau^2$. Positive correlations increase $D_\tau$; negative correlations can reduce it, even to zero.

<div id="42-cumulant-expansion-of-the-correlated-eulerlotka-equation"></div>

#### B.3 Asymptotic cumulants

One sufficient framework for the calculation is a family $\tau_j=\mu+\epsilon u_j$ with positive lifetimes, fixed correlation scales, and a scaled cumulant-generating function admitting a uniform expansion in $\epsilon$ near $s=\ln2/\mu$. Its cumulant rates must exist, with $K_m=O(\epsilon^m)$ and the higher-order remainder controlled there. Differentiating a large-$n$ limit or exchanging the cumulant series with that limit requires this regularity; covariance summability alone does not justify it. The approximation can fail when the correlation time diverges as noise is reduced, even if each lifetime has a small variance.

For small enough $s$ and a process with well-behaved asymptotic cumulants,

$$
\begin{aligned}
\frac1n\ln\mathbb E[e^{-sT_n}]
&=-s\frac{\mathbb E[T_n]}{n}
 +\frac{s^2}{2}\frac{\operatorname{Var}(T_n)}{n}\\
&\quad-\frac{s^3}{6}\frac{\kappa_3(T_n)}{n}+\cdots.
\end{aligned}
$$

Taking $n\to\infty$,

$$
\psi(s)=-s\mu+\frac{s^2}{2}D_\tau
-\frac{s^3}{6}K_3+\cdots,
\tag{B.2}
$$

where $K_3=\lim_{n\to\infty}\kappa_3(T_n)/n$, if the limit exists.

Substitute $s=\Lambda$ into (C.6):

$$
0=\ln2-\Lambda\mu
 +\frac{\Lambda^2}{2}D_\tau+O(\text{third cumulants}).
\tag{B.3}
$$

Using $\Lambda_0=\ln2/\mu$ as the zero-noise solution yields

$$
\boxed{
\Lambda\simeq\frac{\ln2}{\mu}
+\frac{(\ln2)^2}{2\mu^3}D_\tau.
}
$$

The approximation assumes sufficiently weak fluctuations and a regular expansion; it may fail quantitatively for large persistence times or broad lifetime distributions.

**Accumulated noise and population growth:**

$$
\text{Population growth responds to accumulated lineage noise,}
\;\text{not merely one-generation noise.}
$$

It is a *long-time* and *long-wavelength* property of the chronological process.


### Appendix C. Many-to-one identity and correlated Euler–Lotka {#appendix-c}

<div id="31-from-one-lifetime-to-an-entire-lineage"></div>

#### C.1 From one lifetime to an entire lineage

Consider a chronological lineage with correlated division times

$$
\tau_0,\tau_1,\ldots,\tau_{n-1}.
$$

Its $n$-th division/birth generation is reached at time

$$
T_n=\sum_{j=0}^{n-1}\tau_j.\tag{C.1}
$$

A one-generation histogram records how often each lifetime occurs, but not which lifetimes follow each other. For example, let a short lifetime $\tau_s$ and a long lifetime $\tau_\ell$ each occur with probability $1/2$. The following comparison shows how their ordering changes the distribution of their sum.

![Independent and alternating two-cycle histories have identical single-cycle lifetime probabilities but different distributions of total elapsed time.](figures/10_histories_and_total_time.svg)

*Fig 9. A lifetime distribution does not specify a lineage history.* Both examples assign probability $1/2$ to each of the two lifetimes. Panel A lists the probabilities of their ordered pairs. Independent sampling allows all four pairs; strict alternation allows only the two mixed pairs. Consequently, panel B has three possible total durations for independent sampling, while the alternating pair always lasts $\tau_s+\tau_\ell$. The correlated calculation must therefore use the distribution of $T_n$, which depends on the joint history.

After $n$ binary divisions there are $2^n$ possible descendant paths. For a locally defined branching process in which the statistical law along a uniformly sampled daughter path is the chronological-lineage law, linearity of expectation gives

$$
\mathbb E\!\left[\sum_{\text{generation-}n\text{ descendants}}
   e^{-s\,T_{n,\mathrm{path}}}\right]
=2^n\,\mathbb E_{\mathrm{chron}}[e^{-sT_n}].
\tag{C.2}
$$

This is the elementary *many-to-one* identity. It does not require all descendant paths to be statistically independent: it requires the correct marginal law for a randomly chosen daughter path.

Define the expected discounted number of generation-$n$ births,

$$
Z_n(s)=2^n\mathbb E_{\mathrm{chron}}[e^{-sT_n}].
\tag{C.3}
$$

The expected Laplace transform of all birth events is formally

$$
\widehat B(s)=\sum_{n=0}^{\infty}Z_n(s).
\tag{C.4}
$$

The exponential growth rate is the boundary between convergence and divergence of this series. If the correlated process admits an appropriate large-$n$ limit and the usual branching-process regularity conditions, define

$$
\psi(s)=\lim_{n\to\infty}
\frac1n\ln\mathbb E[e^{-sT_n}].
\tag{C.5}
$$

Then $Z_n(s)$ has exponential rate $\ln2+\psi(s)$. Its critical value is determined by

$$
\boxed{\ln2+\psi(\Lambda)=0.}\tag{C.6}
$$

Equivalently,

$$
\ln2+
\lim_{n\to\infty}\frac1n
\ln\mathbb E\!\left[
\exp\left(-\Lambda\sum_{j=0}^{n-1}\tau_j\right)
\right]=0.
\tag{C.7}
$$

This is the correlated-lineage generalization we will use. Its derivation is particularly transparent, and mathematically exact, for the finite-state Markov model of Appendix D.

**Assumption check.** Equation (C.7) is an asymptotic relation, not a universal identity for every imaginable correlated or environmentally forced population. We assume a stationary/ergodic lineage law, an appropriate scaled cumulant-generating-function limit, and a branching population in which a well-defined Malthusian exponent is selected by this discounted-birth condition. Strong common-environment noise, nonergodic mixtures, or unbounded-growth pathologies require separate care.

For clarity, the convergence argument can be made explicit. Assume strictly positive generation times, no explosion, and finite expected birth counts on every bounded time interval. Let $B_{\mathrm{cum}}(t)$ count all births up to time $t$, including the founder. In a binary population without death, each division creates two births and increases the living count by one, so pathwise $N(t)=[B_{\mathrm{cum}}(t)+1]/2$. Tonelli's theorem for nonnegative weights gives, for real $s>0$,

$$
\int_0^\infty e^{-st}M(t)\,dt
=\frac{1+\widehat B(s)}{2s}.\tag{C.8}
$$

Thus expected living counts and the expected discounted-birth series have the same positive convergence boundary. If $r(s)=\ln2+\psi(s)<0$, the root test makes $\sum_n Z_n(s)$ converge; if $r(s)>0$, its terms do not even tend to zero. At $r(s)=0$ the root test is inconclusive about convergence exactly at that point. If $r$ is continuous and crosses zero at a unique positive $s_*$, the boundary is nevertheless $s_*$. When $\lim_{t\to\infty}t^{-1}\ln M(t)$ exists, comparison with exponentials identifies this boundary with that limit, proving (C.6) for **expected** growth. Existence of $\psi$ without continuity, a crossing, or an expected-growth limit does not complete this argument.

To identify the same exponent in a typical realized population requires a branching theorem for the specified model: appropriate irreducibility and integrability, and reproduction governed locally rather than by a shared random environment. For example, finite-state irreducible age-dependent branching with lifetimes bounded above and away from zero, two offspring, and conditionally independent daughter subtrees supplies a standard sufficient setting. Sibling independence is not required for the many-to-one expectation itself; a typical-growth assertion with dependent reproduction needs its own hypotheses. Nonlattice lifetimes are needed for some constant-amplitude limits, not for every logarithmic growth exponent. The stationary chronological law is the input to $\psi$; a different founder law must have the same exponential asymptotics for the same root to apply.

<div id="32-independent-times-reappear-as-a-special-case"></div>

#### C.2 Independent times reappear as a special case

If $\tau_0,\ldots,\tau_{n-1}$ are independent,

$$
\mathbb E[e^{-sT_n}]
=\prod_{j=0}^{n-1}\mathbb E[e^{-s\tau_j}]
=\big(\mathbb E[e^{-s\tau}]\big)^n.
$$

Thus

$$
\psi(s)=\ln\mathbb E[e^{-s\tau}],
$$

and (C.6) reduces immediately to

$$
2\mathbb E[e^{-\Lambda\tau}]=1.
$$

The original Euler–Lotka equation has **not been discarded**; it is the uncorrelated limit of the new calculation.

<div id="33-the-mean-time-inequality-still-holds"></div>

#### C.3 The mean-time inequality still holds

For stationary generation times with mean $\mu=\mathbb E[\tau_n]$, Jensen's inequality applied to the *entire sum* gives

$$
\begin{aligned}
\frac1n\ln\mathbb E[e^{-sT_n}]
&\ge-\frac{s}{n}\mathbb E[T_n]\\
&=-s\mu.
\end{aligned}
$$

Therefore $\psi(s)\ge-s\mu$. At $s=\Lambda$,

$$
-\ln2=\psi(\Lambda)\ge-\Lambda\mu,
$$

so

$$
\Lambda\ge\frac{\ln2}{\mu}.\tag{C.9}
$$

Unlike the independent case, equality can hold **even when each $\tau_n$ has nonzero variance**. The key is whether fluctuations of $T_n$ accumulate extensively as $n$ grows.


<div id="5-exactly-solvable-example-i-a-two-state-markov-model"></div>

### Appendix D. A complete two-state Markov branching model {#appendix-d}

This model isolates the effect of temporal correlation while **holding the one-generation distribution fixed**. It is a phenomenological branching model, not automatically a physically consistent constant-$\lambda$, size-homeostatic cell model.

<div id="51-definition"></div>

#### D.1 Definition

There are two possible generation times:

$$
\tau_s=\mu-\delta,
\qquad
\tau_\ell=\mu+\delta,
\qquad 0\lt \delta\lt \mu.
$$

The cell has a state $s$ (short) or $\ell$ (long). From one generation to the next, the state is retained with probability $p$:

$$
P=
\begin{pmatrix}
p&1-p\\
1-p&p
\end{pmatrix}.
\tag{D.1}
$$

A cell in state $i$ lives for a time $\tau_i$ and then produces two daughters. Each daughter's state follows the transition matrix $P$, conditional on the mother's state. For $0\lt p\lt 1$, the chronological stationary state probabilities are $(1/2,1/2)$, independent of $p$. Therefore

$$
\mathbb P(\tau=\tau_s)=\mathbb P(\tau=\tau_\ell)=\tfrac12,
\quad
\mathbb E[\tau]=\mu,
\quad
\operatorname{Var}(\tau)=\delta^2.
\tag{D.2}
$$

Changing $p$ changes **only the temporal organization** of the division times.

![Two-state generation-time autocorrelations and exact population growth as state persistence increases.](figures/02_correlations_and_growth.svg)

*Fig 10. Correlation changes growth without changing the lifetime histogram.* Panel A plots $\rho_\tau(k)=(2p-1)^k$. Panel B solves the exact Perron-eigenvalue condition (D.7) for the same 2 h and 5 h lifetimes used in the table; dots mark its four rows. The dashed line is $\ln2/\mu$. The curve ends below $p=1$, where the chain becomes reducible. These are predictions of the two-state branching model, not measurements of size-homeostatic cells.

<div id="52-exact-autocorrelation"></div>

#### D.2 Exact autocorrelation

Let $z_n=-1$ in state $s$ and $z_n=+1$ in state $\ell$. Then

$$
\tau_n=\mu+\delta z_n.
$$

Conditional on the current state,

$$
\mathbb E[z_{n+1}\mid z_n]
=p z_n+(1-p)(-z_n)
=(2p-1)z_n.
$$

Set

$$
a=2p-1.
$$

By iterating conditional expectations,

$$
\mathbb E[z_{n+k}\mid z_n]=a^k z_n.
$$

Since $\mathbb E[z_n]=0$ and $z_n^2=1$,

$$
C_\tau(k)=\delta^2 a^{|k|},\qquad a=2p-1.\tag{D.3}
$$

For $0\lt p\lt 1$, $|a|\lt 1$. Equation (4.3) then gives

$$
\begin{aligned}
D_\tau
&=\delta^2\left(1+2\sum_{k=1}^\infty a^k\right)\\
&=\delta^2\frac{1+a}{1-a}\\
&=\delta^2\frac{p}{1-p}.
\end{aligned}
\tag{D.4}
$$

Three regimes follow:

- $p\lt 1/2$: successive generation times tend to alternate; correlations are negative at odd lags.
- $p=1/2$: generation times are independent; $D_\tau=\delta^2$.
- $p>1/2$: states persist; $D_\tau>\delta^2$.

At exactly $p=0$, the states alternate deterministically; the absolute covariance sum in (4.3) does not converge, but directly from $T_n$ one finds $\operatorname{Var}(T_n)=O(1)$ and $D_\tau=0$. At $p=1$, the chain is nonergodic; its two disjoint subpopulations require separate treatment.

<div id="53-exact-population-growth-via-a-transfer-matrix"></div>

#### D.3 Exact population growth via a transfer matrix

The joint probability of a state path $i_0,i_1,\ldots,i_n$ contains a product of transition probabilities. Its discounted birth-time weight contains

$$
\exp\left(-\Lambda\sum_{j=0}^{n-1}\tau_{i_j}\right).
$$

At each generation, binary division contributes a factor of 2. Inheritance supplies the transition probability, and the mother's lifetime supplies the exponential discount, as shown below.

![A two-state directed graph assigns each offspring transition the weight two times its inheritance probability times the discount for the mother's lifetime.](figures/08_discounted_branching.svg)

*Fig 11. Reading the transfer matrix as a branching graph.* Each outgoing arrow from state $i$ carries the expected daughter count $2P_{ij}$ multiplied by $e^{-s\tau_i}$. Multiplying these weights along a path gives its contribution to $Z_n(s)$. The diagram uses a trial rate $s$; at the selected growth rate $s=\Lambda$, the leading eigenvalue is $1$, as in (D.7).

Therefore the propagation matrix for expected discounted descendant weights is

$$
K_\Lambda=
2\begin{pmatrix}
e^{-\Lambda\tau_s}&0\\
0&e^{-\Lambda\tau_\ell}
\end{pmatrix}P.
\tag{D.5}
$$

Define

$$
A=e^{-\Lambda\tau_s},\qquad B=e^{-\Lambda\tau_\ell}.
$$

Then

$$
K_\Lambda=
2\begin{pmatrix}
pA&(1-p)A\\
(1-p)B&pB
\end{pmatrix}.
\tag{D.6}
$$

For an initial-state probability row vector $\boldsymbol\pi$, the generation-$n$ discounted weight is

$$
Z_n(\Lambda)=\boldsymbol\pi K_\Lambda^n\mathbf 1.
$$

For $0\lt p\lt 1$, Perron–Frobenius theory gives its leading growth factor $\rho(K_\Lambda)^n$, where $\rho$ is the spectral radius. Equation (C.6) is therefore equivalent to the **exact eigenvalue condition**

$$
\rho(K_\Lambda)=1.\tag{D.7}
$$

To compute it explicitly, the trace and determinant are

$$
\operatorname{tr}K_\Lambda=2p(A+B),
\qquad
\det K_\Lambda=4AB(2p-1).
$$

The Perron eigenvalue is

$$
\rho(K_\Lambda)
=p(A+B)+\sqrt{p^2(A-B)^2+4(1-p)^2AB}.
\tag{D.8}
$$

Solving $\rho(K_\Lambda)=1$ for $\Lambda$ yields the population growth rate. This calculation is exact for the specified model; it is not a cumulant approximation.

<div id="54-two-analytically-transparent-limits"></div>

#### D.4 Two analytically transparent limits

**Case 1: $p=1/2$, independent generations.** Then

$$
\rho(K_\Lambda)=A+B,
$$

so (D.7) becomes

$$
e^{-\Lambda\tau_s}+e^{-\Lambda\tau_\ell}=1.
\tag{D.9}
$$

This is precisely $2\mathbb E[e^{-\Lambda\tau}]=1$ for the two-point distribution.

**Case 2: $p=0$, strict alternation.** Then

$$
\rho(K_\Lambda)=2\sqrt{AB}
=2e^{-\Lambda(\tau_s+\tau_\ell)/2}.
$$

Therefore

$$
\Lambda=\frac{\ln2}{(\tau_s+\tau_\ell)/2}
      =\frac{\ln2}{\mu}.\tag{D.10}
$$

The single-generation distribution still has variance $\delta^2>0$, yet the growth advantage disappears. Every pair of successive generation times has exactly the same total duration.

For $p\to1^-$, long-lived state persistence lets fast lineages dominate the population; the growth exponent approaches $\ln2/\tau_s$. The point $p=1$ itself is reducible, so the observed exponent depends on which state classes are initially present.

<div id="55-numerical-example-same-marginal-distribution-different-growth"></div>

#### D.5 Numerical example: same marginal distribution, different growth

Choose

$$
\tau_s=2\ \mathrm h,
\qquad
\tau_\ell=5\ \mathrm h,
\qquad
\mu=C.5\ \mathrm h.
$$

Solving (D.8)–(D.7):

| Persistence $p$ | Correlation structure | $\Lambda$ ($\mathrm h^{-1}$) | Doubling time $\ln2/\Lambda$ (h) |
|---:|---|---:|---:|
| 0 | Strict alternation | 0.198042 | 3.500 |
| 0.25 | Antipersistent | 0.202456 | 3.424 |
| 0.50 | Independent | 0.212289 | 3.265 |
| 0.75 | Persistent | 0.244769 | 2.832 |

Every row has **the same** $p(\tau)$, $\bar\tau$, and $\operatorname{Var}(\tau)$. Only correlations differ. This branching model isolates the statistical role of correlations; it does not require all cases to satisfy both fixed size-growth rate and size homeostasis. The rows should therefore not be interpreted directly as cellular controls under identical physiological conditions.

> **Lifetime distributions and correlations.** A one-generation histogram is not enough to predict the proliferation rate. The same random variables, arranged differently along a lineage, can produce different collective growth rates.


### Appendix E. Adder variance, size memory, and covariance cancellation {#appendix-e}

For the i.i.d. added-size model in Section 3, center $v_n$ at $\Delta$ and write

$$
x_{n+1}=a x_n+\xi_n,\qquad a=1/2,\quad\xi_n=\eta_n/2.\tag{E.1}
$$

![Small and large cells add the same size increment before division, and their daughters move closer to the mean birth size.](figures/03_adder_size_control.svg)

*Fig 12. The adder as restoring feedback.* In this noise-free illustration, cells born at $0.70\Delta$ and $1.30\Delta$ both add $\Delta$ and divide symmetrically. Their daughters are born at $0.85\Delta$ and $1.15\Delta$: the deviation from the mean is halved. Capsule lengths encode relative size schematically; they are not microscopy images. Added-size noise produces the stochastic map (E.1).

<div id="62-exact-stationary-variance"></div>

#### E.1 Exact stationary variance

In the stationary state, $x_n$ and $\xi_n$ are independent because $\xi_n$ is newly sampled noise. Thus

$$
\operatorname{Var}(x_{n+1})
=a^2\operatorname{Var}(x_n)+\operatorname{Var}(\xi_n).
$$

Let $\sigma_x^2=\operatorname{Var}(x_n)$. Then

$$
\sigma_x^2=a^2\sigma_x^2+\sigma_\xi^2,
$$

or

$$
\sigma_x^2=\frac{\sigma_\xi^2}{1-a^2}.\tag{E.2}
$$

For $a=1/2$ and $\sigma_\xi^2=\sigma_\eta^2/4$,

$$
\operatorname{Var}(v_n)=\frac{\sigma_\eta^2}{3}.\tag{E.3}
$$

Unlike the timer, the adder has a finite stationary size variance. The deterministic part of the map shrinks deviations by a factor of 2 every generation.

<div id="63-exact-birth-size-autocorrelation"></div>

#### E.2 Exact birth-size autocorrelation

Iterating (E.1) $k$ times,

$$
x_{n+k}
=a^k x_n+
\sum_{j=0}^{k-1}a^{k-1-j}\xi_{n+j}.
\tag{E.4}
$$

The future noises are independent of $x_n$, so

$$
\operatorname{Cov}(x_n,x_{n+k})
=a^k\sigma_x^2.
$$

Dividing by $\sigma_x^2$ gives

$$
\rho_v(k)
\equiv\operatorname{Corr}(v_n,v_{n+k})
=a^k=2^{-k}.
\tag{E.5}
$$

For example,

$$
\rho_v(1)=\frac12,
\qquad
\rho_v(2)=\frac14.
$$

These values agree closely with the correlations $\approx0.49$ and $\approx0.24$ listed on p. 5 of the [Ho laboratory materials](https://ho.lab.westlake.edu.cn/teaching_en.html), which cite Ho et al., *Annual Review of Biophysics* (2018).

The important point is not merely the numerical agreement: a deterministic restoring factor fixes the decay of the correlation function. Noise sets the variance, while feedback sets the memory time.

![Adder predictions and experimental relaxation of initially small and large birth sizes.](figures/taheri_araghi_2015_f3_highres.png)

*Fig 13. Observed return toward the stationary birth size.* Original Figure 3 from [Taheri-Araghi et al. (2015)](https://doi.org/10.1016/j.cub.2014.12.009). Panel A predicts convergence under constant added size; panels B and C compare with *E. coli* and *B. subtilis*. This illustrates the restoring factor in (E.1). Conditional mean relaxation and the stationary ACF (E.5) are related observables, but are distinct measurements.

#### E.3 Generation-time covariances in the small-fluctuation limit

In the adder, nearby birth sizes satisfy

$$
\rho_v(k)=2^{-k}>0.
$$

What about division times?

For small size fluctuations $|x_n|\ll\Delta$, expand the logarithm in (3.4):

$$
\ln v_n=\ln(\Delta+x_n)
\simeq\ln\Delta+\frac{x_n}{\Delta}.
$$

Therefore

$$
\delta\tau_n
\equiv\tau_n-\frac{\ln2}{\lambda}
\simeq\frac{x_{n+1}-x_n}{\lambda\Delta}.
\tag{E.6}
$$

For a general stationary AR(1) process

$$
C_x(k)=\sigma_x^2 a^{|k|},\qquad 0\lt a\lt 1,
$$

we can compute the covariance of the differences explicitly. At lag zero,

$$
\begin{aligned}
C_{\delta\tau}(0)
&\simeq\frac1{(\lambda\Delta)^2}
\mathbb E[(x_{n+1}-x_n)^2]\\
&=\frac{2\sigma_x^2(1-a)}{(\lambda\Delta)^2}.
\end{aligned}
\tag{E.7}
$$

For $k\ge1$,

$$
\begin{aligned}
C_{\delta\tau}(k)
&\simeq\frac1{(\lambda\Delta)^2}
\operatorname{Cov}(x_{n+1}-x_n,\,x_{n+k+1}-x_{n+k})\\
&=\frac{2C_x(k)-C_x(k-1)-C_x(k+1)}{(\lambda\Delta)^2}\\
&=-\frac{\sigma_x^2(1-a)^2 a^{k-1}}{(\lambda\Delta)^2}.
\end{aligned}
\tag{E.8}
$$

Dividing by (E.7) yields

$$
\rho_\tau(k)
\simeq-\frac{1-a}{2}a^{k-1},
\qquad k\ge1.
\tag{E.9}
$$

For the adder, $a=1/2$:

$$
\rho_\tau(k)\simeq-2^{-(k+1)},\qquad k\ge1.\tag{E.10}
$$

Thus

$$
\rho_\tau(1)\simeq-\tfrac14,
\qquad
\rho_\tau(2)\simeq-\tfrac18,
\qquad
\rho_\tau(3)\simeq-\tfrac1{16}.
$$

**Two distinct autocorrelation functions:** birth sizes are positively correlated across generations, while successive *generation-time fluctuations* are negatively correlated in this model. Confusing these two autocorrelation functions would obscure the mechanism.

The negative correlation has a simple interpretation: a cell with a longer-than-usual division time tends to produce a larger daughter, and that daughter typically requires less time to divide next.

<div id="68-an-explicit-check-of-the-cancellation"></div>

#### E.4 An explicit check of the cancellation

Using (E.7)–(E.8),

$$
\begin{aligned}
D_\tau
&\simeq C_{\delta\tau}(0)+2\sum_{k=1}^\infty C_{\delta\tau}(k)\\
&=\frac{2\sigma_x^2(1-a)}{(\lambda\Delta)^2}
-\frac{2\sigma_x^2(1-a)^2}{(\lambda\Delta)^2}
\sum_{k=1}^\infty a^{k-1}\\
&=\frac{2\sigma_x^2(1-a)}{(\lambda\Delta)^2}
-\frac{2\sigma_x^2(1-a)}{(\lambda\Delta)^2}\\
&=0.
\end{aligned}
\tag{E.11}
$$

The sum of the negative nonzero-lag covariances cancels the positive single-generation variance. Equation (3.5) showed that the cancellation is actually exact at the level of cumulative times, even when the small-noise approximation used in (E.6) is not justified.

> **The dynamical constraint from size homeostasis.** Size homeostasis is not just a narrow size distribution. It imposes a dynamical constraint on the entire lineage history. Under constant exponential size growth and symmetric division, this constraint removes the extensive component of accumulated generation-time fluctuations.


<div id="appendix-b-general-linear-size-control-laws-timer-adder-and-sizer"></div>

### Appendix F. General linear size control and complete spectra {#appendix-f}

The strategy comparison is introduced in Section 1.3; here we retain the unified linear description and the full frequency-domain calculations.

#### F.1 Local linear feedback laws

A unified local description is

$$
x_{n+1}=a x_n+\xi_n,
\tag{F.1}
$$

where $x_n$ is deviation of birth size from its stationary mean.

| Strategy | Typical size-map slope $a$ | Physical behavior |
|---|---:|---|
| Sizer | $0$ | Division size independent of birth size; no inherited linear size memory |
| Adder | $1/2$ | Half the birth-size deviation is inherited |
| Timer, linearized near a balanced size | $1$ | No restoring feedback; random-walk/neutral behavior |

For $|a|\lt 1$,

$$
\operatorname{Var}(x)=\frac{\operatorname{Var}(\xi)}{1-a^2},
\qquad
\rho_x(k)=a^{|k|}.
$$

These familiar control strategies differ in **how quickly birth-size perturbations are forgotten**. The observed slope alone constrains the effective dynamics, but not a unique biochemical implementation.

The timer row is a local or additive approximation: the exact multiplicative timer map is the log-size random walk derived in Section 1.

The following calculations give the full spectra behind Section 4.3, using frequency per generation.

<div id="71-autocovariance-and-the-wienerkhinchin-relation"></div>

#### F.2 Autocovariance and the Wiener–Khinchin relation

For a stationary discrete-time process $y_n$, define

$$
C_y(k)=\mathbb E[(y_n-\bar y)(y_{n+k}-\bar y)].
$$

Its two-sided power spectral density is

$$
S_y(\omega)=\sum_{k=-\infty}^{\infty}
C_y(k)e^{-i\omega k},
\qquad -\pi\le\omega\le\pi.
\tag{F.2}
$$

The inverse relation is

$$
C_y(k)=\frac1{2\pi}\int_{-\pi}^{\pi}
S_y(\omega)e^{i\omega k}\,d\omega.
\tag{F.3}
$$

These are the discrete-time Wiener–Khinchin relations. The frequency $\omega$ is angular frequency **per generation**. This distinction matters when comparing lineage-indexed observations with measurements sampled in physical time.

At zero frequency,

$$
S_\tau(0)=\sum_{k=-\infty}^{\infty}C_\tau(k)
=C_\tau(0)+2\sum_{k=1}^\infty C_\tau(k).
$$

Thus, whenever the covariance sum converges,

$$
S_\tau(0)=D_\tau.
$$

This gives equation (4.4) a spectral form:

$$
\Lambda\simeq
\frac{\ln2}{\mu}
+\frac{(\ln2)^2}{2\mu^3}S_\tau(0).
$$

**Growth corrections in the frequency domain:** at quadratic order, long-time population growth is controlled by *zero-frequency generation-time noise*. The full power spectrum additionally tells us whether there are preferred dynamical timescales.

<div id="72-spectrum-of-the-ar1-size-control-model"></div>

#### F.3 Spectrum of the AR(1) size-control model

Consider again

$$
x_{n+1}=a x_n+\xi_n,
\qquad |a|\lt 1,
$$

with independent zero-mean noise of variance $\sigma_\xi^2$. In the stationary regime,

$$
x_n=\sum_{j=0}^\infty a^j\xi_{n-1-j}.
\tag{F.4}
$$

This is a linear filter of white noise. In the frequency domain the transfer function (including the one-generation delay) is

$$
H(\omega)=\frac{e^{-i\omega}}{1-ae^{-i\omega}},
$$

and therefore

$$
\begin{aligned}
S_x(\omega)
&=\sigma_\xi^2|H(\omega)|^2\\
&=\frac{\sigma_\xi^2}
{1+a^2-2a\cos\omega}.
\end{aligned}
\tag{F.5}
$$

For an adder, $a=1/2$. The size spectrum is greatest near $\omega=0$ and decreases smoothly toward $\pi$. There is no intrinsic *finite-frequency* oscillatory resonance.

This is the spectral counterpart of the monotone, positive autocorrelation

$$
C_x(k)=\sigma_x^2(1/2)^{|k|}.
$$

<div id="73-spectrum-of-generation-time-feedback-filters-out-zero-frequency"></div>

#### F.4 Spectrum of generation time: feedback filters out zero frequency

From (E.6),

$$
\delta\tau_n\simeq\frac{x_{n+1}-x_n}{\lambda\Delta}.
$$

Taking a difference multiplies a Fourier mode $e^{i\omega n}$ by $e^{i\omega}-1$. Hence

$$
\begin{aligned}
S_\tau(\omega)
&\simeq\frac{|e^{i\omega}-1|^2}{(\lambda\Delta)^2}S_x(\omega)\\
&=\boxed{
\frac{4\sin^2(\omega/2)}{(\lambda\Delta)^2}S_x(\omega).
}
\end{aligned}
\tag{F.6}
$$

At zero frequency the numerator vanishes:

$$
S_\tau(0)=0.\tag{F.7}
$$

The same fact appeared as the telescoping sum in (3.5) and the covariance cancellation in (E.11). These are three mathematical descriptions of **one physical mechanism**.

Notice that differencing may emphasize high-frequency fluctuations, even when the underlying size process does not have an oscillatory mode. A large spectral value near the Nyquist frequency $\omega=\pi$ is not, by itself, evidence for a biological oscillator.

<div id="74-an-extra-connection-the-markov-models-spectrum"></div>

#### F.5 An extra connection: the Markov model's spectrum

For the two-state Markov model of Appendix D, $C_\tau(k)=\delta^2a^{|k|}$, where $a=2p-1$. Summing the Fourier series for $|a|\lt 1$ gives

$$
S_\tau(\omega)
=\delta^2\frac{1-a^2}{1+a^2-2a\cos\omega}.
\tag{F.8}
$$

In particular,

$$
S_\tau(0)=\delta^2\frac{1+a}{1-a}
=\delta^2\frac{p}{1-p},
$$

which recovers (D.4) without separately summing the covariances.

Thus both models are readily understood through spectra:

- The Markov model can change low-frequency division-time noise by changing persistence.
- The adder converts size fluctuations into differences, suppressing division-time noise at zero frequency.


<div id="appendix-c-why-the-division-time-cancellation-is-more-general-than-the-adder"></div>

### Appendix G. Scope of exact cancellation and sufficient asymptotic conditions {#appendix-g}

#### G.1 Identity, moments, and growth are different levels of assertion

The identity (3.5) holds at every finite $n$ whenever sizes are positive, growth has one constant exponential rate, and divisions are symmetric. It is valid even for a nonstationary size process; homeostasis is needed to bound its endpoint term, not to make it telescope. In a stationary process, finite $\mathbb E|\ln v|$ fixes the mean time, and finite $\operatorname{Var}(\ln v)$ bounds $\operatorname{Var}(T_n)$ and gives $D_\tau=0$. Mixing strengthens this bound to the variance limit in Section 3.3. None of these moment conditions alone controls the exponential tilt in (3.7).

For a stationary positive-size process, put $q=s/\lambda>0$. A useful sufficient condition for that tilt is

$$
1\le\mathbb E[(v_n/v_0)^{-q}]
\le\big(\mathbb E[v^{2q}]\mathbb E[v^{-2q}]\big)^{1/2}<\infty.\tag{G.1}
$$

The lower bound follows from Jensen and $\mathbb E[\ln v_n-\ln v_0]=0$; Cauchy–Schwarz supplies the upper bound without requiring endpoint independence. These moments imply integrability of the logarithm, and the bounds are uniform in $n$. They therefore give $\lim_n n^{-1}\ln\mathbb E[(v_n/v_0)^{-q}]=0$. If they hold for $s$ in a neighborhood of $\lambda$, the expected discounted-birth weights have exponential rate $\ln2-s\ln2/\lambda$ there, selecting $\lambda$ under Appendix C's conditions. Fixed initial sizes can be treated with corresponding uniform endpoint bounds, but are not covered just by assuming stationary marginals.

An adder with added increments bounded above and bounded away from zero provides a simple stationary realization with all these moments. Merely requiring positive increments and finite added-size variance does not automatically establish the inverse moments or the population snapshot condition. The alternative conservation proof uses $t^{-1}\ln\overline v_{\mathrm{snap}}(t)\to0$ and applies directly to a realized population. Chronological stationarity and snapshot homeostasis should not be silently interchanged.

#### G.2 Universality within the stated kinematics

Equation (3.4) used only

$$
v_n^{\mathrm d}=v_n e^{\lambda\tau_n},
\qquad
v_{n+1}=\frac12v_n^{\mathrm d}.
$$

The functional form of the division-size control law did not enter. Therefore any stationary positive-size process satisfying those two relations obeys (3.5). A sizer, adder, or another feedback law can all exhibit this cancellation if they provide sufficient homeostasis.

For an asymmetric division where the followed daughter receives fraction $r_n$ of the mother's size,

$$
v_{n+1}=r_n v_ne^{\lambda\tau_n},
$$

and hence

$$
\tau_n=\frac{-\ln r_n+\ln v_{n+1}-\ln v_n}{\lambda}.
$$

The accumulated time becomes

$$
T_n=\frac1\lambda\sum_{j=0}^{n-1}(-\ln r_j)
+\frac{\ln v_n-\ln v_0}{\lambda}.
\tag{G.2}
$$

The first term is no longer deterministic when partition fractions fluctuate. Thus the statement $D_\tau=0$ does *not* follow from homeostasis alone under random asymmetry. Nevertheless, population-level mass conservation still gives $\Lambda=\lambda$ when every division conserves total size and the population mean snapshot size remains stationary. The appropriate lineage growth calculation must account for how daughter partition fractions and sampling are weighted.

A similar caution applies when individual growth rates fluctuate: if

$$
\dot v=\lambda_n v
$$

during generation $n$, then

$$
\tau_n=\frac{\ln2+\ln v_{n+1}-\ln v_n}{\lambda_n},
$$

and the sum no longer telescopes with one common prefactor. This is a natural route to genuinely new growth-versus-control effects.


<div id="appendix-a-kl-divergence-for-lifetime-and-trajectory-distributions"></div>

### Appendix H. KL divergence for lifetime and trajectory distributions {#appendix-h}

The weighted density $q$ introduced in Section 2.1 gives an information-theoretic expression for the growth increment relative to the mean-generation-time baseline. For independent lifetimes with $\bar\tau=\mathbb E_p[\tau]\lt\infty$, the Kullback–Leibler (KL) divergence between $p$ and $q$ is

$$
\begin{aligned}
D_{\mathrm{KL}}(p\|q)
&=\int_0^\infty p(\tau)\ln\frac{p(\tau)}{q(\tau)}\,d\tau\\
&=\int_0^\infty p(\tau)(\Lambda\tau-\ln2)\,d\tau\\
&=\Lambda\bar\tau-\ln2.
\end{aligned}
\tag{H.1}
$$

Nonnegativity of KL divergence recovers the Jensen bound (2.4), while $D_{\mathrm{KL}}(p\|q)/\bar\tau$ gives the exact excess of $\Lambda$ over $\ln2/\bar\tau$. The same distributional bias therefore has both a sampling interpretation and a KL representation.

For correlated lifetimes, Section 4.1 shows that the growth rate depends on complete division histories. The KL identity therefore uses a joint distribution of lifetimes. Let $\boldsymbol\tau=(\tau_0,\ldots,\tau_{n-1})$, let $p_n(\boldsymbol\tau)$ denote its chronological joint law, and write $T_n=\sum_{j=0}^{n-1}\tau_j$. Define the weighted trajectory distribution

$$
q_n(\boldsymbol\tau)
=\frac{2^n e^{-\Lambda T_n}p_n(\boldsymbol\tau)}{Z_n(\Lambda)},
\quad
Z_n(\Lambda)=2^n\mathbb E_p[e^{-\Lambda T_n}].
\tag{H.2}
$$

This is normalized by construction. For independent lifetimes, $Z_n(\Lambda)=1$ and $q_n(\boldsymbol\tau)=\prod_{j=0}^{n-1}q(\tau_j)$, recovering the single-lifetime weighting in Section 2.1. For a general joint law, the KL divergence is

$$
\begin{aligned}
D_{\mathrm{KL}}(p_n\|q_n)
&=\mathbb E_p\!\left[
\ln\frac{p_n}{q_n}\right]\\
&=\Lambda\mathbb E_p[T_n]-n\ln2+\ln Z_n(\Lambda).
\end{aligned}
\tag{H.3}
$$

If $\mathbb E[T_n]=n\bar\tau$ and $\lim_{n\to\infty}\ln Z_n(\Lambda)/n=0$, dividing by $n$ gives

$$
\lim_{n\to\infty}\frac1nD_{\mathrm{KL}}(p_n\|q_n)
=\Lambda\bar\tau-\ln2.
\tag{H.4}
$$

The nonnegative KL divergence *rate* extends (H.1) to correlated histories and recovers the generalized Jensen bound in Section 4.1. For correlated processes, $q_n$ is an exponentially tilted trajectory ensemble; identifying it with a specific retrospective sampling protocol can involve boundary effects and additional assumptions.

In a size-homeostatic, constant-$\lambda$ model, the extensive KL rate vanishes under the same conditions that yield $\Lambda\bar\tau=\ln2$, even though finite-path distributions need not be identical.


---

<div id="8-ordinary-fluctuations-or-an-oscillatory-dynamical-mode"></div>

## Further Topic: Fluctuations versus Oscillatory Modes {#further-topic}

Zero-frequency noise cancellation and nonzero-frequency oscillatory modes address different questions. The former concerns whether fluctuations accumulate over many generations; the latter concerns whether the dynamics have a characteristic oscillatory timescale.

<div id="81-relaxation-is-not-oscillation"></div>

### FT.1 Relaxation is not oscillation

An AR(1) model with $0\lt a\lt 1$ has a single real positive eigenvalue $a$. Perturbations decay as $a^n$. Its normalized ACF is monotone:

$$
\rho_x(k)=a^{|k|}.
$$

The spectrum has no isolated nonzero-frequency resonant mode. Random sample paths may visibly rise and fall, but this alone does not imply a periodic dynamical mechanism.

If $-1\lt a\lt 0$, the ACF alternates sign, $\rho_x(k)=a^{|k|}$, and the spectral maximum is at $\omega=\pi$, corresponding to a two-generation alternation. This already illustrates why the *sign and eigenstructure of feedback* matter more than the visual regularity of a noisy trace.

<div id="82-a-minimal-model-capable-of-damped-oscillations"></div>

### FT.2 A minimal model capable of damped oscillations

Introduce a delayed feedback term:

$$
x_{n+1}=2r\cos\omega_0\,x_n-r^2x_{n-1}+\xi_n,
\qquad 0\lt r\lt 1,
\quad0\lt \omega_0\lt \pi.
\tag{FT.1}
$$

The deterministic recurrence has characteristic equation

$$
z^2-2r\cos\omega_0\,z+r^2=0.
$$

Its roots are

$$
z_\pm=re^{\pm i\omega_0}.\tag{FT.2}
$$

The modulus $r$ controls damping, while the argument $\omega_0$ fixes the characteristic angular frequency. The deterministic response therefore contains terms of the form

$$
x_n\sim r^n\cos(n\omega_0+\phi).
\tag{FT.3}
$$

White noise continually excites these damped modes, leading to a stationary noisy oscillator.

The stationary autocovariance obeys the same second-order recurrence at nonzero lags:

$$
C_x(k+1)=2r\cos\omega_0\,C_x(k)-r^2C_x(k-1),
\quad k\ge1.
\tag{FT.4}
$$

Hence its normalized ACF has the form

$$
\rho_x(k)=r^{|k|}
\big[\cos(k\omega_0)+b\sin(|k|\omega_0)\big],
$$

where the coefficient $b$ is fixed by stationarity. For this AR(2) model,

$$
b=\frac{1-r^2}{1+r^2}\cot\omega_0.
\tag{FT.5}
$$

The power spectrum is

$$
S_x(\omega)=\frac{\sigma_\xi^2}
{\left|1-2r\cos\omega_0\,e^{-i\omega}
+r^2e^{-2i\omega}\right|^2}.
\tag{FT.6}
$$

For sufficiently weak damping and a frequency away from the boundaries, this exhibits a resonant peak near $\omega_0$. The exact peak location need not equal $\omega_0$ at finite damping.

![Adder spectra suppress zero-frequency time noise, while a delayed-feedback oscillator has a finite-frequency peak and oscillating autocorrelation.](figures/05_spectra_and_oscillations.svg)

*Fig 14. Zero-frequency cancellation and oscillatory modes answer different questions.* Panel A uses the adder's $a=1/2$ and the linearized generation-time relation (F.6); each spectrum is divided by its own process variance. Panel B compares the AR(1) size process with the AR(2) model (FT.1), using $r=0.9$ and $\omega_0=\pi/3$ radians per generation. The dotted line marks the mode angle; the spectral maximum can be shifted by damping. Panel C shows the corresponding analytic autocorrelations. All spectra use generation-indexed frequency, not physical-time frequency.

<div id="83-what-the-spectrum-can-and-cannot-establish"></div>

### FT.3 What the spectrum can and cannot establish

A finite-frequency peak suggests a preferred dynamical timescale, but does not, by itself, identify a unique molecular mechanism. A noisy linearly stable oscillator, a deterministic limit cycle with noise, or suitably structured inherited environmental variations can produce superficially similar spectra.

Similarly, a finite record of independent or correlated noise can show spurious bumps in its estimated power spectrum. Robust claims require sufficiently long records, multiple lineages/replicates, uncertainty estimates, and attention to whether the data were sampled in generation count or in physical time.

> **Zero and nonzero frequencies.** The zero-frequency component controls long-time accumulation of generation-time noise. Finite-frequency structure diagnoses dynamical modes. Both emerge from the same correlation function, but they answer different questions.

---

<div id="problems-for-discussion-or-homework"></div>

## Supplementary Derivations / Further Reading {#problems-further-reading}

### 1. What does it take to remove size memory?

For $x_{n+1}=ax_n+\xi_n$, size memory is described by $\mathbb E[x_{n+k}\mid x_n]$, the stationary variance, and the autocorrelation function. The feedback coefficient $a$ determines how perturbations decay; $a\lt 0$ also makes their signs alternate.

**Result.** $\mathbb E[x_{n+k}\mid x_n]=a^kx_n$, $\operatorname{Var}(x)=\sigma_\xi^2/(1-a^2)$, and $\rho(k)=a^{|k|}$. Finite stationary variance requires $|a|\lt 1$ for independent noise of finite nonzero variance. Negative $a$ leads to alternating-sign correlations.

### 2. Correlation without changing a histogram

The two-state Markov model has the same marginal distribution $\tfrac12\delta_{\tau_s}+\tfrac12\delta_{\tau_\ell}$ for every $0\lt p\lt 1$, but the correlations in (D.3) vary with $p$. Thus the growth correction in (4.4) also varies, and its approximation fails as $p\to1^-$.

**Result.** The transition matrix is symmetric, so its stationary probabilities are $1/2,1/2$. Correlations decay as $(2p-1)^k$. Therefore $D_\tau=\delta^2p/(1-p)$. As $p\to1^-$, the correlation time diverges and weak-noise/truncated-cumulant expansions cease to be uniform; use the exact matrix condition instead.

### 3. Can noisy lifetimes have zero long-run noise?

For a mixing stationary size process with finite $\operatorname{Var}(\ln v)$, equation (3.5) keeps $\operatorname{Var}(T_n)$ bounded as $n\to\infty$, without requiring each $\tau_n$ to be deterministic.

**Result.** $T_n-n\ln2/\lambda=(\ln v_n-\ln v_0)/\lambda$. Its variance tends to $2\operatorname{Var}(\ln v)/\lambda^2$ if endpoints decorrelate. Individual times are noisy, but their integrated correlations cancel the variance accumulation.

### 4. ACF versus PSD

The AR(1) spectrum can also be derived directly from its autocovariance function. Combined with (F.6), this explains the zero of the generation-time spectrum at $\omega=0$.

**Result.** Sum the geometric series $\sum_{k=-\infty}^{\infty}a^{|k|}e^{-i\omega k}=(1-a^2)/(1+a^2-2a\cos\omega)$; multiply by $\sigma_x^2=\sigma_\xi^2/(1-a^2)$. Differencing multiplies each Fourier component by $e^{i\omega}-1$, which vanishes at zero frequency.

### 5. Which part of a spectrum matters for growth? {#problem-5-which-part-of-a-spectrum-matters-for-growth}

Two stochastic lineage processes with the same $\bar\tau$ and $\operatorname{Var}(\tau)$ can have different growth rates because their correlation functions differ. A pronounced nonzero-frequency peak alone does not determine growth; in the quadratic weak-noise approximation, the correction depends on $S_\tau(0)$.

**Result.** Finite-frequency peaks and zero-frequency noise answer different questions. At quadratic order in a valid weak-noise cumulant expansion, the growth correction is proportional to $S_\tau(0)$, not to the height of a nonzero-frequency peak. Beyond this approximation, higher cumulant rates can matter.

### Further reading

For the experimental motivation, start with Taheri-Araghi et al. (2015) and Hashimoto et al. (2016). For the correlated-history formulation, see [Pigolotti (2021)](https://doi.org/10.1103/PhysRevE.103.L060402), read alongside Appendix C's distinction between expected and typical growth. [Lin and Amir (2020)](https://doi.org/10.1103/PhysRevE.101.012401) develops the effect of single-cell growth-rate variability on population growth. The mechanistic references below provide background on why phenomenological size-control slopes do not uniquely identify molecular regulation.

### References

- Taheri-Araghi, S., et al. (2015). *Cell-size control and homeostasis in bacteria*. **Current Biology, 25**(3), 385–391. [doi:10.1016/j.cub.2014.12.009](https://doi.org/10.1016/j.cub.2014.12.009). Original Figures 1–3 reproduced as Figs 2, 4, and 13.
- Ho et al., *Annual Review of Biophysics* (2018): cell-size regulation, lineage correlations, and fluctuations.
- Godin et al., *Nature Methods* (2010): single-cell growth measurements.
- Hashimoto, M., et al. (2016). *Noise-driven growth rate gain in clonal cellular populations*. **PNAS, 113**(12), 3251–3256. [doi:10.1073/pnas.1519412113](https://doi.org/10.1073/pnas.1519412113). Original Figure 3 reproduced as Fig 7.
- Wallden et al., *Cell* (2016); Ho and Amir, *Frontiers in Microbiology* (2015); Micali et al., *Cell Reports* (2018): replication-linked models and the distinction between phenomenological size control and molecular implementation.
- [Ho laboratory materials](https://ho.lab.westlake.edu.cn/teaching_en.html)
- Pigolotti, S. (2021). *Generalized Euler–Lotka equation for correlated cell divisions*. **Physical Review E, 103**, L060402. [doi:10.1103/PhysRevE.103.L060402](https://doi.org/10.1103/PhysRevE.103.L060402).
- Lin, J., and Amir, A. (2020). *From single-cell variability to population growth*. **Physical Review E, 101**, 012401. [doi:10.1103/PhysRevE.101.012401](https://doi.org/10.1103/PhysRevE.101.012401).
