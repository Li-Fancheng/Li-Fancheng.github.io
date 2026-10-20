"""Original analytic figures for the hierarchical-data review.

The curves illustrate the stated analytic/coverage/annealed models; they are not
measurements from trained networks or reproductions of published figure data.
"""

from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq
from scipy.special import gamma

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelcolor": "#24364b", "text.color": "#24364b",
    "axes.edgecolor": "#a5afba", "xtick.color": "#536475",
    "ytick.color": "#536475", "svg.fonttype": "none",
    "figure.facecolor": "white", "axes.facecolor": "white",
})
COLORS = ["#2563a6", "#c16436", "#268575", "#8a5fa7"]
policy_path = ROOT.parents[1] / "data" / "hierarchical_data_figures.json"
policy = json.loads(policy_path.read_text()) if policy_path.exists() else {}
checks = {}


def save(fig, name, data):
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=180, bbox_inches="tight")
    from PIL import Image
    width, height = Image.open(OUT / f"{name}.png").size
    policy[f"figures/{name}"] = {
        "width": 1000, "kind": "original-analytic",
        "intrinsic_width": width, "intrinsic_height": height,
    }
    (OUT / f"{name}.json").write_text(json.dumps(data, indent=2) + "\n")
    plt.close(fig)


# A generated-tree schematic and exact combinatorial counts.
fig, axs = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)
ax = axs[0]
for depth in range(4):
    xs = (np.arange(2**depth) + .5) / 2**depth
    if depth:
        parents = (np.arange(2**(depth-1)) + .5) / 2**(depth-1)
        for j, x in enumerate(xs):
            ax.plot([parents[j//2], x], [4-depth, 3-depth], color="#b6c4cf", lw=1.4, zorder=1)
    ax.scatter(xs, np.full(len(xs), 3-depth), s=330 if depth < 3 else 210,
               color=COLORS[depth], edgecolor="white", zorder=3)
    ax.text(-.03, 3-depth, f"level {3-depth}", va="center", ha="right", fontsize=10)
    for j, x in enumerate(xs):
        label = "$Y$" if depth == 0 else (f"$h_{j+1}$" if depth < 3 else f"$x_{j+1}$")
        ax.text(x, 3-depth, label, ha="center", va="center", color="white", fontsize=9)
ax.text(.5, -.55, r"one local rule: $a\to(b,c)$ or $a\to(d,e)$", ha="center", fontsize=10)
ax.set(xlim=(-.12, 1.05), ylim=(-.85, 3.5), title="A   Reusable local rules")
ax.axis("off")
ax = axs[1]
levels = np.arange(1, 11)
d = 2.**levels
log_samples = np.log10(16) + (d-1)*np.log10(4)
log_dictionary = np.log10(levels*16*4*2)
ax.plot(d, log_samples, "o-", color=COLORS[0], label="complete samples")
ax.plot(d, log_dictionary, "s-", color=COLORS[1], label="rule entries")
ax.set(xscale="log", xlabel=r"input length $d=2^L$", ylabel=r"$\log_{10}$ count",
       title="B   Combinations vs. dictionary size")
ax.legend(frameon=False, loc="upper left")
ax.grid(alpha=.18)
save(fig, "grammar", {"s": 2, "m": 4, "v": 16, "n_c": 16,
    "L": levels.tolist(), "d": d.tolist(), "log10_samples": log_samples.tolist(),
    "log10_rule_entries": log_dictionary.tolist()})


# Published large-v correlation formula, evaluated analytically.
v, m, s, L = 32, 4, 2, 9
rho = m/v**(s-1)
distance = np.arange(1, s**L)
ell = np.floor(np.log(distance)/np.log(s)).astype(int)+1
corr = np.sqrt((1-rho)/(v**3*m**(2.*ell-1)))
fig, axs = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)
ax = axs[0]
ax.step(distance, corr, where="post", color=COLORS[0], lw=2, label="correlation plateaus")
for P, color in zip([1e4, 1e7, 1e10], COLORS[1:]):
    ax.axhline(1/(v*np.sqrt(P)), color=color, ls="--", label=fr"noise: $P=10^{{{int(np.log10(P))}}}$")
ax.set(xscale="log", yscale="log", xlabel=r"token distance $t$", ylabel=r"$\widetilde C(t)$", title="A   Signal and sampling noise")
ax.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=1, fontsize=9)
ax.grid(alpha=.18, which="both")
P_values = np.logspace(2, 14, 500)
thresholds = v*m**(2.*np.arange(1,L+1)-1)/(1-rho)
resolved = (P_values[:,None] >= thresholds).sum(axis=1)
context = s**resolved-1
ax = axs[1]
ax.step(P_values, np.maximum(context, .5), where="post", color=COLORS[0], lw=2)
envelope = s**((np.log(P_values*(1-rho)/v)/np.log(m)+1)/2)-1
ax.plot(P_values, np.clip(envelope, .5, s**L-1), "--", color=COLORS[1], label="continuous envelope")
ax.set(xscale="log", yscale="log", xlabel=r"independent sequences $P$",
       ylabel="resolvable context length", title="B   Effective context at finite depth")
ax.set_yticks([1, 3, 15, 63, 255, 511], ["1", "3", "15", "63", "255", "511"])
ax.legend(frameon=False, fontsize=9)
ax.grid(alpha=.18, which="both")
checks["correlation_threshold_max_relative_error"] = float(np.max(np.abs(
    np.sqrt((1-rho)/(v**3*m**(2.*np.arange(1,L+1)-1))) * v*np.sqrt(thresholds)-1)))
save(fig, "correlations", {"v": v, "m": m, "s": s, "L": L,
    "thresholds": thresholds.tolist(), "distance": distance.tolist(),
    "correlation": corr.tolist(), "P": P_values.tolist(), "effective_context": context.tolist()})


# Exact finite-dictionary missing mass, with asymptotic power laws.
P_values = np.logspace(0, 10, 220)
rank = np.arange(1, 200001, dtype=float)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)
rare_data = {"dictionary_size": len(rank), "P": P_values.tolist(), "curves": {}}
for q, color in zip([1.3, 1.5, 2.0], COLORS):
    freq = rank**(-q)
    freq /= freq.sum()
    missing = np.array([np.sum(freq*np.exp(P*np.log1p(-freq))) for P in P_values])
    exponent = (q-1)/q
    c = freq[0]
    asymptotic = c**(1/q)/q*gamma(1-1/q)*P_values**(-exponent)
    axs[0].loglog(P_values, missing, color=color, label=fr"$q={q}$, $\alpha={exponent:.2f}$")
    axs[0].loglog(P_values, asymptotic, color=color, ls="--", lw=1, alpha=.65)
    rare_data["curves"][str(q)] = {"missing_mass": missing.tolist(), "alpha": exponent}
    assert np.all(np.diff(missing) <= 1e-14)
    if q == 1.5:
        for B, col in zip([1, 8, 64], COLORS):
            error = -np.expm1(B*np.log1p(-missing))
            axs[1].loglog(P_values, error, color=col, label=fr"$B={B}$")
axs[0].set(xlabel=r"samples $P$", ylabel="unseen probability mass", title="A   Rare-rule coverage")
axs[1].set(xlabel=r"samples $P$ (coverage model, $A=1$)", ylabel="at least one uncovered component",
           title=r"B   Accumulated error opportunities ($q=1.5$)")
for ax in axs:
    ax.set_ylim(1e-5, 1.3)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=.18, which="both")
save(fig, "rare_rules", rare_data)


# Finite-alphabet annealed ratio of expected rule weights.
v, m, s = 32, 8, 2
a = (m-1)/(v**s-1)
b = (m*v-1)/(v**s-1)
def F(p):
    return (p**s+a*(1-p**s))/(p**s+b*(1-p**s))
def derivative(p):
    return s*p**(s-1)*(b-a)/(b+(1-b)*p**s)**2
pstar = brentq(lambda p:F(p)-p, 1/v+1e-6, 1-1e-6)
lam = derivative(pstar)
eps_star = (1-pstar)/(1-1/v)
pgrid = np.linspace(1/v, 1, 1000)
epsilon = np.linspace(0, 1, 1201)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)
axs[0].plot(pgrid, F(pgrid), color=COLORS[0], label=r"$F(p)$")
axs[0].plot(pgrid, pgrid, "--", color="#8896a4", label=r"$p$")
axs[0].scatter([1/v,pstar,1], [1/v,pstar,1], color=[COLORS[2],COLORS[1],COLORS[2]], zorder=5)
axs[0].set(xlabel="correct-symbol message mass p", ylabel="updated mass",
           title="A   Stable branches and separator")
axs[0].legend(frameon=False)
depths=[1,3,8,24]
denoising_data={"v":v,"m":m,"s":s,"p_star":pstar,"epsilon_star":eps_star,
                "lambda":lam,"nu":np.log(s)/np.log(lam), "epsilon":epsilon.tolist(),"curves":{}}
for depth, color in zip(depths, COLORS):
    p = 1-epsilon+epsilon/v
    for _ in range(depth): p=F(p)
    axs[1].plot(epsilon,p,color=color,label=fr"$L={depth}$")
    denoising_data["curves"][str(depth)] = p.tolist()
    assert np.all(np.diff(p) <= 1e-12)
axs[1].axvline(eps_star,color="#8896a4",ls="--",lw=1)
axs[1].axhline(1/v,color="#8896a4",ls=":",lw=1)
axs[1].set(xlabel=r"soft-evidence noise $\epsilon$", ylabel="root correct-symbol mass",
           title="B   Repeated coarse-graining sharpens the crossover")
axs[1].legend(frameon=False)
for ax in axs: ax.grid(alpha=.18)
checks.update({
    "F_uniform_fixed_point_error": abs(F(1/v)-1/v),
    "F_correct_fixed_point_error": abs(F(1)-1),
    "F_separator_error": abs(F(pstar)-pstar),
    "F_derivative_at_correct": derivative(1),
    "F_derivative_formula_error": abs(derivative(1)-s*m*(v-1)/(v**s-1)),
    "unstable_derivative": lam,
    "mean_field_nu": np.log(s)/np.log(lam),
})
save(fig,"denoising",denoising_data)


# Independent check: upward/downward BP against exhaustive enumeration on a
# depth-two binary grammar (24 valid labelled derivations).
rng = np.random.default_rng(20261020)
vv, mm = 3, 2
rules = [rng.choice(vv*vv, size=vv*mm, replace=False).reshape(vv,mm) for _ in range(2)]
rules = [np.stack([r//vv,r%vv],axis=-1) for r in rules]
evidence = rng.uniform(.1,1,(4,vv))
low = np.array([np.prod(evidence[2*j+np.arange(2),rules[0][a,k]])
                for j in range(2) for a in range(vv) for k in range(mm)]).reshape(2,vv,mm).mean(-1)
up_root = np.array([np.mean([low[0,b]*low[1,c] for b,c in rules[1][a]]) for a in range(vv)])
bp_root = up_root/up_root.sum()
enumerated_root = np.zeros(vv)
enumerated_leaves = np.zeros((4,vv))
for root in range(vv):
    for b,c in rules[1][root]:
        for left in rules[0][b]:
            for right in rules[0][c]:
                leaves=np.r_[left,right]
                weight=np.prod(evidence[np.arange(4),leaves])/(vv*mm**3)
                enumerated_root[root]+=weight
                enumerated_leaves[np.arange(4),leaves]+=weight
z=enumerated_root.sum()
enumerated_root/=z
enumerated_leaves/=z
down_mid=np.zeros((2,vv))
for root in range(vv):
    for b,c in rules[1][root]:
        down_mid[0,b]+=low[1,c]/(vv*mm)
        down_mid[1,c]+=low[0,b]/(vv*mm)
down_leaf=np.zeros((4,vv))
for j in range(2):
    for a in range(vv):
        for b,c in rules[0][a]:
            down_leaf[2*j,b]+=down_mid[j,a]*evidence[2*j+1,c]/mm
            down_leaf[2*j+1,c]+=down_mid[j,a]*evidence[2*j,b]/mm
bp_leaf=down_leaf*evidence
bp_leaf/=bp_leaf.sum(axis=1,keepdims=True)
checks["BP_root_vs_enumeration_max_error"]=float(np.max(np.abs(bp_root-enumerated_root)))
checks["BP_leaf_vs_enumeration_max_error"]=float(np.max(np.abs(bp_leaf-enumerated_leaves)))
assert checks["BP_root_vs_enumeration_max_error"] < 1e-12
assert checks["BP_leaf_vs_enumeration_max_error"] < 1e-12
assert checks["F_uniform_fixed_point_error"] < 1e-12
assert checks["correlation_threshold_max_relative_error"] < 1e-12
assert derivative(1/v)<1 and derivative(1)<1 and lam>1
(ROOT / "calculation_checks.json").write_text(json.dumps(checks,indent=2)+"\n")
(ROOT.parents[1] / "data" / "hierarchical_data_figures.json").write_text(json.dumps(policy,indent=2)+"\n")
print(json.dumps(checks,indent=2))
