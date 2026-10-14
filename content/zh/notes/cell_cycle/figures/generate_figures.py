"""Reproduce the note's model illustrations: conda run -n env python figures/generate_figures.py."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
from scipy.optimize import brentq

OUT = Path(__file__).resolve().parent
BLUE, ORANGE, GREEN, PURPLE = "#2767a0", "#d47725", "#238579", "#8057a5"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.titlesize": 13, "axes.labelsize": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "legend.frameon": False, "lines.linewidth": 2.4,
})


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=190, bbox_inches="tight")
    plt.close(fig)


def finish(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=.16)


# 1. Exact log-timer variance versus an illustrative restoring log-size map.
fig, axes = plt.subplots(1, 2, figsize=(11, 4.1), layout="constrained")
g = np.arange(61)
rng = np.random.default_rng(14)
for j in range(7):
    eps = rng.normal(0, .10, 60)
    timer = np.r_[0, np.cumsum(eps)]
    controlled = np.zeros(61)
    for n in range(60):
        controlled[n+1] = .5 * controlled[n] + eps[n]
    axes[0].plot(g, timer, color=ORANGE, alpha=.35, lw=1.2)
    axes[0].plot(g, controlled, color=BLUE, alpha=.50, lw=1.2)
axes[0].plot([], [], color=ORANGE, label="Independent timer")
axes[0].plot([], [], color=BLUE, label="Restoring log-size map")
axes[0].axhline(0, color="#777777", lw=.7)
axes[0].set_title("A  |  Illustrative lineage trajectories")
axes[0].legend(loc="lower left", fontsize=9)
finish(axes[0], "Generation n", r"Log-size deviation $y_n-y_*$")
axes[1].plot(g, g, color=ORANGE, label=r"Timer: $n$")
axes[1].plot(g, (1-.5**(2*g))/(1-.5**2), color=BLUE,
             label=r"Restoring map: $(1-a^{2n})/(1-a^2)$")
axes[1].set_title("B  |  Variance from a fixed initial size")
axes[1].legend(fontsize=9)
finish(axes[1], "Generation n", "Variance / innovation variance")
save(fig, "01_timer_and_homeostasis")


# 2. Exact two-state Perron-eigenvalue condition, using the note's parameters.
def growth(p):
    def residual(rate):
        A, B = np.exp(-rate*np.array([2., 5.]))
        return p*(A+B) + np.sqrt(p*p*(A-B)**2 + 4*(1-p)**2*A*B) - 1
    return brentq(residual, 0, 1)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.1), layout="constrained")
k = np.arange(9)
for p, color, label in [( .25, BLUE, "Alternating tendency"),
                         (.5, GREEN, "Independent"), (.75, ORANGE, "Persistent")]:
    axes[0].plot(k, (2*p-1)**k, "o-", color=color, label=f"p = {p:.2f}: {label}")
axes[0].axhline(0, color="#777777", lw=.8)
axes[0].set_title("A  |  Same lifetimes, different memory")
axes[0].legend(fontsize=9)
finish(axes[0], "Lag k (generations)", r"Generation-time correlation $\rho_\tau(k)$")
ps = np.linspace(0, .995, 200)
axes[1].plot(ps, [growth(p) for p in ps], color=PURPLE, label="Exact growth rate")
for p in [0, .25, .5, .75]:
    axes[1].plot(p, growth(p), "o", color=PURPLE)
axes[1].axhline(np.log(2)/3.5, color="#777777", ls="--", lw=1.4,
                label=r"Mean-time baseline $\ln 2/\mu$")
axes[1].set_title(r"B  |  $\tau_s=2$ h, $\tau_\ell=5$ h")
axes[1].legend(fontsize=9)
finish(axes[1], "State-retention probability p", r"Population growth $\Lambda$ (h$^{-1}$)")
save(fig, "02_correlations_and_growth")


# 3. A geometric diagram of the noise-free adder, with lengths encoding size.
fig, ax = plt.subplots(figsize=(10.5, 4.4), layout="constrained")
ax.set_xlim(0, 10.5)
ax.set_ylim(-.1, 3.5)
ax.axis("off")
ax.set_title("The adder: equal added size brings daughters closer to the mean", pad=12)

def cell(x, y, size, color, label):
    ax.add_patch(FancyBboxPatch((x, y), size*.8, .38,
                 boxstyle="round,pad=0,rounding_size=0.19", ec=color, fc=color, alpha=.8))
    ax.text(x+size*.4, y-.13, label, ha="center", va="top", fontsize=10)

for y, v, tag in [(2.5, .7, "Small birth"), (1.05, 1.3, "Large birth")]:
    ax.text(.0, y+.6, tag, color=BLUE, weight="bold")
    cell(.2, y, v, BLUE, f"Birth: {v:.2f} Δ")
    ax.annotate("", (3.2, y+.19), (1.5, y+.19), arrowprops=dict(arrowstyle="->", color="#555555"))
    ax.text(2.35, y+.42, "+ Δ", ha="center", color=ORANGE, weight="bold")
    cell(3.3, y, v+1, ORANGE, f"Division: {v+1:.2f} Δ")
    ax.annotate("", (7.1, y+.19), (5.5, y+.19), arrowprops=dict(arrowstyle="->", color="#555555"))
    ax.text(6.3, y+.42, "Split equally", ha="center", fontsize=10)
    daughter = (v+1)/2
    cell(7.3, y+.25, daughter, GREEN, "")
    cell(7.3, y-.25, daughter, GREEN, f"Each: {daughter:.2f} Δ")
ax.text(5.3, -.05, r"Mean birth size = Δ     •     Deviation is halved: $x_{n+1}=x_n/2$",
        ha="center", fontsize=11)
save(fig, "03_adder_size_control")


# 4. Linearized adder time correlations and exact AR(1)-difference variance.
fig, axes = plt.subplots(1, 2, figsize=(11, 4.1), layout="constrained")
k = np.arange(9)
rho_time = np.r_[1., -2.**(-(k[1:]+1))]
axes[0].plot(k, .5**k, "o-", color=BLUE, label="Birth size (exact adder map)")
axes[0].plot(k, rho_time, "s-", color=ORANGE, label="Generation time (linearized)")
axes[0].axhline(0, color="#777777", lw=.8)
axes[0].set_title("A  |  Size memory creates time anticorrelation")
axes[0].legend(fontsize=9)
finish(axes[0], "Lag k (generations)", "Normalized autocorrelation")
n = np.arange(1, 41)
axes[1].plot(n, n, color=ORANGE, label="Independent times")
axes[1].plot(n, 2*(1-.5**n), color=BLUE, label="Linearized adder times")
axes[1].set_title("B  |  Matched one-generation time variance")
axes[1].legend(fontsize=9)
finish(axes[1], "Number of generations n", r"Accumulated variance $\mathrm{Var}(T_n)/C_\tau(0)$")
save(fig, "04_noise_cancellation")


# 5. Analytic spectra and ACFs; each PSD is scaled to unit process variance.
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), layout="constrained")
w = np.linspace(0, np.pi, 2001)
a = .5
sx = (1-a*a)/(1+a*a-2*a*np.cos(w))
st = 4*np.sin(w/2)**2*sx/(2*(1-a))
axes[0].plot(w/np.pi, sx, color=BLUE, label="Birth size")
axes[0].plot(w/np.pi, st, color=ORANGE, label="Generation time (linearized)")
axes[0].set_title("A  |  Feedback suppresses zero frequency")
axes[0].legend(fontsize=9)
finish(axes[0], r"Angular frequency $\omega/\pi$", "PSD / process variance")
r, w0 = .9, np.pi/3
ar2 = 1/abs(1-2*r*np.cos(w0)*np.exp(-1j*w)+r*r*np.exp(-2j*w))**2
ar2 /= np.trapezoid(ar2, w)/np.pi
axes[1].plot(w/np.pi, sx, color=BLUE, label="AR(1): relaxation")
axes[1].plot(w/np.pi, ar2, color=PURPLE, label="AR(2): damped oscillation")
axes[1].axvline(w0/np.pi, color=PURPLE, ls=":", lw=1, label=r"Mode angle $\omega_0/\pi$")
axes[1].set_title("B  |  Delayed feedback creates a peak")
axes[1].legend(fontsize=9)
finish(axes[1], r"Angular frequency $\omega/\pi$", "PSD / process variance")
k = np.arange(21)
b = (1-r*r)/(1+r*r)/np.tan(w0)
rho2 = r**k*(np.cos(k*w0)+b*np.sin(k*w0))
axes[2].plot(k, a**k, "o-", color=BLUE, ms=4, label="AR(1)")
axes[2].plot(k, rho2, "o-", color=PURPLE, ms=4, label="AR(2)")
axes[2].axhline(0, color="#777777", lw=.8)
axes[2].set_title("C  |  The corresponding memories")
axes[2].legend(fontsize=9)
finish(axes[2], "Lag k (generations)", "Normalized autocorrelation")
save(fig, "05_spectra_and_oscillations")

print("Generated five model figures in", OUT)
