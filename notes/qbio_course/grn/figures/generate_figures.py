"""Reproduce the GRN note's nine figures (model calculations, not data).

Requires Python 3.10+, NumPy, SciPy and Matplotlib. Run from any directory:
    python notes/qbio_course/grn/figures/generate_figures.py

SVG is used in the note; PNG is exported for inspection and reuse. Random
trajectories use fixed seeds. All rates and times are dimensionless.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import roots_jacobi
from scipy.stats import poisson


OUT = Path(__file__).resolve().parent
BLUE, ORANGE, GREEN, PURPLE = "#2767a0", "#d47725", "#238579", "#8057a5"
INK, MUTED = "#243747", "#637181"
plt.rcParams.update({
    "font.family": "STIXGeneral", "mathtext.fontset": "stix",
    "font.size": 12, "axes.titlesize": 14, "axes.titleweight": "bold",
    "axes.labelsize": 12, "legend.fontsize": 10.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "text.color": INK, "axes.labelcolor": INK,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "legend.frameon": False, "lines.linewidth": 2.3,
    "svg.fonttype": "path", "svg.hashsalt": "grn-theory",
})


def save(fig, name):
    svg = OUT / f"{name}.svg"
    fig.savefig(svg, bbox_inches="tight", metadata={"Date": None})
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines())+"\n")
    fig.savefig(OUT / f"{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def title(ax, text):
    ax.set_title(text, loc="left", pad=14)


def arrow(ax, start, end, color=INK, inhibit=False, **kwargs):
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle="-[,widthB=0.45,lengthB=0" if inhibit else "-|>",
        mutation_scale=12 if inhibit else 14, color=color, lw=1.8, **kwargs))


def node(ax, xy, label, color=BLUE, radius=.22):
    ax.add_patch(Circle(xy, radius, facecolor="white", edgecolor=color, lw=2))
    ax.text(*xy, label, ha="center", va="center", color=color, fontsize=14)


def fixed(ax, x, y, stable=True, color=INK):
    ax.scatter(x, y, s=64, zorder=7, edgecolors=color, linewidths=1.8,
               facecolors=color if stable else "white", clip_on=False)


def regulation():
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), layout="constrained")
    ax = axes[0]
    title(ax, "A  |  Feedback motifs")
    ax.set(xlim=(-.6, 3.8), ylim=(-.6, 3.4))
    ax.axis("off")
    node(ax, (.2, 2.5), "$X$")
    arrow(ax, (.32, 2.69), (.05, 2.71), BLUE, connectionstyle="arc3,rad=2.5")
    ax.text(1.25, 2.5, "Self-activation (+)", va="center")
    node(ax, (.0, 1.05), "$X$")
    node(ax, (1.35, 1.05), "$Y$", ORANGE)
    arrow(ax, (.2, 1.18), (1.12, 1.18), BLUE, inhibit=True)
    arrow(ax, (1.12, .92), (.2, .92), ORANGE, inhibit=True)
    ax.text(2.0, 1.05, "Toggle\n(−) × (−) = (+)", va="center", fontsize=11)
    pts = [(0, -.25), (1.35, -.25), (.68, .35)]
    for p, lab, c in zip(pts, ["$X_1$", "$X_2$", "$X_3$"], [BLUE, ORANGE, GREEN]):
        node(ax, p, lab, c, radius=.19)
    for i in range(3):
        p, q = np.array(pts[i]), np.array(pts[(i+1) % 3])
        d = (q-p) / np.linalg.norm(q-p)
        arrow(ax, p+.25*d, q-.25*d, MUTED, inhibit=True)
    ax.text(2.0, -.02, "Repression ring\n(−)³ = (−)", va="center", fontsize=11)
    x = np.linspace(0, 3, 500)
    for ax, positive, heading in zip(axes[1:], [True, False],
                                     ["B  |  Activation", "C  |  Repression"]):
        title(ax, heading)
        for n, c in [(1, BLUE), (2, GREEN), (4, ORANGE)]:
            h = x**n/(1+x**n)
            ax.plot(x, h if positive else 1-h, color=c, label=f"$n={n}$")
        ax.axvline(1, ls=":", lw=1.2, color=MUTED)
        ax.axhline(.5, ls=":", lw=1.2, color=MUTED)
        ax.set(xlim=(0, 3), ylim=(0, 1.05), xlabel="$x/K$",
               ylabel=r"$(r-r_{\min})/(r_{\max}-r_{\min})$")
        ax.legend(loc="upper left" if positive else "upper right")
        ax.grid(alpha=.12)
    save(fig, "01_regulation_and_hill")


def activation_parameters():
    basal, alpha = .02, 2.5
    roots = np.sort(np.roots([1, -(basal+alpha), 1, -basal])).real
    return basal, alpha, roots


def self_activation():
    basal, alpha, (low, threshold, high) = activation_parameters()
    drift = lambda x: basal+alpha*x*x/(1+x*x)-x
    fig, (ax, time) = plt.subplots(1, 2, figsize=(11.8, 4.4), layout="constrained")
    title(ax, "A  |  Production meets removal")
    x = np.linspace(0, 2.75, 600)
    ax.plot(x, basal+alpha*x*x/(1+x*x), color=BLUE, label=r"Production $F(x)$")
    ax.plot(x, x, color=MUTED, label="Removal $x$")
    for root, stable in [(low, True), (threshold, False), (high, True)]:
        fixed(ax, root, root, stable)
    ax.set(xlim=(0, 2.75), ylim=(-.43, 2.8), xlabel="Expression $x$", ylabel="Rate")
    ax.legend(loc="upper left")
    ax.text(.54, .25, r"Threshold $x_U$", fontsize=11)
    ax.text(2.7, .14, r"$b=0.02,\ \alpha=2.5,\ n=2$", ha="right", fontsize=11)
    ax.plot([0, 2.65], [-.25, -.25], color=MUTED, lw=.7)
    for l, r, right in [(low, threshold, False), (threshold, high, True), (high, 2.65, False)]:
        a, b = l+.05, r-.05
        arrow(ax, (a if right else b, -.25), (b if right else a, -.25))
    title(time, "B  |  Same parameters, different histories")
    ts = np.linspace(0, 16, 400)
    for x0 in [.08, .3, .44, .49, .9, 2.6]:
        sol = solve_ivp(lambda t, v: drift(v), (0, 16), [x0],
                        t_eval=ts, rtol=1e-8, atol=1e-10)
        time.plot(ts, sol.y[0], color=BLUE if x0 < threshold else ORANGE, alpha=.85)
    time.axhline(threshold, color=MUTED, ls="--", lw=1.3, label="Unstable threshold")
    time.text(15.5, 2.15, "High attractor", ha="right", color=ORANGE)
    time.text(15.5, .13, "Low attractor", ha="right", color=BLUE)
    time.set(xlim=(0, 16), ylim=(-.08, 2.75), xlabel="Time $t$", ylabel="$x(t)$")
    time.legend(loc="upper right")
    time.grid(alpha=.12)
    save(fig, "02_self_activation")


def hysteresis():
    fig, (ax, phase) = plt.subplots(1, 2, figsize=(11.8, 4.6), layout="constrained")
    title(ax, "A  |  Same self-activation model: hysteresis")
    basal = .02
    zc = 1/np.sqrt(3)
    fold_b = lambda x: x*(1-x*x)/2
    fold_a = lambda x: (1+x*x)**2/(2*x)
    xup = brentq(lambda x: fold_b(x)-basal, 1e-10, zc)
    xdown = brentq(lambda x: fold_b(x)-basal, zc, 1)
    aup, adown = fold_a(xup), fold_a(xdown)
    alpha_of = lambda x: (x-basal)*(1+x*x)/(x*x)
    ax.axvspan(adown, aup, color=GREEN, alpha=.08)
    for l, r, stable, label in [(basal+1e-8, xup, True, "Stable"),
                                (xup, xdown, False, "Unstable"),
                                (xdown, 15, True, None)]:
        xs = np.geomspace(l, r, 700)
        ax.plot(alpha_of(xs), xs, color=INK, ls="-" if stable else "--", label=label)
    high_at_up = np.sort(np.roots([1, -(basal+aup), 1, -basal])).real[-1]
    low_at_down = np.sort(np.roots([1, -(basal+adown), 1, -basal])).real[0]
    for alpha, start, end, color in [(aup, xup, high_at_up, ORANGE),
                                     (adown, xdown, low_at_down, BLUE)]:
        ax.axvline(alpha, color=MUTED, ls=":", lw=1)
        arrow(ax, (alpha, start), (alpha, end), color)
    for xs, color in [(np.array([.024, .029]), ORANGE), (np.array([8., 5.]), BLUE)]:
        arrow(ax, (alpha_of(xs[0]), xs[0]), (alpha_of(xs[1]), xs[1]), color)
    ax.text(7, .25, "History-dependent state", ha="center", fontsize=11)
    ax.set(xlim=(1.5, 14), ylim=(.017, 16), yscale="log",
           xlabel=r"Activation strength $\alpha$", ylabel=r"Expression $x^*$ (log scale)")
    ax.text(.70, .60, r"$b=0.02,\ n=2$", transform=ax.transAxes, ha="right", va="top")
    ax.legend(loc="upper left")
    title(phase, "B  |  Self-activation: parameter regions")
    bc = 1/(3*np.sqrt(3))
    bs = np.linspace(.0001, bc-1e-8, 700)
    xs_low = np.array([brentq(lambda x: fold_b(x)-b, 1e-10, zc) for b in bs])
    xs_high = np.array([brentq(lambda x: fold_b(x)-b, zc, 1) for b in bs])
    lower, upper = fold_a(xs_high), fold_a(xs_low)
    phase.fill_between(bs, lower, upper, color=GREEN, alpha=.17)
    phase.plot(bs, lower, color=GREEN, label="Saddle-node boundaries")
    phase.plot(bs, upper, color=GREEN)
    phase.plot(bc, fold_a(zc), "o", color=GREEN)
    phase.annotate("Cusp", (bc, fold_a(zc)), xytext=(.216, 1.1),
                   arrowprops={"arrowstyle": "->", "color": MUTED}, fontsize=11)
    phase.text(.035, 2.6, "Bistable", color=GREEN)
    phase.text(.17, 3.5, "One high state", ha="center")
    phase.text(.08, .8, "One low state", ha="center")
    phase.scatter(.02, 2.5, color=PURPLE, zorder=6)
    phase.annotate(r"$(0.02,\ 2.5)$", (.02, 2.5), xytext=(.03, 3.35),
                   arrowprops={"arrowstyle": "->", "color": PURPLE}, color=PURPLE)
    phase.set(xlim=(0, .26), ylim=(.5, 4.1), xlabel="Basal production $b$",
              ylabel=r"Activation strength $\alpha$")
    phase.legend(loc="upper right")
    save(fig, "03_hysteresis_and_phase_diagram")


def toggle():
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), layout="constrained")
    for ax, alpha, label in zip(axes[:2], [1.5, 3],
                                ["A  |  One attractor", "B  |  Two attractors"]):
        title(ax, label)
        end = alpha+.2
        vals = np.linspace(0, end, 38)
        x, y = np.meshgrid(vals, vals)
        if alpha > 2:
            ax.contourf(x, y, x-y, levels=[-end, 0, end],
                        colors=[ORANGE, BLUE], alpha=.07)
            ax.plot([0, end], [0, end], ls="--", color=PURPLE, lw=1.4,
                    label="Basin boundary")
        ax.streamplot(vals, vals, alpha/(1+y*y)-x, alpha/(1+x*x)-y,
                      color="#c5cdd3", density=.8, linewidth=.8, arrowsize=.8)
        xx = np.linspace(0, end, 500)
        ax.plot(alpha/(1+xx*xx), xx, color=BLUE, label=r"$\dot x=0$")
        ax.plot(xx, alpha/(1+xx*xx), color=ORANGE, label=r"$\dot y=0$")
        s = brentq(lambda s: s*(1+s*s)-alpha, 0, alpha)
        fixed(ax, s, s, alpha < 2)
        if alpha > 2:
            hi, lo = (alpha+np.sqrt(alpha**2-4))/2, (alpha-np.sqrt(alpha**2-4))/2
            fixed(ax, [hi, lo], [lo, hi])
        for initial, color in [([.2, end-.2], ORANGE), ([end-.2, .2], BLUE),
                               ([end-.4, end-.7], BLUE)]:
            sol = solve_ivp(lambda t, v: [alpha/(1+v[1]**2)-v[0],
                                          alpha/(1+v[0]**2)-v[1]],
                            (0, 25), initial, rtol=1e-8, atol=1e-10)
            ax.plot(*sol.y, color=color, lw=1.5, alpha=.8)
        ax.set(xlim=(0, end), ylim=(0, end), xlabel="$x$", ylabel="$y$")
        ax.set_aspect("equal")
        ax.text(.96, .95, rf"$\alpha={alpha:g},\ n=2$", transform=ax.transAxes,
                ha="right", va="top")
        ax.legend(loc="upper right", bbox_to_anchor=(1, -.17), ncol=2, fontsize=10)
    ax = axes[2]
    title(ax, "C  |  Symmetry breaking")
    a = np.linspace(2, 4, 400)
    ax.plot([.5, 2], [0, 0], color=INK, label="Stable")
    ax.plot([2, 4], [0, 0], color=INK, ls="--", label="Unstable")
    ax.plot(a, np.sqrt(a*a-4), color=BLUE)
    ax.plot(a, -np.sqrt(a*a-4), color=ORANGE)
    ax.axvline(2, ls=":", lw=1, color=MUTED)
    ax.text(2.1, .28, r"$\alpha_c=2$")
    ax.text(3.85, 2.3, "$X$ high", color=BLUE, ha="right")
    ax.text(3.85, -2.3, "$Y$ high", color=ORANGE, ha="right")
    ax.set(xlim=(.5, 4), ylim=(-3.8, 3.8), xlabel=r"$\alpha$", ylabel="$x^*-y^*$")
    ax.legend(loc="upper left")
    save(fig, "04_toggle_phase_portraits")


def oscillator():
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8), layout="constrained")
    solutions = {}
    for ax, alpha, label in zip(axes[0], [3, 6],
                                ["A  |  Damped response", "B  |  Sustained oscillation"]):
        title(ax, label)
        ts = np.linspace(0, 90, 4500)
        sol = solve_ivp(lambda t, v: alpha/(1+np.roll(v, 1)**3)-v,
                        (0, 90), [1.1, 1.2, .9], t_eval=ts, rtol=1e-9, atol=1e-11)
        solutions[alpha] = sol
        for j, c in enumerate([BLUE, ORANGE, GREEN]):
            ax.plot(ts, sol.y[j], color=c, lw=1.8, label=rf"$x_{j+1}$")
        ax.set(xlim=(0, 45), ylim=(.5, 3.1), xlabel="Time $t$", ylabel="Expression")
        ax.text(.98, .96, rf"$n=3,\ \alpha={alpha}$", transform=ax.transAxes,
                ha="right", va="top")
        ax.legend(loc="upper right", bbox_to_anchor=(1, .86), ncol=3)
        ax.grid(alpha=.12)
    ax = axes[1, 0]
    title(ax, "C  |  Trajectories projected onto $(x_1,x_2)$")
    for alpha, color in [(3, BLUE), (6, ORANGE)]:
        sol = solutions[alpha]
        mask = sol.t >= 60 if alpha == 6 else np.ones(len(sol.t), dtype=bool)
        ax.plot(sol.y[0, mask], sol.y[1, mask], color=color,
                label=rf"$\alpha={alpha}$", lw=2)
        s = brentq(lambda s: s*(1+s**3)-alpha, 0, alpha)
        fixed(ax, s, s, alpha == 3, color)
    ax.set(xlim=(.5, 3.1), ylim=(.5, 3.1), xlabel="$x_1$", ylabel="$x_2$")
    ax.set_aspect("equal")
    ax.legend(loc="upper left")
    ax.text(.97, .04, "3D dynamics; 2D projection", transform=ax.transAxes,
            ha="right", fontsize=10, color=MUTED)
    phase = axes[1, 1]
    title(phase, "D  |  Linear stability boundary")
    ns = np.linspace(2.001, 6, 900)
    ah = ns/(ns-2)*(2/(ns-2))**(1/ns)
    phase.axvspan(1, 2, color=BLUE, alpha=.07)
    phase.fill_between(ns, 0, np.minimum(ah, 12), color=BLUE, alpha=.07)
    phase.fill_between(ns, ah, 12, where=ah <= 12, color=ORANGE, alpha=.13)
    phase.plot(ns, ah, color=PURPLE, label=r"$\alpha_H(n)$")
    phase.scatter([3, 3], [3, 6], color=[BLUE, ORANGE], zorder=6)
    phase.text(3.13, 3, "A", va="center", color=BLUE)
    phase.text(3.13, 6, "B", va="center", color=ORANGE)
    phase.text(4.2, 8, "Unstable fixed point", ha="center", color=ORANGE)
    phase.text(2.9, .6, "Stable fixed point", ha="center", color=BLUE)
    phase.set(xlim=(1, 6), ylim=(0, 12), xlabel="Hill exponent $n$",
              ylabel=r"Repression strength $\alpha$")
    phase.legend(loc="upper right")
    save(fig, "05_repressilator")
    # A numerical example of attraction, not a general proof of Hopf criticality.
    assert np.ptp(solutions[3].y[:, -500:], axis=1).max() < 1e-3
    assert np.ptp(solutions[6].y[:, -1000:], axis=1).min() > .5


def telegraph_trace(total=45, seed=28):
    rng = np.random.default_rng(seed)
    t, state, count = 0., 0, 0
    times, states, counts = [t], [state], [count]
    while t < total:
        # Switching rate 0.2 in both directions; r=20, gamma=1.
        rates = np.array([.2, 20*state, count], dtype=float)
        t += rng.exponential(1/rates.sum())
        if t > total:
            break
        event = rng.choice(3, p=rates/rates.sum())
        if event == 0:
            state = 1-state
        else:
            count += 1 if event == 1 else -1
        times.append(t)
        states.append(state)
        counts.append(count)
    times.append(total)
    states.append(state)
    counts.append(count)
    return np.array(times), np.array(states), np.array(counts)


def expression_noise():
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), layout="constrained")
    ax = axes[0]
    title(ax, "A  |  Promoter switching leaves a trace")
    t, s, m = telegraph_trace()
    ax.fill_between(t, 0, 33, where=s == 1, step="post", color=GREEN, alpha=.1,
                    label="Promoter ON")
    ax.step(t, m, where="post", color=BLUE, lw=1.3, label="mRNA count")
    ax.set(xlim=(0, 45), ylim=(0, 33), xlabel="Time $t$", ylabel="$M(t)$")
    ax.legend(loc="upper left")
    ax = axes[1]
    title(ax, "B  |  Same mean, different noise")
    count = np.arange(0, 41)
    ax.plot(count, poisson.pmf(count, 10), color=INK, label="Constant rate: Poisson")
    for k, color, label in [(.2, ORANGE, "Slow switching"), (20, BLUE, "Fast switching")]:
        # Exact stationary telegraph distribution: a beta mixture of Poissons.
        # Q~Beta(k_on/gamma,k_off/gamma), M|Q~Poisson((r/gamma)Q).
        nodes, weights = roots_jacobi(160, k-1, k-1)
        weights /= weights.sum()
        pmf = poisson.pmf(count[:, None], 20*(nodes[None, :]+1)/2) @ weights
        ax.plot(count, pmf, color=color, label=label)
        assert abs(np.dot(count, pmf)-10) < 1e-3
    ax.set(xlim=(0, 34), ylim=(0, .3), xlabel="mRNA number $M$", ylabel="Probability")
    ax.legend(loc="upper right", fontsize=9.5)
    ax = axes[2]
    title(ax, "C  |  Slow switching amplifies variance")
    kappa = np.logspace(-2, 3, 500)
    ax.semilogx(kappa, 1+10/(1+kappa), color=PURPLE)
    ax.axhline(1, color=INK, ls="--", lw=1.3, label="Poisson baseline")
    ax.scatter([.4, 40], [1+10/1.4, 1+10/41], color=[ORANGE, BLUE], zorder=6)
    ax.text(.97, .94, r"$\pi=1/2,\ r/\gamma=20$", transform=ax.transAxes,
            ha="right", va="top")
    ax.set(xlim=(.01, 1000), ylim=(0, 12), xlabel=r"$\kappa/\gamma$",
           ylabel=r"Fano factor $\mathrm{Var}(M)/\langle M\rangle$")
    ax.legend(loc="center right")
    ax.grid(alpha=.12)
    save(fig, "06_expression_noise")


def landscape():
    basal, alpha, (low, threshold, high) = activation_parameters()
    drift = lambda x: basal+alpha*x*x/(1+x*x)-x
    slope = lambda x: 2*alpha*x/(1+x*x)**2-1
    potential = lambda x: x*x/2-basal*x-alpha*(x-np.arctan(x))
    barrier = potential(threshold)-potential(high)
    prefactor = 2*np.pi/np.sqrt(slope(threshold)*-slope(high))
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), layout="constrained")
    ax = axes[0]
    title(ax, "A  |  Same model, two escape barriers")
    x = np.linspace(0, 3, 700)
    ax.plot(x, potential(x), color=INK, label="$U(x)$")
    for root, stable in [(low, True), (threshold, False), (high, True)]:
        fixed(ax, root, potential(root), stable)
    ax.axhline(potential(threshold), color=MUTED, ls=":", lw=1)
    ax.annotate("", xy=(high, potential(threshold)), xytext=(high, potential(high)),
                arrowprops={"arrowstyle": "<->", "color": ORANGE, "lw": 1.7})
    ax.text(1.88, -.12, r"$\Delta U_H\simeq0.297$", color=ORANGE, fontsize=11, ha="right")
    ax.annotate(r"$\Delta U_L\simeq0.0247$", (threshold, potential(threshold)),
                xytext=(.7, .14), fontsize=10,
                arrowprops={"arrowstyle": "->", "color": MUTED})
    ax.set(xlim=(-.08, 3), ylim=(-.34, .3), xlabel="Expression $x$", ylabel="Potential $U$")
    ax.text(.96, .95, r"$b=0.02,\ \alpha=2.5$", transform=ax.transAxes, ha="right", va="top", fontsize=11)
    ax = axes[1]
    title(ax, "B  |  Switching with a reflecting boundary")
    dt, duration, diffusion = .005, 700, .1
    rng = np.random.default_rng(74)
    ts = np.arange(0, duration+dt/2, dt)
    xs = np.empty(len(ts))
    xs[0] = high
    noise = rng.normal(0, np.sqrt(2*diffusion*dt), len(ts)-1)
    for i in range(len(ts)-1):
        # Reflected Euler approximation on the nonnegative concentration axis.
        xs[i+1] = abs(xs[i]+drift(xs[i])*dt+noise[i])
    ax.plot(ts[::8], xs[::8], color=BLUE, lw=.8)
    ax.axhline(threshold, color=MUTED, ls="--", lw=1)
    ax.axhline(high, color=ORANGE, ls=":", lw=1)
    ax.text(.97, .96, r"$D=0.1,\ x\geq0$", transform=ax.transAxes, ha="right", va="top")
    ax.set(xlim=(0, duration), ylim=(0, 3.6), xlabel="Time $t$", ylabel="$x(t)$")
    ax = axes[2]
    title(ax, "C  |  High-state lifetime: weak-noise limit")
    inverse_d = np.linspace(40, 160, 250)
    tau = prefactor*np.exp(barrier*inverse_d)
    ax.semilogy(inverse_d, tau, color=PURPLE)
    ax.text(.06, .94, r"$\log\tau_K=\log C+\Delta U_H/D$",
            transform=ax.transAxes, va="top", fontsize=12)
    ax.set(xlim=(40, 160), xlabel="Inverse noise strength $1/D$", ylabel="Kramers escape time")
    ax.grid(alpha=.12, which="both")
    save(fig, "07_landscape_and_switching")


def boolean_network():
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), layout="constrained")
    ax = axes[0]
    title(ax, "A  |  The Derrida map")
    d = np.linspace(0, .55, 500)
    ax.plot(d, d, color=INK, ls="--", lw=1.3, label="$d_{t+1}=d_t$")
    for k, c in [(1, BLUE), (2, GREEN), (3, ORANGE)]:
        ax.plot(d, .5*(1-(1-d)**k), color=c, label=rf"$K={k}$")
    ax.text(.97, .04, "$p=1/2$", transform=ax.transAxes, ha="right")
    ax.set(xlim=(0, .55), ylim=(0, .55), xlabel="$d_t$", ylabel="$d_{t+1}$")
    ax.set_aspect("equal")
    ax.legend(loc="upper left")
    ax = axes[1]
    title(ax, "B  |  A small perturbation evolves")
    for k, c in [(1, BLUE), (2, GREEN), (3, ORANGE)]:
        ds = [.05]
        for _ in range(35):
            ds.append(.5*(1-(1-ds[-1])**k))
        ax.plot(np.arange(len(ds)), ds, color=c, label=rf"$K={k}$")
    ax.set(xlim=(0, 35), ylim=(0, .42), xlabel="Update step $t$", ylabel="Hamming distance $d_t$")
    ax.legend(loc="upper right")
    ax.grid(alpha=.12)
    ax = axes[2]
    title(ax, "C  |  Annealed phase diagram")
    p = np.linspace(.001, .999, 1000)
    kc = 1/(2*p*(1-p))
    ax.fill_between(p, 0, np.minimum(kc, 8), color=BLUE, alpha=.09)
    ax.fill_between(p, kc, 8, where=kc <= 8, color=ORANGE, alpha=.16)
    ax.plot(p, kc, color=PURPLE, label=r"$2Kp(1-p)=1$")
    ax.scatter(.5, 2, color=PURPLE, zorder=5)
    ax.text(.5, 6, "Perturbations spread\n" + r"$\Lambda>1$", ha="center", color=ORANGE)
    ax.text(.5, .65, "Perturbations shrink\n" + r"$\Lambda<1$", ha="center", color=BLUE)
    ax.annotate("$K_c=2$", (.5, 2), xytext=(.62, 2.7),
                arrowprops={"arrowstyle": "->", "color": MUTED})
    ax.set(xlim=(0, 1), ylim=(0, 8), xlabel="Output bias $p$", ylabel="Input degree $K$")
    ax.legend(loc="upper right")
    save(fig, "08_boolean_phase_diagram")


def feedforward():
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), layout="constrained")
    ax = axes[0]
    title(ax, "A  |  Coherent feed-forward circuit")
    ax.axis("off")
    ax.set(xlim=(-.6, 3.2), ylim=(-.2, 3.5))
    node(ax, (0, 1.2), "$u$", BLUE, .28)
    node(ax, (1.1, 2.7), "$Y$", GREEN, .28)
    node(ax, (2.5, 1.2), "$Z$", ORANGE, .28)
    arrow(ax, (.2, 1.46), (.88, 2.44), BLUE)
    arrow(ax, (1.34, 2.5), (2.27, 1.47), GREEN)
    arrow(ax, (.32, 1.2), (2.18, 1.2), BLUE)
    ax.text(1.8, 1.85, "AND", color=ORANGE, fontsize=11)
    ax.text(1.25, .35, "Direct input AND accumulated $Y$", ha="center", fontsize=11)
    yc = .6
    ton = -np.log(1-yc)
    for ax, duration, heading in zip(axes[1:], [.5, 2.5],
                                     ["B  |  Short pulse is rejected", "C  |  Long pulse is passed"]):
        title(ax, heading)
        t = np.linspace(0, 5, 1600)
        u = (t < duration).astype(float)
        y = np.where(t <= duration, 1-np.exp(-t), (1-np.exp(-duration))*np.exp(-(t-duration)))
        gate = ((u > .5) & (y > yc)).astype(float)
        ax.plot(t, u, color=BLUE, ls="--", lw=1.5, label="Input $u$")
        ax.plot(t, y, color=GREEN, label="$y(t)$")
        ax.plot(t, gate, color=ORANGE, label="Production gate", lw=2.2)
        ax.axhline(yc, color=MUTED, ls=":", lw=1.2)
        ax.text(4.9, yc+.04, "$y_c=0.6$", ha="right", fontsize=10)
        if duration > ton:
            ax.axvline(ton, color=MUTED, ls=":", lw=1)
            ax.text(ton+.05, .37, r"$t_{\rm on}\simeq0.916$", fontsize=10)
        ax.text(.98, .35, rf"$T={duration}$", transform=ax.transAxes, ha="right")
        ax.set(xlim=(0, 5), ylim=(-.06, 1.25), xlabel="Time $t$", ylabel="Normalized signal")
        ax.legend(loc="upper right", bbox_to_anchor=(1, 1.01), fontsize=9.5)
    save(fig, "09_feedforward_filter")


if __name__ == "__main__":
    for make in [regulation, self_activation, hysteresis, toggle, oscillator,
                 expression_noise, landscape, boolean_network, feedforward]:
        make()
        print(f"Generated {make.__name__}")
