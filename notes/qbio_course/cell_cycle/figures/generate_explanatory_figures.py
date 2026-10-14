"""Generate five explanatory PNG/SVG figures with NumPy and Matplotlib.

Run from any directory: python notes/qbio_course/cell_cycle/figures/generate_explanatory_figures.py
These are diagrams and analytic model calculations, not experimental data.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path as PlotPath
import numpy as np
from figure_style import configure_style, save_figure


OUT = Path(__file__).resolve().parent
BLUE, ORANGE, GREEN, PURPLE = "#2767a0", "#d47725", "#238579", "#8057a5"
INK, MUTED = "#243747", "#637181"
configure_style(line_width=2.7, text_color=INK)


def save(fig, name):
    save_figure(fig, OUT, name)
    plt.close(fig)


def arrow(ax, start, end, color=MUTED, **kwargs):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>",
                 mutation_scale=16, linewidth=1.8, color=color, **kwargs))


def cell(ax, x, y, color):
    ax.add_patch(FancyBboxPatch((x-.25, y-.14), .5, .28,
                 boxstyle="round,pad=0,rounding_size=.14",
                 facecolor=color, edgecolor=color, linewidth=1.5))


def renewal_and_sampling():
    fig, (ax, density) = plt.subplots(1, 2, figsize=(12, 5),
                            gridspec_kw={"width_ratios": [1.25, 1]},
                            layout="constrained")
    ax.set(xlim=(0, 9.6), ylim=(-.5, 6.5))
    ax.axis("off")
    ax.set_title("A  |  Two possible histories", loc="left", pad=14)

    ax.text(.35, 5.8, "Not yet divided", weight="bold", color=BLUE)
    for x in (.7, 8.6):
        cell(ax, x, 5.0, BLUE)
    arrow(ax, (1.1, 5.0), (8.2, 5.0), BLUE)
    ax.text(4.65, 5.35, r"Probability $S(t)$", ha="center")
    ax.text(8.6, 4.55, "1 cell", ha="center", color=BLUE)
    ax.text(.7, 3.9, "OR", color=MUTED, fontsize=10, weight="bold")

    ax.text(.35, 3.25, r"First division at $\tau$", weight="bold", color=GREEN)
    cell(ax, .7, 1.85, BLUE)
    arrow(ax, (1.1, 1.85), (3.1, 1.85), BLUE)
    ax.text(2.1, 2.15, r"$p(\tau)\,d\tau$", ha="center")
    ax.plot([3.2, 3.7, 3.2], [2.0, 1.85, 1.7], color=GREEN, lw=2)
    for y in (2.55, 1.05):
        arrow(ax, (3.7, 1.85), (4.7, y), GREEN)
        ax.add_patch(FancyBboxPatch((4.85, y-.34), 3.75, .68,
                     boxstyle="round,pad=.06,rounding_size=.14",
                     facecolor="#e9f4f1", edgecolor=GREEN, lw=1.4))
        ax.text(6.72, y, r"Descendants: $M(t-\tau)$",
                ha="center", va="center", fontsize=11, color=GREEN)
    ax.text(6.7, .28, "Two daughter subtrees", ha="center", color=GREEN)

    arrow(ax, (.5, -.13), (9.1, -.13))
    for x, label in ((.7, "0"), (3.7, r"$\tau$"), (8.6, "$t$")):
        ax.plot([x, x], [-.2, -.06], color=MUTED, lw=1)
        ax.text(x, -.42, label, ha="center")

    # Gamma(shape=4, rate=4): mean 1; Euler-Lotka has an explicit root.
    time = np.linspace(0, 2.7, 650)
    growth = 4 * (2**.25 - 1)
    p = 4**4 / 6 * time**3 * np.exp(-4*time)
    q = 2 * np.exp(-growth*time) * p
    density.fill_between(time, p, color=BLUE, alpha=.10)
    density.fill_between(time, q, color=ORANGE, alpha=.10)
    density.plot(time, p, color=BLUE, label=r"Chronological $p(\tau)$")
    density.plot(time, q, color=ORANGE, label=r"Growth-weighted $q(\tau)$")
    density.axvline(1, color=BLUE, ls=":", lw=1.4)
    density.axvline(4/(4+growth), color=ORANGE, ls=":", lw=1.4)
    density.set_title("B  |  Short lifetimes get more weight", loc="left", pad=14)
    density.set(xlim=(0, 2.7), ylim=(0, 1.55),
                xlabel=r"Lifetime $\tau$ (mean of $p$ = 1)", ylabel="Probability density")
    density.text(.95, .95, r"$q/p=2e^{-\Lambda\tau}$",
                 transform=density.transAxes, ha="right", va="top", fontsize=13)
    density.legend(loc="upper right", bbox_to_anchor=(1, .83), fontsize=10)
    density.grid(alpha=.14)
    save(fig, "06_renewal_and_sampling")


def accumulated_time_covariances():
    fig, (matrix, variance) = plt.subplots(1, 2, figsize=(12, 5),
                              gridspec_kw={"width_ratios": [1, 1.25]},
                              layout="constrained")
    n_matrix = 8
    lag = np.abs(np.arange(n_matrix)[:, None] - np.arange(n_matrix)[None, :])
    covariance = (-.5)**lag
    cmap = LinearSegmentedColormap.from_list("signed_covariance", [BLUE, "white", ORANGE])
    image = matrix.imshow(covariance, cmap=cmap, vmin=-1, vmax=1)
    matrix.set_title("A  |  Every square enters the variance", loc="left", pad=14)
    matrix.set(xticks=np.arange(n_matrix), yticks=np.arange(n_matrix),
               xticklabels=np.arange(1, n_matrix+1), yticklabels=np.arange(1, n_matrix+1),
               xlabel=r"Generation $j$", ylabel=r"Generation $i$")
    for j in range(n_matrix):
        matrix.text(j, j, "1", ha="center", va="center", color="white", fontsize=11)
    for x, y in ((2, 3), (3, 2)):
        matrix.add_patch(Circle((x, y), .42, fill=False, edgecolor=PURPLE, lw=2.2))
    fig.colorbar(image, ax=matrix, orientation="horizontal", shrink=.88, pad=.13,
                 label=r"$\mathrm{Cov}(\tau_i,\tau_j)/\sigma_\tau^2$", ticks=[-1, 0, 1])

    ns = np.arange(1, 41)
    for r, color, label in ((-.5, BLUE, r"Alternating memory: $r=-1/2$"),
                             (0., GREEN, r"Independent: $r=0$"),
                             (.5, ORANGE, r"Persistent memory: $r=1/2$")):
        values = np.array([n + 2*np.sum((n-np.arange(1, n))*r**np.arange(1, n))
                           for n in ns])
        variance.plot(ns, values, color=color, label=label)
    variance.set_title("B  |  Same variance per generation", loc="left", pad=14)
    variance.set(xlim=(0, 40), ylim=(0, 120), xlabel=r"Number of generations $n$",
                 ylabel=r"Accumulated variance $\mathrm{Var}(T_n)/\sigma_\tau^2$")
    variance.legend(loc="upper left", fontsize=10)
    variance.grid(alpha=.14)
    save(fig, "07_accumulated_time_covariances")


def discounted_branching():
    fig, ax = plt.subplots(figsize=(12, 4.7), layout="constrained")
    ax.set(xlim=(0, 10.8), ylim=(-.7, 4.6))
    ax.axis("off")
    ax.set_title("Each arrow carries daughters × inheritance × time discount", pad=14)
    positions = [(2.1, BLUE, r"Short state", r"$\tau_s=\mu-\delta$"),
                 (8.7, ORANGE, "Long state", r"$\tau_\ell=\mu+\delta$")]
    for x, color, label, lifetime in positions:
        ax.add_patch(Circle((x, 1.9), .8, facecolor=color, edgecolor=color, alpha=.10))
        ax.add_patch(Circle((x, 1.9), .8, fill=False, edgecolor=color, lw=2))
        ax.text(x, 2.12, label, ha="center", va="center", color=color, weight="bold")
        ax.text(x, 1.67, lifetime, ha="center", va="center", color=color, fontsize=14)
        loop = PlotPath([(x+.43, 2.56), (x+1.45, 4.10), (x-1.45, 4.10), (x-.43, 2.56)],
                        [PlotPath.MOVETO, PlotPath.CURVE4, PlotPath.CURVE4, PlotPath.CURVE4])
        ax.add_patch(FancyArrowPatch(path=loop, arrowstyle="-|>", mutation_scale=18,
                                    color=color, lw=2))
    ax.text(2.1, 4.08, r"$2p\,e^{-s\tau_s}$", ha="center", color=BLUE, fontsize=15)
    ax.text(8.7, 4.08, r"$2p\,e^{-s\tau_\ell}$", ha="center", color=ORANGE, fontsize=15)
    arrow(ax, (2.9, 2.02), (7.9, 2.02), BLUE, connectionstyle="arc3,rad=-.24")
    arrow(ax, (7.9, 1.76), (2.9, 1.76), ORANGE, connectionstyle="arc3,rad=-.24")
    ax.text(5.4, 3.09, r"$2(1-p)\,e^{-s\tau_s}$", ha="center", color=BLUE, fontsize=15)
    ax.text(5.4, .62, r"$2(1-p)\,e^{-s\tau_\ell}$", ha="center", color=ORANGE, fontsize=15)
    ax.text(5.4, -.18, r"$K_s(i,j)=2e^{-s\tau_i}P_{ij}$",
            ha="center", fontsize=16)
    save(fig, "08_discounted_branching")


def telescoping_identity():
    fig, ax = plt.subplots(figsize=(12, 4.5), layout="constrained")
    ax.set(xlim=(0, 12.2), ylim=(0, 4.9))
    ax.axis("off")
    ax.set_title("Adjacent cycles share a size term with opposite signs", pad=14)
    left_positions = [.4, 4.3, 8.2]
    for j, left in enumerate(left_positions):
        center = left+1.8
        ax.text(center, 4.35, f"Cycle {j}", ha="center", weight="bold")
        ax.text(center, 3.79, rf"$\lambda\tau_{j}-\ln2$", ha="center", fontsize=15)
        ax.text(center, 3.50, "=", ha="center", va="center", fontsize=14)
        ax.add_patch(FancyBboxPatch((left, 2.32), 3.6, 1.03,
                     boxstyle="round,pad=.05,rounding_size=.13",
                     facecolor="#f7f9fb", edgecolor="#cdd6df", lw=1.2))
    terms = [(1.3, r"$-\ln v_0$", ORANGE, False),
             (3.1, r"$+\ln v_1$", BLUE, True),
             (5.2, r"$-\ln v_1$", BLUE, True),
             (7.0, r"$+\ln v_2$", PURPLE, True),
             (9.1, r"$-\ln v_2$", PURPLE, True),
             (10.9, r"$+\ln v_3$", GREEN, False)]
    for x, label, color, cancels in terms:
        ax.text(x, 2.87, label, ha="center", va="center", color=color, fontsize=16)
        if cancels:
            ax.plot([x-.55, x+.55], [2.72, 3.02], color=color, alpha=.60, lw=1.7)
    for start, end, color in ((3.1, 5.2, BLUE), (7., 9.1, PURPLE)):
        curve = PlotPath([(start, 2.45), (start, 1.63), (end, 1.63), (end, 2.45)],
                         [PlotPath.MOVETO, PlotPath.CURVE4, PlotPath.CURVE4, PlotPath.CURVE4])
        ax.add_patch(FancyArrowPatch(path=curve, arrowstyle="-", color=color, lw=1.8))
        ax.text((start+end)/2, 1.69, "cancel", ha="center", va="top", fontsize=11, color=color)
    ax.text(6.1, 1.06, "Only the endpoint sizes remain — for any number of cycles",
            ha="center", weight="bold", fontsize=12)
    ax.text(6.1, .41,
            r"$T_n=\dfrac{n\ln2}{\lambda}+\dfrac{\ln v_n-\ln v_0}{\lambda}$",
            ha="center", va="center", fontsize=17)
    save(fig, "09_telescoping_identity")


def histories_not_histograms():
    fig, (history, total) = plt.subplots(1, 2, figsize=(12, 4.7),
                            gridspec_kw={"width_ratios": [1.1, 1]},
                            layout="constrained")
    history.set(xlim=(0, 6.2), ylim=(-.3, 5.8))
    history.axis("off")
    history.set_title("A  |  Same single-cycle distribution", loc="left", pad=14)
    for x, label in ((1.2, "Two-cycle history"), (3.45, "Independent"), (5.25, "Alternating")):
        history.text(x, 5.0, label, ha="center", weight="bold", fontsize=11)
    for j, pair in enumerate((("s", "s"), ("s", "l"), ("l", "s"), ("l", "l"))):
        y = 4.15-j
        for column, state in enumerate(pair):
            x, color = .64+.87*column, BLUE if state == "s" else ORANGE
            history.add_patch(FancyBboxPatch((x-.30, y-.21), .60, .42,
                              boxstyle="round,pad=0,rounding_size=.21",
                              edgecolor=color, facecolor=color, alpha=.14))
            history.text(x, y, r"$\tau_s$" if state == "s" else r"$\tau_\ell$",
                         color=color, ha="center", va="center", fontsize=12)
        history.text(3.45, y, "1/4", ha="center", va="center", fontsize=14)
        history.text(5.25, y, "1/2" if pair[0] != pair[1] else "0",
                     ha="center", va="center", fontsize=14,
                     color=GREEN if pair[0] != pair[1] else MUTED)
    history.text(3.05, .22, r"Both have $P(\tau_s)=P(\tau_\ell)=1/2$",
                 ha="center", fontsize=12)

    x = np.array([1., 2., 3.])
    probabilities = np.array([.25, .5, .25])
    total.vlines(x, 0, probabilities, color=BLUE, linewidth=2.3)
    total.scatter(x, probabilities, color=BLUE, s=65, zorder=4, label="Independent")
    total.vlines(2, 0, 1, color=GREEN, linewidth=2.0)
    total.scatter([2], [1], color=GREEN, marker="s", s=65, zorder=4, label="Strict alternation")
    total.set_title(r"B  |  Distribution of $T_2=\tau_0+\tau_1$", loc="left", pad=14)
    total.set(xlim=(.55, 3.45), ylim=(0, 1.14),
              xticks=x, xticklabels=[r"$2\tau_s$", r"$\tau_s+\tau_\ell$", r"$2\tau_\ell$"],
              xlabel="Total duration of two generations", ylabel="Probability")
    total.legend(loc="upper right", fontsize=10)
    total.grid(alpha=.14)
    save(fig, "10_histories_and_total_time")


if __name__ == "__main__":
    renewal_and_sampling()
    accumulated_time_covariances()
    discounted_branching()
    telescoping_identity()
    histories_not_histograms()
    print("Generated five explanatory figures in", OUT)
