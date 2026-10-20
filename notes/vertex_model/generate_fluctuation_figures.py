"""Reproducible calculations for fluctuations/body.md; no experimental fits."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch

ROOT = Path(__file__).resolve().parent / "fluctuations" / "figures"
BLUE, ORANGE, GREEN, INK = "#277da8", "#d67837", "#429582", "#263e50"
for filename in ["Arial.ttf", "Arial Bold.ttf"]:
    path = Path("/System/Library/Fonts/Supplemental") / filename
    if path.exists():
        font_manager.fontManager.addfont(path)
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14,
    "legend.fontsize": 12, "mathtext.fontset": "cm", "svg.fonttype": "path",
    "axes.spines.top": False, "axes.spines.right": False,
    "text.color": INK, "axes.labelcolor": INK, "axes.titlecolor": INK,
    "lines.linewidth": 2.2, "svg.hashsalt": "vertex-fluctuations",
    "savefig.facecolor": "white",
})


def save(fig, name):
    ROOT.mkdir(parents=True, exist_ok=True)
    fig.savefig(ROOT / (name + ".svg"), bbox_inches="tight", pad_inches=.2,
                metadata={"Date": None})
    fig.savefig(ROOT / (name + ".png"), bbox_inches="tight", pad_inches=.2, dpi=300)
    plt.close(fig)


def correlation_scales():
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 4.1), layout="constrained")
    rng = np.random.default_rng(20261020)
    dt = .01
    time = np.arange(0, 15 + dt / 2, dt)
    lag = np.linspace(0, 6, 500)
    omega = np.geomspace(.02, 100, 500)
    for tau, color in [(.2, BLUE), (2, ORANGE)]:
        a = np.exp(-dt / tau)
        x = np.empty(time.size)
        x[0] = rng.normal()
        noise = rng.normal(size=time.size - 1)
        for i in range(1, time.size):
            x[i] = a * x[i - 1] + np.sqrt(1 - a * a) * noise[i - 1]
        axs[0].plot(time, x, color=color, lw=1, alpha=.85, label=rf"$\tau={tau:g}$")
        axs[1].plot(lag, np.exp(-lag / tau), color=color)
        axs[2].loglog(omega, 2 * tau / (1 + omega**2 * tau**2), color=color)
    axs[0].axhline(0, color=INK, lw=.8, alpha=.4)
    axs[0].set(title="A  Sample paths, same variance", xlabel="Time", ylabel=r"$x(t)$")
    axs[0].legend(loc="upper right")
    axs[1].set(title="B  Exact autocorrelation", xlabel="Lag time", ylabel=r"$C_x(t)$",
               ylim=(-.03, 1.04))
    axs[2].set(title="C  Exact power spectrum", xlabel=r"Frequency, $\omega$",
               ylabel=r"$S_x(\omega)$")
    save(fig, "correlation_scales")


def spectrum_response():
    fig, axs = plt.subplots(2, 2, figsize=(10.8, 7.3), layout="constrained")
    omega = np.geomspace(.01, 100, 600)
    gain = 1 / (1 + omega**2)
    axs[0, 0].loglog(omega, gain, color=INK)
    axs[0, 0].set(title="A  Mechanical filtering", xlabel=r"$\omega$",
                  ylabel=r"$|\chi(\omega)|^2$")
    for tau, color in [(.1, BLUE), (1, ORANGE), (10, GREEN)]:
        spectrum = 2 * tau / (1 + omega**2 * tau**2)
        axs[0, 1].loglog(omega, gain * spectrum, color=color, label=rf"$\tau={tau:g}$")
        axs[1, 1].loglog(omega, spectrum / 2, color=color, label=rf"$\tau={tau:g}$")
    axs[0, 1].set(title="B  Displacement spectrum", xlabel=r"$\omega$",
                  ylabel=r"$S_u(\omega)$")
    axs[0, 1].legend()
    tau = np.geomspace(.01, 100, 500)
    axs[1, 0].semilogx(tau, tau / (1 + tau), color=BLUE, label=r"Fixed $s_f^2=1$")
    axs[1, 0].semilogx(tau, 1 / (1 + tau), color=ORANGE, label=r"Fixed $D_f=1$")
    axs[1, 0].set(title="C  What is held fixed?", xlabel=r"$\tau$", ylabel=r"$\langle u^2\rangle$",
                  ylim=(-.04, 1.04))
    axs[1, 0].legend(loc="center right")
    axs[1, 1].set(title="D  Fluctuation / response ratio", xlabel=r"$\omega$",
                  ylabel=r"$\mathcal{E}_{\rm FR}(\omega)$")
    axs[1, 1].legend()
    save(fig, "spectrum_response")


def spring(ax, start, end, y):
    x = np.linspace(start, end, 25)
    z = np.zeros(25)
    z[2:-2] = .075 * (-1.)**np.arange(21)
    ax.plot(x, y + z, color=INK, lw=1.2)


def correlated_forces():
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.8), layout="constrained")
    ax = axs[0]
    for y, c, color in [(.8, 1, BLUE), (-.8, -1, ORANGE)]:
        spring(ax, -1.8, -.6, y)
        spring(ax, -.6, .6, y)
        spring(ax, .6, 1.8, y)
        ax.scatter([-.6, .6], [y, y], s=105, color=INK, zorder=3)
        for x, direction in [(-.6, c), (.6, 1)]:
            ax.add_patch(FancyArrowPatch((x, y + .23), (x + direction * .48, y + .23),
                                        arrowstyle="-|>", mutation_scale=16,
                                        color=color, lw=2.3))
        ax.text(0, y + .65, "Same force: c = +1" if c == 1 else "Opposite force: c = -1",
                ha="center", fontsize=14)
        for x, label in [(-1.2, r"$k_0$"), (0, r"$k_c$"), (1.2, r"$k_0$")]:
            ax.text(x, y - .28, label, ha="center")
        ax.plot([-1.8, -1.8], [y - .15, y + .15], color=INK)
        ax.plot([1.8, 1.8], [y - .15, y + .15], color=INK)
    ax.set(xlim=(-2.05, 2.05), ylim=(-1.4, 1.65), title="A  Two correlated driving forces")
    ax.axis("off")
    c = np.linspace(-1, 1, 401)
    axs[1].plot(c, (1 + c) / 2, color=BLUE, label=r"Common: $\langle q_+^2\rangle$")
    axs[1].plot(c, (1 - c) / 12, color=ORANGE, label=r"Relative: $\langle q_-^2\rangle$")
    axs[1].axvline(0, color=INK, alpha=.3, lw=1)
    axs[1].set(title="B  Mode variances", xlabel=r"Force correlation, $c$",
               ylabel="Displacement variance", ylim=(-.025, 1.04))
    axs[1].legend(loc="upper left")
    save(fig, "correlated_forces")


def window_variances(radius):
    """2D disk-window variance for S=1, 1+4 exp(-q²/2), 1-exp(-q²/2).

    rho=xi=1. Integrate the disk overlap against r exp(-r²/2).
    The overlap loss is evaluated directly to avoid cancellation at large R.
    """
    nodes, weights = np.polynomial.legendre.leggauss(180)
    area = np.pi * radius**2
    loss_integral = np.empty_like(radius)
    for i, rwin in enumerate(radius):
        limit = min(2 * rwin, 12.)
        r = (nodes + 1) * limit / 2
        z = r / (2 * rwin)
        loss = 2 * rwin**2 * np.arcsin(z) + r * rwin * np.sqrt(1 - z * z)
        loss_integral[i] = limit / 2 * np.sum(weights * r * np.exp(-r * r / 2) * loss)
    suppressed = area * np.exp(-2 * radius**2) + loss_integral
    overlap_integral = area - suppressed
    enhanced = area + 4 * overlap_integral
    return area, enhanced, suppressed


def density_windows():
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 4.3), layout="constrained")
    q = np.geomspace(.003, 20, 500)
    specs = [np.ones_like(q), 1 + 4 * np.exp(-q*q/2), -np.expm1(-q*q/2)]
    labels = ["Poisson", "Enhanced, finite S(0)", "Suppressed, S(0) = 0"]
    colors = [INK, ORANGE, BLUE]
    radius = np.geomspace(.03, 200, 350)
    variances = window_variances(radius)
    mean = np.pi * radius**2
    for spec, var, label, color in zip(specs, variances, labels, colors):
        axs[0].loglog(q, spec, color=color, label=label)
        axs[1].semilogx(radius, var / mean, color=color)
        axs[2].loglog(mean, np.sqrt(var), color=color)
    axs[0].set(title="A  Specified structure factors", xlabel=r"Wave number, $q$",
               ylabel=r"$S(q)$", ylim=(2e-6, 7))
    axs[0].legend(loc="lower right", fontsize=10.5)
    axs[1].set(title="B  Disk-window count variance", xlabel=r"Window radius, $R$",
               ylabel=r"$\mathrm{Var}(N_R)/\langle N_R\rangle$", ylim=(-.12, 5.25))
    axs[2].set(title="C  Standard deviation", xlabel=r"Mean count, $\langle N_R\rangle$",
               ylabel=r"$\Delta N_R$")
    save(fig, "density_windows")


if __name__ == "__main__":
    correlation_scales()
    spectrum_response()
    correlated_forces()
    density_windows()
    print("Generated four fluctuation figures (SVG and PNG).")
