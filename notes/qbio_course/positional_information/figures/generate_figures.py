"""Reproduce explanatory figures for the positional-information note.

These are model calculations, not embryo data or a reproduction of the full
network optimization. Requires NumPy, SciPy and Matplotlib. Run this file from
any directory; SVG and PNG outputs are written beside it.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.integrate import simpson
from scipy.special import logsumexp


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
    "legend.frameon": False, "lines.linewidth": 2.2,
    "svg.fonttype": "path", "svg.hashsalt": "positional-information",
})


def save(fig, name):
    target = OUT / f"{name}.svg"
    fig.savefig(target, bbox_inches="tight", metadata={"Date": None})
    target.write_text("\n".join(s.rstrip() for s in target.read_text().splitlines()) + "\n")
    fig.savefig(OUT / f"{name}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def title(ax, text):
    ax.set_title(text, loc="left", pad=13)


def gradient():
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.2), layout="constrained")
    x, ell = np.linspace(0, 1, 500), .2
    c = np.exp(-x / ell)
    ax = axes[0]
    title(ax, "A  |  A gradient and a threshold")
    ax.plot(x, c, color=BLUE, label=r"$c_0 e^{-x/\lambda}$")
    ax.plot(x, 1.1*c, color=ORANGE, label=r"$1.1\,c_0 e^{-x/\lambda}$")
    threshold = np.exp(-.4/ell)
    ax.axhline(threshold, color=MUTED, ls=":", lw=1.1)
    for boundary, color in [(.4, BLUE), (.4+ell*np.log(1.1), ORANGE)]:
        ax.plot([boundary, boundary], [0, threshold], color=color, ls="--", lw=1.3)
    ax.set(xlabel="$x/L$", ylabel="$c/c_0$", xlim=(0, 1), ylim=(0, 1.13))
    ax.legend()
    ax = axes[1]
    title(ax, "B  |  Finite-interval boundary")
    finite = np.cosh((1-x)/ell)/np.sinh(1/ell)
    ax.semilogy(x, c, color=BLUE, label="Semi-infinite approximation")
    ax.semilogy(x, finite, color=GREEN, ls="--", label="Finite interval, zero flux at $L$")
    ax.set(xlabel="$x/L$", ylabel=r"$c/(R\tau/\lambda)$", xlim=(0, 1))
    ax.legend(loc="upper right", fontsize=9)
    ax = axes[2]
    title(ax, "C  |  Poisson counting limit")
    for m0, color in [(400, ORANGE), (1600, BLUE), (6400, GREEN)]:
        error = ell/np.sqrt(m0*c)
        ax.plot(x, 100*error, color=color, label=f"$m_0={m0}$")
    ax.axhline(1, color=MUTED, ls=":", lw=1.2, label="1% EL")
    ax.set(xlabel="$x/L$", ylabel=r"$100\,\sigma_x/L$ (%)", xlim=(0, 1), ylim=(0, 13))
    ax.legend()
    for ax in axes:
        ax.grid(alpha=.12)
    save(fig, "01_gradient_and_counting")


def correlated_noise():
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.1), layout="constrained")
    for ax, tangent, heading in zip(axes[:2], [(1, 1), (1, -1)],
                                    ["A  |  Equal slopes", "B  |  Opposite slopes"]):
        title(ax, heading)
        cov = np.array([[1., .8], [.8, 1.]])*.025
        eigen, vectors = np.linalg.eigh(cov)
        angle = np.degrees(np.arctan2(vectors[1, 1], vectors[0, 1]))
        for center in [(.35, .35), (.65, .65) if tangent[1] > 0 else (.65, .05)]:
            ax.add_patch(Ellipse(center, 2*np.sqrt(eigen[1]), 2*np.sqrt(eigen[0]),
                                 angle=angle, facecolor=BLUE, alpha=.17, edgecolor=BLUE))
        start = np.array([.35, .35])
        end = start + .3*np.array(tangent)
        ax.annotate("", xy=end, xytext=start,
                    arrowprops=dict(arrowstyle="->", color=ORANGE, lw=2.5))
        ax.text(.04, .91, r"$\rho=0.8$; ellipses: 1 SD", transform=ax.transAxes)
        ax.set(xlabel="$g_1$", ylabel="$g_2$", xlim=(0, 1), ylim=(-.2, .9), aspect="equal")
        ax.grid(alpha=.12)
    ax = axes[2]
    title(ax, "C  |  Correlation depends on direction")
    rho = np.linspace(-.9, .9, 400)
    ax.plot(rho, 2/(1+rho), color=BLUE, label=r"$\mu'=(a,a)$")
    ax.plot(rho, 2/(1-rho), color=ORANGE, label=r"$\mu'=(a,-a)$")
    ax.axhline(2, color=MUTED, ls=":", lw=1.1)
    ax.set(xlabel=r"Noise correlation $\rho$", ylabel=r"$\mathcal{J}/(a^2/\sigma^2)$", ylim=(0, 12))
    ax.legend()
    ax.grid(alpha=.12)
    save(fig, "02_noise_geometry")


def decoding_map(mean, sigma=.06, order=13):
    """Noise-averaged Bayes posterior on equal-probability spatial bins.

    Columns are true positions and rows are inferred positions. A tensor-product
    Gauss-Hermite quadrature averages over independent Gaussian expression noise.
    Returned values are densities per unit x/L, not bin probabilities.
    """
    n, d = mean.shape
    nodes, weights = hermgauss(order)
    nodes, weights = np.sqrt(2)*nodes, weights/np.sqrt(np.pi)
    result = np.zeros((n, n))
    for index in np.ndindex(*(order,)*d):
        shift = np.array([nodes[k] for k in index])*sigma
        weight = np.prod([weights[k] for k in index])
        observations = mean + shift
        exponent = -np.sum((mean[:, None, :]-observations[None, :, :])**2, axis=2)/(2*sigma**2)
        posterior = np.exp(exponent-logsumexp(exponent, axis=0, keepdims=True))
        result += weight*posterior
    assert np.allclose(result.sum(axis=0), 1, atol=1e-10)
    return n*result


def ambiguity():
    n = 180
    x = (np.arange(n)+.5)/n
    hump = np.exp(-.5*((x-.5)/.19)**2)
    mean = np.column_stack([hump, x])
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.1), layout="constrained")
    ax = axes[0]
    title(ax, "A  |  A second gene labels the branch")
    ax.plot(x, hump, color=BLUE, label=r"$\mu_1(x)$: symmetric peak")
    ax.plot(x, x, color=ORANGE, label=r"$\mu_2(x)=x/L$")
    ax.axvline(.3, ls=":", color=MUTED, lw=1)
    ax.axvline(.7, ls=":", color=MUTED, lw=1)
    ax.set(xlabel="$x/L$", ylabel="Mean expression", xlim=(0, 1), ylim=(0, 1.05))
    ax.legend(loc="upper left", fontsize=9)
    maps = [decoding_map(mean[:, :1]), decoding_map(mean)]
    for ax, matrix, heading in zip(axes[1:], maps,
                                  ["B  |  One gene: two branches", "C  |  Two genes: a single ridge"]):
        title(ax, heading)
        im = ax.imshow(matrix, origin="lower", extent=(0, 1, 0, 1),
                       aspect="equal", cmap="Blues", vmin=0, vmax=12)
        ax.plot([0, 1], [0, 1], color=ORANGE, ls=":", lw=1)
        ax.set(xlabel="True position $x/L$", ylabel="Inferred position $x^*/L$")
    fig.colorbar(im, ax=list(axes[1:]), shrink=.77, label="Posterior density")
    save(fig, "03_global_ambiguity")


def scalar_information(mean, sigma, x):
    y = np.linspace(-7*sigma, 1+7*sigma, 1801)
    likelihood = np.exp(-.5*((y[:, None]-mean[None, :])/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
    marginal = simpson(likelihood, x=x, axis=1)
    assert abs(simpson(marginal, x=y)-1) < 1e-7
    entropy = -simpson(marginal*np.log2(np.maximum(marginal, 1e-300)), x=y)
    return entropy-.5*np.log2(2*np.pi*np.e*sigma**2)


def optimal_gradient():
    x = np.linspace(0, 1, 801)
    steep = (np.tanh(7*(x-.5))+np.tanh(3.5))/(2*np.tanh(3.5))
    slope = 7/(2*np.tanh(3.5)*np.cosh(7*(x-.5))**2)
    sigma = .025
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.2), layout="constrained")
    ax = axes[0]
    title(ax, "A  |  Same expression range")
    for mean, color, label in [(x, BLUE, "Linear"), (steep, ORANGE, "Steep sigmoid")]:
        info = scalar_information(mean, sigma, x)
        ax.plot(x, mean, color=color, label=f"{label}: {info:.2f} bits")
    ax.set(xlabel="$x/L$", ylabel=r"$\mu(x)$", xlim=(0, 1), ylim=(0, 1.05))
    ax.legend(loc="upper left", fontsize=9)
    ax = axes[1]
    title(ax, "B  |  Where precision is spent")
    ax.semilogy(x, np.full_like(x, sigma), color=BLUE, label="Linear")
    ax.semilogy(x, sigma/slope, color=ORANGE, label="Sigmoid (local estimate)")
    ax.set(xlabel="$x/L$", ylabel=r"$\sigma_x/L=\sigma_g/|\mu'|$", xlim=(0, 1))
    ax.legend()
    ax = axes[2]
    title(ax, "C  |  Small-noise approximation")
    noise = np.geomspace(.009, .2, 18)
    exact = np.array([scalar_information(x, s, x) for s in noise])
    approximation = -np.log2(np.sqrt(2*np.pi*np.e)*noise)
    ax.semilogx(noise, exact, color=BLUE, label="Numerical mutual information")
    ax.semilogx(noise, approximation, color=MUTED, ls="--", label="Small-noise formula")
    ax.set(xlabel=r"Expression noise $\sigma_g$", ylabel="Information (bits)")
    ax.legend(fontsize=9)
    for ax in axes:
        ax.grid(alpha=.12)
    save(fig, "04_optimal_gradient")


def diffusion_tradeoff():
    """A periodic linear lattice with additive, independent source noise.

    dg_i/dt = s_i - g_i + D(g_{i+1}+g_{i-1}-2g_i) + sqrt(Q) xi_i.
    Q and gamma=1 are fixed; diffusion is deterministic coupling, with no
    additional hopping noise. This is an excess-source-noise model.
    """
    n, mode = 64, 4
    eigenvalues = 4*np.sin(np.pi*np.arange(n)/n)**2
    d = np.geomspace(.001, 1000, 500)
    variance = np.mean(1/(1+d[:, None]*eigenvalues), axis=1)
    amplitude = 1/(1+d*eigenvalues[mode])
    relative_j = amplitude**2/variance
    best = d[np.argmax(relative_j)]
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.2), layout="constrained")
    title(axes[0], "A  |  Signal and local noise")
    axes[0].loglog(d, amplitude, color=BLUE, label="Mean mode amplitude")
    axes[0].loglog(d, np.sqrt(variance), color=ORANGE, label="Local noise SD")
    axes[0].set(xlabel=r"$D/(\gamma a^2)$", ylabel="Relative to zero diffusion")
    axes[0].legend()
    title(axes[1], "B  |  A finite precision optimum")
    axes[1].semilogx(d, relative_j, color=GREEN)
    axes[1].axvline(best, color=MUTED, ls=":", lw=1.3)
    axes[1].text(.04, .93, f"Optimum near {best:.2f}", transform=axes[1].transAxes, va="top")
    axes[1].set(xlabel=r"$D/(\gamma a^2)$", ylabel=r"$\mathcal{J}_{\rm slope}(D)/\mathcal{J}_{\rm slope}(0)$")
    title(axes[2], "C  |  Strong coupling erases structure")
    x = np.linspace(0, 1, 500)
    for dv, color in [(0, BLUE), (best, GREEN), (100, ORANGE)]:
        amp = .35/(1+dv*eigenvalues[mode])
        axes[2].plot(x, 1+amp*np.cos(2*np.pi*mode*x), color=color, label=rf"$D/(\gamma a^2)={dv:.2g}$")
    axes[2].set(xlabel="Position / system length", ylabel="Mean expression", xlim=(0, 1), ylim=(.58, 1.42))
    axes[2].legend(fontsize=9)
    for ax in axes:
        ax.grid(alpha=.12)
    save(fig, "05_diffusion_tradeoff")


def mean_and_worst():
    n = 1000
    x = (np.arange(n)+.5)/n
    uniform = np.full(n, .015)
    defective = np.full(n, .015*(.015/.06)**(.1/.9))
    defective[(x >= .45) & (x < .55)] = .06
    geometric = lambda s: np.exp(np.mean(np.log(s)))
    assert np.isclose(geometric(uniform), geometric(defective), rtol=1e-12)
    fig, axes = plt.subplots(1, 2, figsize=(11.7, 4.2), layout="constrained")
    title(axes[0], "A  |  Equal geometric mean error")
    axes[0].plot(x, 100*uniform, color=BLUE, label="Uniform precision")
    axes[0].plot(x, 100*defective, color=ORANGE, label="A local defect")
    axes[0].set(xlabel="$x/W$", ylabel=r"$100\,\sigma_x/W$ (%)", xlim=(0, 1), ylim=(0, 7))
    axes[0].legend()
    title(axes[1], "B  |  Average and worst-case objectives")
    categories = ["Geometric mean", "RMS error", "Worst position"]
    for offset, arr, color, label in [(-.18, uniform, BLUE, "Uniform"), (.18, defective, ORANGE, "Local defect")]:
        vals = [geometric(arr), np.sqrt(np.mean(arr**2)), np.max(arr)]
        axes[1].bar(np.arange(3)+offset, 100*np.array(vals), width=.34, color=color, label=label)
    axes[1].set(xticks=np.arange(3), xticklabels=categories, ylabel="Error (% of interval)", ylim=(0, 7))
    axes[1].legend()
    axes[1].tick_params(axis="x", labelsize=10)
    for ax in axes:
        ax.grid(axis="y", alpha=.12)
    save(fig, "06_average_and_defects")


if __name__ == "__main__":
    for draw in [gradient, correlated_noise, ambiguity, optimal_gradient,
                 diffusion_tradeoff, mean_and_worst]:
        draw()
        print(draw.__name__)
