"""Reproduce explanatory figures for the positional-information note.

These are model calculations, not embryo data or a reproduction of the full
network optimization. Requires NumPy, SciPy and Matplotlib. Run this file from
any directory; SVG and PNG outputs are written beside it.
"""

from pathlib import Path
import sys
from textwrap import fill

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, Rectangle
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.integrate import simpson
from scipy.special import logsumexp, ndtr


OUT = Path(__file__).resolve().parent
BLUE, ORANGE, GREEN, PURPLE = "#2767a0", "#d47725", "#238579", "#8057a5"
INK, MUTED = "#243747", "#637181"
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from plot_style import configure_style, save_figure

configure_style(line_width=2.2, text_color=INK)


def save(fig, name):
    save_figure(fig, OUT, name)
    plt.close(fig)


def title(ax, text):
    ax.set_title(fill(text, width=34, break_long_words=False), loc="left", pad=16)


def french_flag():
    """Two fixed thresholds, then independent Gaussian input-reading noise."""
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.9), layout="constrained")
    z = np.linspace(0, 1, 700)
    c = np.exp(-z/.4)
    high, low = np.exp(-np.array([1/3, 2/3])/.4)
    colors = ["#2767a0", "#a8b2bc", "#d45150"]
    ax = axes[0]
    title(ax, "A  |  One input, two thresholds")
    ax.plot(z, c, color=INK, label=r"$c/c_0=e^{-z/0.4}$")
    for value, boundary, label in [(high, 1/3, r"$c_H$"), (low, 2/3, r"$c_L$")]:
        ax.axhline(value, color=MUTED, ls=":", lw=1.2)
        ax.plot([boundary, boundary], [0, value], color=MUTED, ls="--", lw=1.1)
        ax.text(.87, value+.025, label)
    ax.set(xlabel=r"Relative position $z=x/L$", ylabel=r"Concentration $c/c_0$",
           xlim=(0, 1), ylim=(0, 1.08))
    ax.legend(loc="upper right", fontsize=13)
    ax = axes[1]
    title(ax, "B  |  A coarse fate decision")
    for j, (face, name) in enumerate(zip([colors[0], "white", colors[2]], ["Blue", "White", "Red"])):
        ax.add_patch(Rectangle((j/3, .38), 1/3, .26, facecolor=face,
                               edgecolor=INK, lw=1))
        ax.text((j+.5)/3, .51, name, ha="center", va="center",
                color=INK if j==1 else "white", fontsize=16)
    ax.text(.5, .87, "Noiseless threshold readout", ha="center")
    ax.text(.5, .16, "Same fate within each domain;\nposition inside it remains unresolved.", ha="center", fontsize=14)
    ax.set(xlabel=r"Relative position $z=x/L$", xlim=(0, 1), ylim=(0, 1), yticks=[])
    ax.spines["left"].set_visible(False)
    ax = axes[2]
    title(ax, "C  |  Noise blurs the decisions")
    sigma = .045
    blue = ndtr((c-high)/sigma)
    red = ndtr((low-c)/sigma)
    white = 1-blue-red
    assert np.all(white >= -1e-12)
    for p, color, name in zip([blue, white, red], colors, ["Blue", "White", "Red"]):
        ax.plot(z, p, color=color, label=name)
    ax.text(.04, .51, r"$\sigma_c/c_0=0.045$", transform=ax.transAxes)
    ax.set(xlabel=r"Relative position $z=x/L$", ylabel="Fate probability over repeated reads",
           xlim=(0, 1), ylim=(0, 1.04))
    ax.legend(loc="center right", fontsize=13)
    for ax in (axes[0], axes[2]):
        ax.grid(alpha=.12)
    save(fig, "00_french_flag")


def gradient():
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 5.0), layout="constrained")
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
    ax.legend(loc="upper right", fontsize=12)
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


def molecular_readout():
    """One-site binding cartoon and reproducible two-state reaction histories.

    A free site binds at rate k_on*c; an occupied site releases at k_off.
    Only these two reactions are simulated, not transcription or hb regulation.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.8), layout="constrained")
    ax = axes[0]
    title(ax, "A  |  A DNA site senses local Bcd")
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    ax.text(.5, .93, "More local Bcd means\na shorter wait to bind", ha="center", fontsize=14)
    for center, occupied in [(.23, False), (.77, True)]:
        for y in [.56, .59]:
            ax.plot([center-.17, center+.17], [y, y], color=MUTED, lw=1.6)
        for x in np.linspace(center-.17, center+.17, 15):
            ax.plot([x, x], [.56, .59], color=MUTED, lw=.7)
        ax.add_patch(Rectangle((center-.045, .535), .09, .08,
                               edgecolor=BLUE, facecolor="#e6eef5", lw=1.5, zorder=3))
        if occupied:
            ax.add_patch(Circle((center, .655), .062, color=BLUE, zorder=4))
            ax.text(center, .655, "Bcd", ha="center", va="center", color="white", fontsize=13, zorder=5)
        ax.text(center, .43, "Occupied (1)" if occupied else "Empty (0)", ha="center")
        ax.text(center, .30, "Same DNA\nbinding site", ha="center", va="center", fontsize=13, color=MUTED)
    for x, y in [(.13, .77), (.32, .82), (.38, .70)]:
        ax.add_patch(Circle((x, y), .018, color=BLUE, alpha=.7))
    ax.annotate("", xy=(.60, .64), xytext=(.40, .64),
                arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.6))
    ax.text(.5, .69, "bind", ha="center", fontsize=13, color=BLUE)
    ax.annotate("", xy=(.40, .51), xytext=(.60, .51),
                arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.6))
    ax.text(.5, .45, "release", ha="center", fontsize=13, color=ORANGE)
    ax.text(.5, .06, "Binding states influence transcription\nand thus RNA production", ha="center", va="center", fontsize=15)
    ax = axes[1]
    title(ax, "B  |  Concentration changes the history")
    rng = np.random.default_rng(734)
    duration, off_rate = 120., .5
    for on_rate, offset, color, label in [(.125, 0., BLUE, "Low c: longer empty intervals"),
                                          (.5, 2.3, ORANGE, "High c: shorter empty intervals")]:
        times, states = [0.], [0]
        while times[-1] < duration:
            rate = on_rate if states[-1] == 0 else off_rate
            next_time = times[-1]+rng.exponential(1/rate)
            times.append(min(next_time, duration))
            states.append(1-states[-1] if next_time < duration else states[-1])
        ax.step(times, np.array(states)+offset, where="post", color=color, lw=1.8)
        ax.text(0, offset+1.28, label, color=color, fontsize=14)
        ax.axhline(offset, color=MUTED, alpha=.2, lw=.8)
    ax.set(xlim=(0, duration), ylim=(-.3, 4.1), xlabel="Time (s)",
           ylabel="Site occupancy (0 = empty, 1 = bound)",
           yticks=[0, 1, 2.3, 3.3], yticklabels=["0", "1", "0", "1"])
    ax.grid(axis="x", alpha=.12)
    save(fig, "09_molecular_readout")


def reading_noise():
    """A local linear teaching example, not measured gene-expression data."""
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.2), layout="constrained")
    g = np.linspace(35, 175, 900)
    for ax, sigma, heading in zip(axes[0], [5, 20],
                                  ["A  |  Small reading noise", "B  |  Large reading noise"]):
        title(ax, heading)
        for mean, color, label in [(100, BLUE, "$x=0.40L$"),
                                    (110, ORANGE, "$x=0.41L$")]:
            density = np.exp(-.5*((g-mean)/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
            ax.plot(g, density, color=color, label=label)
            ax.fill_between(g, density, color=color, alpha=.12)
            ax.axvline(mean, color=color, ls=":", lw=1)
        ax.text(.05, .82, rf"$\sigma_g={sigma}$", transform=ax.transAxes, fontsize=16)
        ax.set(xlabel="One expression reading $g$\n(arbitrary units)",
               ylabel="Probability density\nof repeated readings", xlim=(35, 175), ylim=(0, .088))
        ax.legend(loc="upper right")
        ax.grid(alpha=.12)
    z = np.linspace(.37, .43, 400)
    for ax, slope, heading in zip(axes[1], [1000, 250],
                                  ["C  |  Steep spatial response", "D  |  Shallow spatial response"]):
        title(ax, heading)
        sigma_g, z0, mean = 5, .4, 100
        sigma_z = sigma_g/slope
        left, right = z0-sigma_z, z0+sigma_z
        ax.axhspan(mean-sigma_g, mean+sigma_g, color=ORANGE, alpha=.18)
        ax.fill_betweenx([65, mean-sigma_g], left, right, color=GREEN, alpha=.13)
        ax.plot(z, mean+slope*(z-z0), color=BLUE)
        for edge, value in [(left, mean-sigma_g), (right, mean+sigma_g)]:
            ax.plot([edge, edge], [65, value], color=GREEN, ls="--", lw=1.3)
            ax.scatter([edge], [value], color=GREEN, s=24, zorder=4)
        ax.scatter([z0], [mean], color=INK, s=24, zorder=4)
        ax.annotate("", xy=(right, 76), xytext=(left, 76),
                    arrowprops=dict(arrowstyle="<->", color=GREEN, lw=1.4))
        ax.text(z0, 69, rf"Half-width: $\sigma_x/L={sigma_z:g}$", ha="center", color=GREEN, fontsize=14)
        ax.text(.05, .89, r"Same reading noise: $\sigma_g=5$", transform=ax.transAxes, fontsize=14)
        ax.set(xlabel=r"Relative position $x/L$", ylabel="Mean expression $\\mu(x)$\n(arbitrary units)",
               xlim=(.37, .43), ylim=(65, 135), xticks=[.37, .39, .4, .41, .43])
        ax.grid(alpha=.12)
    save(fig, "07_reading_noise")


def information_budget():
    """Four equiprobable positions, readout of half with symmetric errors."""
    def entropy(p):
        p = np.asarray(p, dtype=float)
        return -np.sum(p[p > 0]*np.log2(p[p > 0]))
    prior = np.full(4, .25)
    channel = np.array([[.9, .1], [.9, .1], [.1, .9], [.1, .9]])
    joint = prior[:, None]*channel
    output = joint.sum(axis=0)
    posterior = joint/output
    remaining = sum(output[j]*entropy(posterior[:, j]) for j in range(2))
    gained = entropy(prior)-remaining
    assert np.allclose(posterior[:, 0], [.45, .45, .05, .05])
    assert np.isclose(gained, .5310044064107188)
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 5.0), layout="constrained")
    ax = axes[0]
    title(ax, "A  |  Update after an ON reading")
    ax.bar(np.arange(4)-.18, prior, width=.36, color=BLUE, label="Before reading")
    ax.bar(np.arange(4)+.18, posterior[:, 0], width=.36, color=GREEN, label="After ON")
    ax.set(xticks=np.arange(4), xticklabels=["A", "B", "C", "D"],
           xlabel="Candidate position", ylabel="Probability", ylim=(0, .62))
    ax.legend(loc="upper right", fontsize=13)
    ax = axes[1]
    title(ax, "B  |  Count uncertainty in bits")
    ax.bar(0, 2, color="#b6c0c9", width=.62)
    ax.bar(1, remaining, color="#b6c0c9", width=.62, label="Remaining uncertainty")
    ax.bar(1, gained, bottom=remaining, color=GREEN, width=.62, label="Information gained")
    ax.text(0, 1, "$H(X)=2$", ha="center", va="center")
    ax.text(1, remaining/2, f"$H(X|G)$\n$={remaining:.3f}$", ha="center", va="center")
    ax.text(1, remaining+gained/2, f"$I={gained:.3f}$", ha="center", va="center", color="white", fontsize=14)
    ax.set(xticks=[0, 1], xticklabels=["Before", "After"], ylabel="Bits", ylim=(0, 2.8))
    ax.legend(loc="upper center", fontsize=12.5)
    ax = axes[2]
    title(ax, "C  |  Error rate and information")
    error = np.linspace(0, .5, 400)
    binary_entropy = np.array([entropy([p, 1-p]) for p in error])
    ax.plot(error, 1-binary_entropy, color=BLUE)
    ax.scatter([.1], [gained], color=GREEN, zorder=3)
    ax.annotate("10% error:\n0.531 bits", xy=(.1, gained), xytext=(.21, .68),
                fontsize=14, arrowprops=dict(arrowstyle="->", color=MUTED))
    ax.text(.48, .94, "At most 1 bit:\nonly the half is encoded", ha="right", va="top", fontsize=13)
    ax.set(xlabel=r"Symmetric readout error $\epsilon$", ylabel="$I(X;G)$ (bits)",
           xlim=(0, .5), ylim=(0, 1.06))
    for ax in axes:
        ax.grid(axis="y", alpha=.12)
    save(fig, "08_information_budget")


def correlated_noise():
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.9), layout="constrained")
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
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.9), layout="constrained")
    ax = axes[0]
    title(ax, "A  |  A second gene labels the branch")
    ax.plot(x, hump, color=BLUE, label=r"$\mu_1(x)$: symmetric peak")
    ax.plot(x, x, color=ORANGE, label=r"$\mu_2(x)=x/L$")
    ax.axvline(.3, ls=":", color=MUTED, lw=1)
    ax.axvline(.7, ls=":", color=MUTED, lw=1)
    ax.set(xlabel="$x/L$", ylabel="Mean expression", xlim=(0, 1), ylim=(0, 1.05))
    ax.legend(loc="upper left", fontsize=12)
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
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 5.0), layout="constrained")
    ax = axes[0]
    title(ax, "A  |  Same expression range")
    for mean, color, label in [(x, BLUE, "Linear"), (steep, ORANGE, "Steep sigmoid")]:
        info = scalar_information(mean, sigma, x)
        ax.plot(x, mean, color=color, label=f"{label}: {info:.2f} bits")
    ax.set(xlabel="$x/L$", ylabel=r"$\mu(x)$", xlim=(0, 1), ylim=(0, 1.05))
    ax.legend(loc="upper left", fontsize=12)
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
    ax.legend(fontsize=12)
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
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 5.0), layout="constrained")
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
    axes[2].legend(fontsize=12)
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
    fig, axes = plt.subplots(1, 2, figsize=(11.7, 4.9), layout="constrained")
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
    axes[1].tick_params(axis="x", labelsize=13)
    for ax in axes:
        ax.grid(axis="y", alpha=.12)
    save(fig, "06_average_and_defects")


if __name__ == "__main__":
    for draw in [french_flag, gradient, molecular_readout, reading_noise, correlated_noise, ambiguity, information_budget, optimal_gradient,
                 diffusion_tradeoff, mean_and_worst]:
        draw()
        print(draw.__name__)
