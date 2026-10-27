"""Reproduce the stochastic-dynamics figures and independent numerical checks.

Run from any directory with NumPy, SciPy, and Matplotlib available. All numerical
examples are models rather than fitted experimental data. Random seeds are fixed.
"""

import json
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/qbio-stochastic-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/qbio-stochastic-cache")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
import numpy as np
from scipy.integrate import quad, solve_bvp
from scipy.stats import binom, poisson

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
from plot_style import configure_style, save_figure

FIGURES = ROOT / "figures"
BLUE, ORANGE, GREEN, PURPLE = "#3978a8", "#dc8540", "#4b9478", "#8b6bb0"


def save(fig, name):
    save_figure(fig, FIGURES, name)
    plt.close(fig)


def gillespie(k, gamma, stop, seed):
    rng = np.random.default_rng(seed)
    times, counts = [0.0], [0]
    while True:
        n = counts[-1]
        rate = k + gamma * n
        t = times[-1] + rng.exponential(1 / rate)
        if t > stop:
            times.append(stop)
            counts.append(n)
            break
        times.append(t)
        counts.append(n + (1 if rng.random() < k / rate else -1))
    return np.asarray(times), np.asarray(counts)


def walk_and_diffusion():
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0))
    m = 80
    x = np.arange(-m, m + 1, 2)
    prob = binom.pmf((x + m) // 2, m, 0.5)
    axes[0].bar(x, prob / 2, width=1.7, color=BLUE, alpha=0.75, label="Exact lattice density")
    fine = np.linspace(-35, 35, 401)
    axes[0].plot(fine, np.exp(-fine**2 / (2 * m)) / np.sqrt(2 * np.pi * m),
                 color=ORANGE, label="Diffusion limit")
    axes[0].set(xlim=(-32, 32), ylim=(0, 0.061), xlabel="Position $x/a$", ylabel="Probability density",
                title="A  Many paths, one distribution")
    axes[0].legend(loc="upper right", fontsize=12)
    steps = np.arange(0, 101, 10)
    axes[1].plot(steps, steps, color=ORANGE, label="$2Dt$, $D=1/2$")
    axes[1].scatter(steps, steps, color=BLUE, s=30, zorder=3, label="Exact discrete variance")
    axes[1].set(xlabel=r"Time $t/\tau$", ylabel="Mean squared displacement / $a^2$",
                title="B  Matching the spreading rate")
    axes[1].legend(loc="upper left", fontsize=12)
    save(fig, "01_walk_and_diffusion")


def hidden_state_memory():
    t = np.linspace(0, 12, 400)
    w, r1, r2 = 0.5, 0.2, 2.0
    s1, s2 = w * np.exp(-r1 * t), (1 - w) * np.exp(-r2 * t)
    survival = s1 + s2
    hazard = (r1 * s1 + r2 * s2) / survival
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0))
    axes[0].semilogy(t, survival, color=BLUE, label="Unresolved closed state")
    axes[0].semilogy(t, s1, "--", color=GREEN, label="Slow contribution")
    axes[0].semilogy(t, s2, ":", color=ORANGE, label="Fast contribution")
    axes[0].set(ylim=(0.02, 1.1), xlabel="Time since entry", ylabel="Survival probability",
                title="A  A mixture of lifetimes")
    axes[0].legend(fontsize=12)
    axes[1].plot(t, hazard, color=PURPLE)
    axes[1].axhline(r1, color=GREEN, linestyle="--", label="Slow exit rate $r_1$")
    axes[1].set(xlabel="Time already spent closed", ylabel="Opening hazard / time$^{-1}$",
                title="B  Observed memory without changing rates", ylim=(0, 1.2))
    axes[1].legend(fontsize=12)
    save(fig, "02_hidden_state_memory")


def birth_death_dynamics():
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2))
    t, n = gillespie(12, 1, 8, 204)
    axes[0].step(t, n, where="post", color=BLUE, linewidth=1.5, label="One reaction trajectory")
    fine = np.linspace(0, 8, 300)
    axes[0].plot(fine, 12 * (1 - np.exp(-fine)), color=ORANGE, label="Ensemble mean")
    axes[0].set(xlabel=r"Time $\gamma t$", ylabel="Molecule count $N$",
                title="A  Events and their mean", ylim=(0, max(n) + 3))
    axes[0].legend(fontsize=12)
    values = np.arange(30)
    for time, color in zip((0.2, 1, 4), (GREEN, ORANGE, BLUE)):
        axes[1].plot(values, poisson.pmf(values, 12 * (1 - np.exp(-time))),
                     "o-", markersize=3.5, color=color, label=f"$t={time:g}$")
    axes[1].set(xlabel="Molecule count $n$", ylabel="$P_n(t)$",
                title="B  The probability distribution evolves")
    axes[1].legend(fontsize=12)
    save(fig, "03_birth_death_dynamics")


def different_clocks():
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.0))
    mu = 12
    n = np.arange(30)
    axes[0].bar(n, poisson.pmf(n, mu), color=BLUE, alpha=0.75)
    axes[0].set(xlabel="Count $n$", ylabel="Stationary probability", title="A  Identical occupancy")
    t = np.linspace(0, 10, 300)
    omega = np.logspace(-2, 2, 400)
    for gamma, color in ((0.25, PURPLE), (2, GREEN)):
        label = rf"$\gamma={gamma:g}$"
        axes[1].plot(t, np.exp(-gamma * t), color=color, label=label)
        axes[2].loglog(omega, 2 * mu * gamma / (gamma**2 + omega**2), color=color, label=label)
    axes[1].set(xlabel="Lag time", ylabel=r"$C_N(t)/\mu$", title="B  Different memory")
    axes[2].set(xlabel=r"Angular frequency $\omega$", ylabel=r"$S_N(\omega)$",
                title="C  Different spectra")
    axes[1].legend(fontsize=12)
    axes[2].legend(fontsize=12)
    save(fig, "04_same_distribution_different_clocks")


def noise_filtering():
    rng = np.random.default_rng(306)
    dt, duration, tau_y = 0.005, 40, 1
    t = np.arange(0, duration + dt / 2, dt)
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.0))
    for tau, color, label in ((0.25, GREEN, "Fast input"), (4, PURPLE, "Slow input")):
        x, y = np.empty_like(t), np.empty_like(t)
        x[0], y[0] = rng.normal(), 0.0
        decay = np.exp(-dt / tau)
        for i in range(1, len(t)):
            x[i] = decay * x[i - 1] + np.sqrt(1 - decay**2) * rng.normal()
            y[i] = np.exp(-dt / tau_y) * y[i - 1] + (1 - np.exp(-dt / tau_y)) * x[i - 1]
        # Suppress initialization of the downstream filter in the plotted window.
        shown = t >= 8
        axes[0].plot(t[shown] - 8, x[shown], color=color, alpha=0.85, linewidth=1, label=label)
        axes[1].plot(t[shown] - 8, y[shown], color=color, linewidth=1.6, label=label)
    axes[0].set(xlabel="Time", ylabel="Input deviation $x$", title="A  Equal input variance")
    axes[1].set(xlabel="Time", ylabel="Output deviation $y$", title="B  One downstream timescale")
    axes[0].legend(fontsize=11, loc="upper right")
    ratios = np.logspace(-2, 2, 300)
    axes[2].semilogx(ratios, ratios / (1 + ratios), color=BLUE)
    axes[2].scatter([0.25, 4], [0.2, 0.8], color=[GREEN, PURPLE], s=55, zorder=4)
    axes[2].set(xlabel=r"$\tau_x/\tau_y$", ylabel=r"$\mathrm{Var}(y)/(g^2\sigma_x^2)$",
                title="C  Slow fluctuations accumulate", ylim=(0, 1.05))
    save(fig, "05_noise_filtering")


def potential(x):
    return (np.asarray(x)**2 - 1)**2 / 4


def first_passage_integral(D):
    def outer(y):
        inside = quad(lambda z: np.exp(-potential(z) / D), -2, y,
                      epsabs=1e-11, epsrel=1e-10)[0]
        return np.exp(potential(y) / D) * inside / D
    return quad(outer, -1, 0.8, epsabs=1e-8, epsrel=1e-9)[0]


def barrier_and_first_passage():
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2))
    q = np.linspace(-1.7, 1.6, 400)
    axes[0].plot(q, potential(q), color=BLUE)
    axes[0].scatter([-1, 0, 1], potential(np.array([-1, 0, 1])), color=[GREEN, ORANGE, GREEN])
    axes[0].axvline(0.8, color=PURPLE, linestyle="--", label="Absorbing target $c=0.8$")
    axes[0].annotate(r"Barrier $\Delta U=1/4$", xy=(0, 0.25), xytext=(-1.5, 0.55),
                     arrowprops={"arrowstyle": "->", "color": ORANGE}, color=ORANGE)
    axes[0].set(xlabel="State coordinate $q$", ylabel="Effective potential $U(q)$",
                title="A  Recovery and escape", ylim=(-0.035, 0.85))
    axes[0].legend(fontsize=11, loc="upper right")
    D = np.linspace(0.035, 0.17, 35)
    exact = np.array([first_passage_integral(d) for d in D])
    approximation = np.sqrt(2) * np.pi * np.exp(0.25 / D)
    axes[1].semilogy(D, exact, color=BLUE, label="First-passage integral")
    axes[1].semilogy(D, approximation, "--", color=ORANGE, label="Kramers approximation")
    axes[1].set(xlabel="Noise strength $D$", ylabel="Mean first-passage time",
                title="B  Exponentially long waiting")
    axes[1].legend(fontsize=12)
    save(fig, "06_barrier_and_first_passage")


def stationary_currents():
    fig, axes = plt.subplots(1, 2, figsize=(10.3, 4.8))
    angles = np.deg2rad([90, -30, 210])
    positions = np.column_stack((np.cos(angles), np.sin(angles)))
    for ax, a, b, title in zip(axes, (1, 1.8), (1, 0.2),
                                ("A  Detailed balance", "B  A directed stationary cycle")):
        for i, p in enumerate(positions):
            ax.add_patch(Circle(p, 0.35, facecolor="#eef4f7", edgecolor=BLUE, linewidth=1.7))
            ax.text(*p, f"${i}$\n$p^*=1/3$", ha="center", va="center", fontsize=13)
        for i in range(3):
            p, q = positions[i], positions[(i + 1) % 3]
            for start, end, rate, color in ((p, q, a, ORANGE), (q, p, b, GREEN)):
                direction = end - start
                norm = np.linalg.norm(direction)
                tangent = direction / norm
                normal = np.array([-tangent[1], tangent[0]])
                # Equal positive curvature separates the two opposite arrows.
                begin, finish = start + 0.38 * tangent, end - 0.38 * tangent
                ax.add_patch(FancyArrowPatch(begin, finish, arrowstyle="-|>",
                                             connectionstyle="arc3,rad=0.19",
                                             mutation_scale=15, color=color,
                                             linewidth=1.3 + 0.8 * rate))
                label = (start + end) / 2 - 0.23 * normal
                ax.text(*label, f"{rate:g}", color=color, ha="center", va="center", fontsize=14,
                        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
        ax.text(0, -1.14, f"Clockwise net current $J={(a-b)/3:.3g}$",
                ha="center", va="center", fontsize=13, color=PURPLE)
        ax.set(title=title, xlim=(-1.45, 1.45), ylim=(-1.35, 1.43), aspect="equal")
        ax.axis("off")
    save(fig, "07_same_occupancy_different_currents")


def numerical_checks():
    """Cross-check formulas using different representations of the process."""
    mu, gamma = 12, 0.7
    n = np.arange(120)
    p = poisson.pmf(n, mu)
    flux = mu * gamma * p[:-1] - gamma * n[1:] * p[1:]
    flux_error = float(np.max(np.abs(flux)))
    assert flux_error < 1e-12
    assert abs(p.sum() - 1) < 1e-12
    assert abs(np.sum((n - mu)**2 * p) - mu) < 1e-10

    # Integrate the output spectrum rather than reusing the moment expression.
    tau_x, tau_y, sigma2 = 4, 2, 12
    spectral_variance = quad(
        lambda omega: 2 * sigma2 * tau_x / ((1 + (omega * tau_x)**2)
                                           * (1 + (omega * tau_y)**2)) / np.pi,
        0, np.inf, epsabs=1e-10)[0]
    predicted_variance = sigma2 * tau_x / (tau_x + tau_y)
    assert abs(spectral_variance - predicted_variance) < 1e-9

    # Independent backward boundary-value solution versus the double integral.
    D = 0.08
    mesh = np.linspace(-2, 0.8, 700)
    def ode(q, y):
        return np.vstack((y[1], (-1 + (q**3 - q) * y[1]) / D))
    solution = solve_bvp(ode, lambda ya, yb: np.array([ya[1], yb[0]]), mesh,
                         np.zeros((2, len(mesh))), tol=1e-8, max_nodes=12000)
    assert solution.success, solution.message
    integral = first_passage_integral(D)
    backward_value = float(solution.sol(-1)[0])
    passage_error = abs(backward_value - integral) / integral
    assert passage_error < 1e-7

    # Solve the finite-state backward linear system independently of recursion.
    K, beta, delta = 24, 1.5, 1.0
    states = np.arange(1, K + 1)
    births = beta * states * (1 - states / K)
    deaths = delta * states
    matrix = np.diag(-(births + deaths))
    matrix += np.diag(births[:-1], 1) + np.diag(deaths[1:], -1)
    times_linear = np.linalg.solve(matrix, -np.ones(K))
    differences = np.zeros(K)
    differences[-1] = 1 / deaths[-1]
    for i in range(K - 2, -1, -1):
        differences[i] = (1 + births[i] * differences[i + 1]) / deaths[i]
    times_recursive = np.cumsum(differences)
    extinction_error = float(np.max(np.abs(times_linear - times_recursive)))
    assert extinction_error < 1e-9

    # Same occupancy and exit clock; distinct current and path irreversibility.
    rings = []
    for a, b in ((1, 1), (1.8, 0.2)):
        Q = np.array([[-a-b, a, b], [b, -a-b, a], [a, b, -a-b]])
        stationary_residual = float(np.max(np.abs(np.ones(3) @ Q / 3)))
        assert stationary_residual < 1e-14
        rings.append({"clockwise_rate": a, "counterclockwise_rate": b,
                      "probability_per_state": 1/3, "mean_residence_time": 1/(a+b),
                      "edge_current": (a-b)/3,
                      "entropy_rate_in_kB": (a-b) * np.log(a/b),
                      "stationary_residual": stationary_residual})

    report = {
        "seed_birth_death": 204, "seed_noise_filter": 306,
        "birth_death": {"mean": mu, "gamma": gamma, "max_stationary_edge_flux": flux_error},
        "filter": {"input_variance": sigma2, "tau_x": tau_x, "tau_y": tau_y,
                   "variance_from_spectrum": spectral_variance,
                   "variance_from_moments": predicted_variance},
        "first_passage": {"D": D, "reflecting_boundary": -2, "target": 0.8,
                          "start": -1, "integral": integral, "backward_BVP": backward_value,
                          "relative_difference": passage_error},
        "extinction": {"K": K, "beta": beta, "delta": delta,
                       "max_recursion_vs_linear_system_difference": extinction_error,
                       "time_from_n8": float(times_linear[7])},
        "rings": rings,
    }
    (ROOT / "numerical_checks.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    FIGURES.mkdir(exist_ok=True)
    configure_style()
    checks = numerical_checks()
    for build in (walk_and_diffusion, hidden_state_memory, birth_death_dynamics,
                  different_clocks, noise_filtering, barrier_and_first_passage, stationary_currents):
        build()
    print(json.dumps({"figures": 7, "checks": checks}, indent=2))


if __name__ == "__main__":
    main()
