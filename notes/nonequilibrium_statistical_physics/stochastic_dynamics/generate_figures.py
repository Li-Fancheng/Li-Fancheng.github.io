"""Reproduce the stochastic-dynamics figures and quantitative checks.

Run from any directory with the packages in requirements.txt. All model
parameters are dimensionless. Figures are original calculations, not copied
from the reference literature. SVG text is outlined for portable rendering.
"""

import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "neq-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "neq-font-cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid, quad, simpson
from scipy.linalg import eigh_tridiagonal


ROOT = Path(__file__).resolve().parent
FIGURES = ROOT / "figures"
BLUE, ORANGE, GREEN, INK = "#2563a6", "#d66b32", "#198477", "#273b4b"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 12,
    "axes.titlesize": 13, "axes.labelsize": 12,
    "axes.spines.top": False, "axes.spines.right": False,
    "text.color": INK, "axes.labelcolor": INK, "axes.titlecolor": INK,
    "legend.frameon": False, "legend.fontsize": 10,
    "lines.linewidth": 2, "svg.fonttype": "path",
    "svg.hashsalt": "nonequilibrium-stochastic-dynamics",
    "savefig.facecolor": "white",
})


def potential(x):
    return (np.asarray(x) ** 2 - 1) ** 2 / 4


def gradient(x):
    return np.asarray(x) ** 3 - np.asarray(x)


def save(fig, name):
    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / f"{name}.svg", bbox_inches="tight", metadata={"Date": None})
    fig.savefig(FIGURES / f"{name}.png", bbox_inches="tight", dpi=155)
    plt.close(fig)


def gaussian(x, mean, variance):
    return np.exp(-(x - mean) ** 2 / (2 * variance)) / np.sqrt(2 * np.pi * variance)


def ou_figures():
    gamma, diffusion, h = 1.0, 0.4, 0.008
    t = np.arange(0, 6 + h / 2, h)
    mean = 2 * np.exp(-gamma * t)
    variance = diffusion / gamma * (1 - np.exp(-2 * gamma * t))
    rng = np.random.default_rng(20261028)
    paths = np.empty((3, len(t)))
    paths[:, 0] = 2
    decay = np.exp(-gamma * h)
    amplitude = np.sqrt(diffusion / gamma * (1 - decay ** 2))
    for n in range(len(t) - 1):
        paths[:, n + 1] = decay * paths[:, n] + amplitude * rng.standard_normal(3)
    fig, ax = plt.subplots(2, 2, figsize=(10, 7.3), layout="constrained")
    for path in paths:
        ax[0, 0].plot(t, path, alpha=0.6, lw=0.9)
    ax[0, 0].fill_between(t, mean - np.sqrt(variance), mean + np.sqrt(variance), color=BLUE, alpha=0.12)
    ax[0, 0].plot(t, mean, color=INK, label="ensemble mean")
    ax[0, 0].set(title="A  Sample paths and ensemble statistics", xlabel="time", ylabel="$X_t$")
    ax[0, 0].legend()
    ax[0, 1].plot(t, variance, color=BLUE, label="$v(t)$")
    ax[0, 1].axhline(diffusion / gamma, color=ORANGE, ls="--", label="$D/\\gamma$")
    ax[0, 1].set(title="B  Continuous injection, finite variance", xlabel="time", ylabel="variance", ylim=(0, 0.47))
    ax[0, 1].legend(loc="lower right")
    lag = np.linspace(-4, 4, 600)
    ax[1, 0].plot(lag, diffusion / gamma * np.exp(-gamma * abs(lag)), color=BLUE)
    ax[1, 0].set(title="C  Stationary correlation", xlabel="lag", ylabel="$C(\\tau)$")
    ax[1, 1].plot(lag[lag < 0], np.zeros(sum(lag < 0)), color=ORANGE)
    ax[1, 1].plot(lag[lag > 0], np.exp(-gamma * lag[lag > 0]), color=ORANGE)
    ax[1, 1].plot([0, 0], [0, 1], color=ORANGE, ls=":", lw=1)
    ax[1, 1].set(title="D  Causal response to drift forcing", xlabel="lag", ylabel="$R(\\tau)$")
    save(fig, "01_ou_fluctuations")

    x = np.linspace(-2.4, 2.4, 800)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.9), layout="constrained")
    for time, color in zip([0, 0.12, 0.5, 2], [INK, ORANGE, GREEN, BLUE]):
        v = 0.04 * np.exp(-2 * time) + diffusion * (1 - np.exp(-2 * time))
        ax[0].plot(x, gaussian(x, 0, v), color=color, label=f"$t={time:g}$")
    ax[0].plot(x, gaussian(x, 0, diffusion), color=INK, ls="--", label="stationary")
    ax[0].set(title="A  Probability spreads out from a narrow peak", xlabel="$x$", ylabel="$p(x,t)$")
    ax[0].legend()
    time = 0.12
    v = 0.04 * np.exp(-2 * time) + diffusion * (1 - np.exp(-2 * time))
    density = gaussian(x, 0, v)
    drift = -gamma * x * density
    spreading = diffusion * x / v * density
    ax[1].plot(x, drift, color=ORANGE, label="drift current")
    ax[1].plot(x, spreading, color=BLUE, label="diffusive current")
    ax[1].plot(x, drift + spreading, color=GREEN, label="net current")
    ax[1].axhline(0, color=INK, alpha=0.3, lw=1)
    ax[1].set(title="B  Competing currents at $t=0.12$", xlabel="$x$", ylabel="$J(x,t)$")
    ax[1].legend()
    save(fig, "02_ou_transport")


def reversible_grid(diffusion, size):
    """Reflecting birth-death discretization and its symmetric similarity.

    Neighbour rates D/dx² exp[-(U_j-U_i)/(2D)] obey exact discrete detailed
    balance. The continuum limit has drift -U' and diffusion D. Omitted jumps
    beyond the endpoint implement reflecting boundaries.
    """
    edges = np.linspace(-2, 2, size + 1)
    dx = edges[1] - edges[0]
    x = (edges[1:] + edges[:-1]) / 2
    u = potential(x)
    right = diffusion / dx ** 2 * np.exp(-np.diff(u) / (2 * diffusion))
    left = diffusion / dx ** 2 * np.exp(np.diff(u) / (2 * diffusion))
    diagonal = -np.r_[right, 0] - np.r_[0, left]
    off = np.sqrt(right * left)
    eigenvalues, vectors = eigh_tridiagonal(diagonal, off)
    pi = np.exp(-u / diffusion)
    pi /= pi.sum()
    return x, dx, pi, eigenvalues, vectors


def bistable_figure():
    diffusion = 0.07
    x, dx, pi, ev, vec = reversible_grid(diffusion, 401)
    initial = gaussian(x, -0.55, 0.06 ** 2) * dx
    initial /= initial.sum()
    root_pi = np.sqrt(pi)
    coefficients = vec.T @ (initial / root_pi)

    def evolve(time):
        return root_pi * (vec @ (np.exp(np.minimum(ev, 0) * time) * coefficients))

    fig, ax = plt.subplots(3, 1, figsize=(9, 9), layout="constrained")
    ax[0].plot(x, potential(x), color=INK)
    ax[0].set(title="A  Bistable potential: fast local confinement, rare crossing", ylabel="$U(x)$", ylim=(-0.02, 0.65), xlim=(-1.7, 1.7))
    for time, color in zip([0.3, 2, 40, 200, 600], [ORANGE, GREEN, BLUE, "#8d62ad", INK]):
        ax[1].plot(x, evolve(time) / dx, color=color, label=f"$t={time:g}$")
    ax[1].plot(x, pi / dx, color=INK, ls="--", label="stationary")
    ax[1].set(title="B  Local shape forms before the two wells equilibrate", ylabel="$p(x,t)$", xlim=(-1.7, 1.7))
    ax[1].legend(ncol=3)
    times = np.geomspace(0.05, 800, 180)
    left_mass = np.array([evolve(time)[x < 0].sum() for time in times])
    k = np.sqrt(2) / (2 * np.pi) * np.exp(-0.25 / diffusion)
    ax[2].semilogx(times, left_mass, color=BLUE, label="left-well mass (grid)")
    ax[2].semilogx(times, 1 - left_mass, color=ORANGE, label="right-well mass (grid)")
    ax[2].semilogx(times, 0.5 * (1 + np.exp(-2 * k * times)), color=INK, ls="--", label="two-state / Kramers")
    ax[2].set(title="C  The well label is a much slower variable", xlabel="time (log scale)", ylabel="probability", ylim=(-0.03, 1.03))
    ax[2].legend()
    save(fig, "03_bistable_transport")
    masses = [evolve(time) for time in [0, 0.3, 40, 600]]
    return {
        "probability_mass_max_error": float(max(abs(p.sum() - 1) for p in masses)),
        "minimum_probability": float(min(p.min() for p in masses)),
        "stationary_generator_residual": float(np.max(abs(
            (-np.r_[diffusion / dx ** 2 * np.exp(-np.diff(potential(x)) / (2 * diffusion)), 0]
             - np.r_[0, diffusion / dx ** 2 * np.exp(np.diff(potential(x)) / (2 * diffusion))]) * root_pi
            + np.r_[diffusion / dx ** 2 * root_pi[1:], 0]
            + np.r_[0, diffusion / dx ** 2 * root_pi[:-1]]))),
    }


def path_figure():
    # Reverse deterministic flow from -0.95 to -0.05; neither endpoint is fixed.
    start, end = -0.95, -0.05
    initial_ratio = 1 / start ** 2 - 1
    duration = 0.5 * np.log((1 / end ** 2 - 1) / initial_ratio)
    t = np.linspace(0, duration, 5001)
    xstar = -1 / np.sqrt(1 + initial_ratio * np.exp(2 * t))
    velocity = gradient(xstar)
    wiggle = np.sin(np.pi * t / duration)
    wiggle_derivative = np.pi / duration * np.cos(np.pi * t / duration)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.1), layout="constrained")
    x = np.linspace(-1.3, 1.3, 700)
    ax[0].plot(x, potential(x), color=INK)
    ax[0].scatter([-1, 0, 1], potential(np.array([-1, 0, 1])), color=[BLUE, ORANGE, BLUE], zorder=3)
    for left, right, color, offset in [(-0.88, -0.34, ORANGE, 0.03), (-0.34, -0.88, BLUE, -0.03), (0.32, 0.85, BLUE, 0.03)]:
        ax[0].annotate("", (right, potential(right) + offset), (left, potential(left) + offset), arrowprops={"arrowstyle": "->", "color": color, "lw": 2})
    ax[0].text(-0.84, 0.35, "reverse drift: uphill", color=ORANGE, fontsize=10)
    ax[0].text(-1.19, 0.11, "drift: downhill", color=BLUE, fontsize=10)
    ax[0].set(title="A  Minimum-cost escape climbs against drift", xlabel="$x$", ylabel="$U(x)$", ylim=(-0.06, 0.42))
    actions = []
    for amplitude, color, label in [(0, INK, "reverse drift"), (0.13, BLUE, "positive perturbation"), (-0.13, ORANGE, "negative perturbation")]:
        candidate = xstar + amplitude * wiggle
        derivative = velocity + amplitude * wiggle_derivative
        action = simpson((derivative + gradient(candidate)) ** 2 / 4, x=t)
        actions.append(float(action))
        ax[1].plot(t, candidate, color=color, label=f"{label}, $I={action:.3f}$")
    difference = float(potential(end) - potential(start))
    ax[1].set(title=f"B  Same endpoints and duration; $\\Delta U={difference:.3f}$", xlabel="time", ylabel="candidate path $x(t)$")
    ax[1].legend(loc="lower right", fontsize=9)
    save(fig, "04_optimal_paths")
    return {"endpoint_potential_difference": difference, "candidate_actions": actions,
            "reverse_path_action_error": abs(actions[0] - difference)}


def mfpt(diffusion, start=-1, target=0.8, left=-2):
    # A fine cumulative inner integral makes evaluating the exact formula cheap.
    z = np.linspace(left, target, 24001)
    inner = cumulative_trapezoid(np.exp(-potential(z) / diffusion), z, initial=0)
    value, _ = quad(lambda y: np.exp(potential(y) / diffusion) * np.interp(y, z, inner) / diffusion,
                    start, target, epsabs=1e-7, epsrel=2e-8, limit=150)
    return value


def spectrum_figure():
    diffusions = np.linspace(0.025, 0.12, 28)
    exact = np.array([mfpt(d) for d in diffusions])
    approximation = np.sqrt(2) * np.pi * np.exp(0.25 / diffusions)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), layout="constrained")
    order = np.argsort(1 / diffusions)
    ax[0].semilogy(1 / diffusions[order], exact[order], color=BLUE, label="first-passage integral")
    ax[0].semilogy(1 / diffusions[order], approximation[order], color=ORANGE, ls="--", label="Kramers asymptotics")
    ax[0].set(title="A  The barrier sets the exponential scale", xlabel="$1/D$", ylabel="mean first-passage time")
    ax[0].legend()
    spectra = {}
    for size in [201, 401, 801]:
        _, _, _, eigenvalues, _ = reversible_grid(0.05, size)
        spectra[size] = -eigenvalues[-2:-7:-1]
    rates = spectra[801]
    modes = np.arange(1, 6)
    ax[1].semilogy(modes, modes, "o-", color=BLUE, label="OU, $\\gamma=1$")
    ax[1].semilogy(modes, rates, "o-", color=ORANGE, label="double well, $D=0.05$")
    ax[1].set(title="B  A slow inter-well mode separates out", xlabel="nonzero mode number", ylabel="relaxation rate", xticks=modes)
    ax[1].legend(loc="lower right")
    save(fig, "05_first_passage_and_spectrum")
    top_time = mfpt(0.025, target=0)
    crossed_time = mfpt(0.025)
    return {
        "bistable_rates_D_005": {str(n): v.tolist() for n, v in spectra.items()},
        "slow_rate_relative_change_401_to_801": float(abs(spectra[401][0] / spectra[801][0] - 1)),
        "kramers_relative_error_D_0025": float(abs(approximation[0] / exact[0] - 1)),
        "target_top_vs_crossing_time_ratio_D_0025": float(crossed_time / top_time),
    }


def analytical_checks():
    gamma, diffusion, step = 1.3, 0.4, 0.1
    discrete_variance = 2 * diffusion * step / (1 - (1 - gamma * step) ** 2)
    predicted_variance = diffusion / (gamma * (1 - gamma * step / 2))
    x0, xf, duration = 0.3, 1.2, 1.7
    t = np.linspace(0, duration, 6001)
    delta = xf - x0 * np.exp(-gamma * duration)
    bridge = x0 * np.exp(-gamma * t) + delta * np.sinh(gamma * t) / np.sinh(gamma * duration)
    derivative = -gamma * x0 * np.exp(-gamma * t) + delta * gamma * np.cosh(gamma * t) / np.sinh(gamma * duration)
    action = simpson((derivative + gamma * bridge) ** 2 / 4, x=t)
    predicted_action = gamma * delta ** 2 / (2 * (1 - np.exp(-2 * gamma * duration)))
    return {"euler_ou_variance_identity_error": abs(discrete_variance - predicted_variance),
            "ou_bridge_action_identity_error": float(abs(action - predicted_action))}


if __name__ == "__main__":
    ou_figures()
    checks = {"bistable_transport": bistable_figure(), "path_action": path_figure(),
              "first_passage_and_spectrum": spectrum_figure(), "analytical": analytical_checks()}
    assert checks["bistable_transport"]["probability_mass_max_error"] < 1e-8
    assert checks["bistable_transport"]["minimum_probability"] > -1e-10
    assert checks["bistable_transport"]["stationary_generator_residual"] < 1e-8
    assert checks["path_action"]["reverse_path_action_error"] < 1e-8
    assert min(checks["path_action"]["candidate_actions"][1:]) > checks["path_action"]["candidate_actions"][0]
    assert checks["first_passage_and_spectrum"]["slow_rate_relative_change_401_to_801"] < 0.002
    assert checks["analytical"]["ou_bridge_action_identity_error"] < 1e-8
    assert checks["analytical"]["euler_ou_variance_identity_error"] < 1e-12
    (ROOT / "numerical_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(json.dumps(checks, indent=2))
