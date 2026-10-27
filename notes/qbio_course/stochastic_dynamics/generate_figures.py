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
from scipy.linalg import expm, solve_continuous_lyapunov
from scipy.stats import binom, norm, poisson

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


def marginals_and_memory():
    rng = np.random.default_rng(810)
    count, r = 180, 0.94
    independent = rng.normal(size=count)
    correlated = np.empty(count)
    correlated[0] = rng.normal()
    for i in range(1, count):
        correlated[i] = r * correlated[i - 1] + np.sqrt(1 - r*r) * rng.normal()
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.1))
    axes[0].plot(independent, color=GREEN, linewidth=1.1, label="Independent")
    axes[0].plot(correlated, color=PURPLE, linewidth=1.7, label="Correlated")
    axes[0].set(xlabel="Time step", ylabel="State value", title="A  Different paths")
    axes[0].legend(fontsize=11)
    x = np.linspace(-3.5, 3.5, 300)
    axes[1].plot(x, norm.pdf(x), color=BLUE)
    axes[1].set(xlabel="State value", ylabel="One-time density",
                title="B  Both have the same Gaussian density")
    lag = np.arange(41)
    axes[2].plot(lag, r**lag, color=PURPLE, label=r"$r^{|\ell|}$, $r=0.94$")
    axes[2].plot(lag, (lag == 0).astype(float), "o-", color=GREEN,
                 markersize=3, label="Independent")
    axes[2].set(xlabel="Lag in time steps", ylabel="Normalized correlation",
                title="C  A distribution does not fix memory")
    axes[2].legend(fontsize=11)
    save(fig, "08_marginals_and_memory")


def persistence_and_diffusion():
    v, rate = 1.0, 0.5
    D = v*v / (2*rate)
    t = np.logspace(-2, 2, 350)
    msd = v*v/rate * (t + np.expm1(-2*rate*t)/(2*rate))
    fig, axes = plt.subplots(1, 2, figsize=(10.7, 4.2))
    axes[0].loglog(t, msd, color=BLUE, label="Persistent motion")
    short, long = t < 0.6, t > 2
    axes[0].loglog(t[short], v*v*t[short]**2, "--", color=ORANGE,
                   label=r"Ballistic: $v^2t^2$")
    axes[0].loglog(t[long], 2*D*t[long], ":", color=GREEN,
                   label=r"Diffusive: $2D_{\rm eff}t$")
    axes[0].set(xlabel="Time", ylabel="Mean squared displacement",
                title="A  Losing direction changes the scaling")
    axes[0].legend(fontsize=11)
    s = np.linspace(0, 7, 300)
    covariance = v*v*np.exp(-2*rate*s)
    axes[1].plot(s, covariance, color=PURPLE)
    axes[1].fill_between(s, 0, covariance, color=PURPLE, alpha=0.2)
    axes[1].text(2.5, 0.5, r"Area $=D_{\rm eff}=1$", fontsize=14, color=PURPLE)
    axes[1].set(xlabel="Time lag", ylabel="Velocity correlation",
                title="B  Memory sets the diffusion coefficient", ylim=(0, 1.08))
    save(fig, "09_persistence_and_diffusion")


def discrete_and_gaussian():
    fig, axes = plt.subplots(1, 2, figsize=(10.9, 4.1))
    for ax, mu, letter in zip(axes, (2, 30), ("A", "B")):
        n = np.arange(-5 if mu == 2 else -3, int(mu + 5*np.sqrt(mu)) + 1)
        bins = norm.cdf((n + .5 - mu)/np.sqrt(mu)) - norm.cdf((n - .5 - mu)/np.sqrt(mu))
        valid = n >= 0
        ax.bar(n[valid], poisson.pmf(n[valid], mu), width=.72, color=BLUE,
               alpha=.7, label="Exact Poisson")
        ax.plot(n, bins, "o-", color=ORANGE, markersize=3, label="Binned Gaussian")
        ax.axvspan(n[0]-.5, -.5, color="#e9dce8", alpha=.8)
        ax.set(xlabel="Count bin $n$", ylabel="Probability per unit bin",
               title=f"{letter}  Mean count $\\mu={mu}$", xlim=(n[0]-.5, n[-1]+.5))
        ax.legend(fontsize=11)
    save(fig, "10_discrete_and_gaussian")


def population_model(K=60, beta=1.5, delta=1):
    n = np.arange(K+1)
    births, deaths = beta*n*(1-n/K), delta*n
    Q = np.diag(-births-deaths) + np.diag(births[:-1], 1) + np.diag(deaths[1:], -1)
    eigenvalues, eigenvectors = np.linalg.eig(Q[1:, 1:].T)
    index = np.argmax(eigenvalues.real)
    assert abs(eigenvalues[index].imag) < 1e-10
    kappa = -float(eigenvalues[index].real)
    q = eigenvectors[:, index].real
    if q.sum() < 0:
        q = -q
    assert q.min() > -1e-10
    q = np.maximum(q, 0)
    q /= q.sum()
    return Q, q, kappa


def quasistationary_extinction():
    Q, q, kappa = population_model()
    initial = np.zeros(61)
    initial[20] = 1
    scaled = np.linspace(0, 4, 100)
    survival = np.array([(expm(Q.T*t/kappa) @ initial)[1:].sum() for t in scaled])
    fig, axes = plt.subplots(1, 2, figsize=(10.9, 4.2))
    axes[0].semilogy(scaled, survival, color=BLUE, label="Finite-state evolution")
    axes[0].semilogy(scaled, np.exp(-scaled), "--", color=ORANGE,
                     label="Quasistationary exponential")
    axes[0].set(xlabel=r"Time / quasistationary lifetime $\kappa^{-1}$",
                ylabel="Survival probability", title="A  Probability continues to leave")
    axes[0].legend(fontsize=11)
    states = np.arange(1, 61)
    for time, color in zip((.1, 1, 3), (GREEN, ORANGE, PURPLE)):
        p = (expm(Q.T*time/kappa) @ initial)[1:]
        axes[1].plot(states, p/p.sum(), color=color, linewidth=1.6,
                     label=rf"$\kappa t={time:g}$")
    axes[1].plot(states, q, "--", color=BLUE, linewidth=2.2, label="QSD eigenvector")
    axes[1].set(xlabel="Population $n$", ylabel="Conditional probability",
                title="B  Surviving samples keep a stable shape", xlim=(0, 45))
    axes[1].legend(fontsize=10)
    save(fig, "11_quasistationary_extinction")


def gaussian_density_and_circulation():
    fig, axes = plt.subplots(1, 3, figsize=(12.9, 4.3))
    grid = np.linspace(-2.8, 2.8, 151)
    X, Y = np.meshgrid(grid, grid)
    p = np.exp(-(X*X+Y*Y)/2)/(2*np.pi)
    for ax, omega, letter in zip(axes[:2], (0, 2), ("A", "B")):
        ax.contourf(X, Y, p, levels=9, cmap="Blues", alpha=.8)
        ax.contour(X, Y, p, levels=[.02, .05, .1, .14], colors="white", linewidths=.7)
        if omega:
            coarse = np.linspace(-2.3, 2.3, 9)
            u, v = np.meshgrid(coarse, coarse)
            density = np.exp(-(u*u+v*v)/2)/(2*np.pi)
            ax.quiver(u, v, -omega*v*density, omega*u*density,
                      color=ORANGE, angles="xy", scale_units="xy", scale=.5, width=.009)
        else:
            ax.text(0, -.1, r"$\mathbf{J}_*=0$", color=ORANGE, fontsize=20,
                    ha="center", va="center")
        ax.set(xlabel="$x$", ylabel="$y$", aspect="equal",
               title=f"{letter}  Same density, $\\omega_0={omega}$")
    time = np.linspace(0, 6, 300)
    cross = np.exp(-time)*np.sin(2*time)
    axes[2].plot(time, cross, color=PURPLE, label="$C_{yx}(t)$")
    axes[2].plot(time, -cross, color=GREEN, label="$C_{xy}(t)$")
    axes[2].axhline(0, color="#aab4bf", linewidth=.8)
    axes[2].set(xlabel="Time lag", ylabel="Cross-correlation",
                title="C  Time order reveals circulation")
    axes[2].legend(fontsize=11)
    save(fig, "12_gaussian_density_and_circulation")


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

    # Persistence: integrate velocity covariance independently of the MSD formula.
    time, speed, rate = 1.7, 1.0, .5
    msd_integral = 2*quad(lambda s: (time-s)*speed**2*np.exp(-2*rate*s), 0, time)[0]
    msd_formula = speed**2/rate*(time-(1-np.exp(-2*rate*time))/(2*rate))
    assert abs(msd_integral-msd_formula) < 1e-12

    # Bursts: solve a truncated exact stationary generator, not moment equations.
    nmax, batch, arrival = 150, 3, 4.0
    burst_Q = np.zeros((nmax+1, nmax+1))
    for n_i in range(nmax+1):
        if n_i + batch <= nmax:
            burst_Q[n_i, n_i+batch] = arrival
        if n_i:
            burst_Q[n_i, n_i-1] = n_i
        burst_Q[n_i, n_i] = -burst_Q[n_i].sum()
    stationarity = burst_Q.T.copy()
    stationarity[-1] = 1
    rhs = np.zeros(nmax+1)
    rhs[-1] = 1
    burst_p = np.linalg.solve(stationarity, rhs)
    burst_mean = float(np.arange(nmax+1) @ burst_p)
    burst_variance = float((np.arange(nmax+1)-burst_mean)**2 @ burst_p)
    assert abs(burst_mean-12) < 1e-9 and abs(burst_variance-24) < 1e-8

    # Quasistationary lifetime: eigenvalue versus weighted backward times.
    population_Q, quasi, kappa = population_model()
    backward_times = np.linalg.solve(population_Q[1:, 1:], -np.ones(60))
    quasi_mean_time = float(quasi @ backward_times)
    assert abs(kappa*quasi_mean_time-1) < 1e-9

    # Stable shear and rotating OU: solve matrix Lyapunov equations independently.
    lam, D_ou, shear, omega = 1.0, 1.0, 3.0, 2.0
    A_shear = np.array([[-lam, shear], [0, -lam]])
    cov_shear = solve_continuous_lyapunov(A_shear, -2*D_ou*np.eye(2))
    predicted_shear = np.array([[D_ou/lam+shear**2*D_ou/(2*lam**3), shear*D_ou/(2*lam**2)],
                               [shear*D_ou/(2*lam**2), D_ou/lam]])
    assert np.max(np.abs(cov_shear-predicted_shear)) < 1e-12
    rotation = np.array([[0., -1.], [1., 0.]])
    A_rotation = -lam*np.eye(2)+omega*rotation
    cov_rotation = solve_continuous_lyapunov(A_rotation, -2*D_ou*np.eye(2))
    assert np.max(np.abs(cov_rotation-D_ou/lam*np.eye(2))) < 1e-12
    lag = .7
    C = expm(A_rotation*lag) @ cov_rotation
    assert abs(C[1, 0]-D_ou/lam*np.exp(-lam*lag)*np.sin(omega*lag)) < 1e-12

    # Splitting probability: an absorbing/absorbing BVP versus its quadrature.
    split_D = .1
    split_mesh = np.linspace(-1, 1, 400)
    def split_ode(q, y):
        return np.vstack((y[1], (q**3-q)*y[1]/split_D))
    splitting = solve_bvp(split_ode, lambda ya, yb: np.array([ya[0], yb[0]-1]),
                          split_mesh, np.vstack(((split_mesh+1)/2, np.full_like(split_mesh,.5))),
                          tol=1e-8)
    assert splitting.success
    split_integral = quad(lambda y: np.exp(potential(y)/split_D), -1, .2)[0]
    split_integral /= quad(lambda y: np.exp(potential(y)/split_D), -1, 1)[0]
    split_backward = float(splitting.sol(.2)[0])
    assert abs(split_integral-split_backward) < 1e-8
    R0 = 1.5
    action_integral = quad(lambda x: np.log(R0*(1-x)), 0, 1-1/R0)[0]
    action_formula = np.log(R0)-1+1/R0
    assert abs(action_integral-action_formula) < 1e-12

    # Promoter switching: stationary joint master equation versus moment formulas.
    promoter_checks = []
    cutoff, synthesis, removal = 100, 20.0, 1.0
    for on_rate, off_rate in ((.2, .2), (20.0, 20.0)):
        joint_Q = np.zeros((2*(cutoff+1), 2*(cutoff+1)))
        for number in range(cutoff+1):
            for state in (0, 1):
                index = 2*number+state
                joint_Q[index, 2*number+1-state] = on_rate if state == 0 else off_rate
                if number:
                    joint_Q[index, index-2] = removal*number
                if state and number < cutoff:
                    joint_Q[index, index+2] = synthesis
                joint_Q[index, index] = -joint_Q[index].sum()
        constraints = joint_Q.T.copy()
        constraints[-1] = 1
        normalizer = np.zeros(len(joint_Q))
        normalizer[-1] = 1
        joint_p = np.linalg.solve(constraints, normalizer)
        counts = np.repeat(np.arange(cutoff+1), 2)
        promoter_states = np.tile([0, 1], cutoff+1)
        mean = float(counts @ joint_p)
        variance = float((counts-mean)**2 @ joint_p)
        pi = on_rate/(on_rate+off_rate)
        covariance = float((counts-mean)*(promoter_states-pi) @ joint_p)
        predicted_covariance = synthesis*pi*(1-pi)/(removal+on_rate+off_rate)
        predicted_fano = 1+synthesis*(1-pi)/(removal+on_rate+off_rate)
        assert abs(mean-synthesis*pi/removal) < 1e-9
        assert abs(covariance-predicted_covariance) < 1e-9
        assert abs(variance/mean-predicted_fano) < 1e-8
        promoter_checks.append({"on_rate": on_rate, "off_rate": off_rate,
                                "mean_from_master_equation": mean,
                                "variance_from_master_equation": variance,
                                "fano_from_master_equation": variance/mean,
                                "fano_from_moments": predicted_fano})

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
        "persistence": {"MSD_from_covariance_integral": msd_integral, "MSD_formula": msd_formula},
        "burst_process": {"batch_size": batch, "arrival_rate": arrival,
                          "mean_from_stationary_generator": burst_mean,
                          "variance_from_stationary_generator": burst_variance},
        "quasistationarity": {"K": 60, "kappa": kappa, "lifetime": 1/kappa,
                              "lifetime_from_backward_times": quasi_mean_time},
        "linear_OU": {"shear_covariance": cov_shear.tolist(),
                      "rotating_covariance": cov_rotation.tolist(),
                      "rotating_cross_covariance_at_t0p7": C.tolist()},
        "splitting_probability": {"D": split_D, "start": .2,
                                  "quadrature": split_integral, "backward_BVP": split_backward},
        "extinction_action": {"R0": R0, "quadrature": action_integral,
                              "closed_form": action_formula},
        "promoter_switching": promoter_checks,
    }
    (ROOT / "numerical_checks.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    FIGURES.mkdir(exist_ok=True)
    configure_style()
    checks = numerical_checks()
    for build in (walk_and_diffusion, hidden_state_memory, birth_death_dynamics,
                  different_clocks, noise_filtering, barrier_and_first_passage, stationary_currents,
                  marginals_and_memory, persistence_and_diffusion, discrete_and_gaussian,
                  quasistationary_extinction, gaussian_density_and_circulation):
        build()
    print(json.dumps({"figures": 12, "checks": checks}, indent=2))


if __name__ == "__main__":
    main()
