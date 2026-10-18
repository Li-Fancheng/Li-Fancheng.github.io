"""Figures 01--03: unbiased and biased random walks.

Run independently:
    python code/random_walk.py
"""
import numpy as np
import matplotlib.pyplot as plt
from common import plot_style, finish, normal_pdf


def make_random_walk_figures():
    plot_style()
    rng = np.random.default_rng(2026)
    a, tau, n_steps, n_walkers = 1.0, 1.0, 400, 12000
    steps = np.where(rng.random((n_walkers, n_steps)) < .5, -1, 1).astype(np.int8)
    x = np.column_stack([np.zeros(n_walkers, dtype=np.int16),
                         np.cumsum(steps, axis=1, dtype=np.int16)])
    t = np.arange(n_steps+1)*tau
    D = a*a/(2*tau)

    # FIG 01: single-particle trajectories. No smoothing or imposed drift.
    fig, ax = plt.subplots(figsize=(9.5, 4.0))
    for k in range(16):
        ax.plot(t, x[k], lw=1.05, alpha=.72)
    ax.axhline(0, color='black', lw=.8, alpha=.25)
    ax.set(xlabel='Time (steps)', ylabel='Position $x/a$',
           title='Unbiased random walk: individual trajectories', xlim=(0, n_steps))
    ax.text(.02, .96, '$p_+=p_-=1/2$', transform=ax.transAxes, va='top',
            bbox={'boxstyle':'round,pad=.4', 'fc':'white', 'ec':'0.84'})
    finish(fig, '01_random_walk_trajectories')

    # FIG 02: mean squared displacement, and spatial distribution.
    msd = np.mean(x.astype(float)**2, axis=0)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.9, 4.0),
                                    gridspec_kw={'width_ratios':[1, 1.04]})
    ax1.plot(t, msd, label='Simulation', lw=2)
    ax1.plot(t, 2*D*t, '--', color='black', label='$2Dt$', lw=1.6)
    ax1.set(xlabel='Time $t$', ylabel=r'$\langle x^2(t)\rangle$',
            title='Mean squared displacement')
    ax1.legend(frameon=False)
    end = x[:, -1].astype(float)
    bins = np.arange(-75, 77, 4)
    ax2.hist(end, bins=bins, density=True, alpha=.45,
             label='12,000 walkers', color='tab:blue')
    gx = np.linspace(-80, 80, 400)
    ax2.plot(gx, normal_pdf(gx, 0, 2*D*n_steps*tau), 'k--', lw=1.8,
             label='Diffusion Gaussian')
    ax2.set(xlabel='Position $x$', ylabel='Probability density',
            title=rf'Distribution at $t={n_steps}$')
    ax2.legend(frameon=False)
    fig.tight_layout(w_pad=2)
    finish(fig, '02_random_walk_to_diffusion')

    # FIG 03: keep step length and time unchanged, only change p_right.
    n_steps, n_walkers = 240, 14000
    p_right = .58
    rng = np.random.default_rng(2027)
    steps_u = np.where(rng.random((n_walkers, n_steps)) < .5, 1, -1)
    steps_b = np.where(rng.random((n_walkers, n_steps)) < p_right, 1, -1)
    xu = steps_u.sum(axis=1)
    xb = steps_b.sum(axis=1)
    v = (2*p_right-1)*a/tau
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.9, 4.0))
    edges = np.arange(-75, 100, 4)
    ax1.hist(xu, bins=edges, density=True, alpha=.48, label='$p_+=0.50$')
    ax1.hist(xb, bins=edges, density=True, alpha=.45, label='$p_+=0.58$')
    ax1.axvline(v*n_steps, color='black', ls='--', lw=1.2,
                label=r'$\langle x\rangle=vt$')
    ax1.set(xlabel='Position $x$', ylabel='Probability density',
            title='Same randomness, different bias')
    ax1.legend(frameon=False)
    tpoints = np.arange(n_steps+1)
    ax2.plot(tpoints, np.zeros_like(tpoints), lw=2, label='Unbiased mean')
    ax2.plot(tpoints, v*tpoints, lw=2, label='Biased mean: $vt$')
    ax2.set(xlabel='Time $t$', ylabel=r'Mean displacement $\langle x\rangle$',
            title='Drift emerges from unequal jump probabilities')
    ax2.legend(frameon=False)
    fig.tight_layout(w_pad=2)
    finish(fig, '03_biased_random_walk')
    print(f'[RW] diffusion coefficient D={D:.3f}; biased velocity v={v:.3f}')


if __name__ == '__main__':
    make_random_walk_figures()
