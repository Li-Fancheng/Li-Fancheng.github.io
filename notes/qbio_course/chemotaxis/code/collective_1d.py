"""Figures 06 and 08: minimal 1D Keller--Segel model on a periodic ring.

Cell density evolves by a finite-volume conservative flux. The chemotactic
flux uses an upwind discretization; diffusion uses centered differences.
Chemical field is updated by an explicit reaction-diffusion step.
An adaptive time step respects positivity / CFL bounds.

The snapshots are not presented as equilibria. Figure 06 shows the
*transient amplification* of a small density fluctuation.
"""
import numpy as np
import matplotlib.pyplot as plt
from common import (plot_style, finish, advdiff_flux_1d,
                    periodic_helmholtz_1d)

PARAMS = dict(L=20., nx=192, D=.3, Dc=1., alpha=1., kappa=1., rho0=1.)


def initial_state(par=PARAMS):
    L, nx = par['L'], par['nx']
    x = np.arange(nx)*(L/nx)
    # Explicit low-amplitude Fourier-mode perturbations: same in both runs.
    perturb = (.025*np.cos(2*np.pi*2*x/L+.1)
               + .025*np.cos(2*np.pi*3*x/L+1.7)
               + .015*np.cos(2*np.pi*4*x/L-1.2))
    rho = par['rho0']*(1+perturb)
    c = periodic_helmholtz_1d(rho, L, par['Dc'], par['alpha'], par['kappa'])
    return x, rho, c


def simulate_collective_1d(chi, total_time=55., nframes=221, par=PARAMS):
    x, rho, c = initial_state(par)
    dx = par['L']/par['nx']
    D, Dc, alpha, kappa = par['D'], par['Dc'], par['alpha'], par['kappa']
    times = np.linspace(0, total_time, nframes)
    density_snapshots = [rho.copy()]
    c_snapshots = [c.copy()]
    t = 0.
    nsteps = 0
    # Chemical positivity / diffusion bound for forward Euler.
    c_limit = .80/(2*Dc/dx**2+kappa)
    for target in times[1:]:
        while t < target-1e-11:
            flux, u = advdiff_flux_1d(rho, c, dx, D, chi)
            # Conservative safe upper bound on outgoing cell flux rate.
            rho_limit = .8/(2*D/dx**2 + 2*np.max(np.abs(u))/dx + 1e-14)
            dt = min(.0027, c_limit, rho_limit, target-t)
            rho_next = rho-dt*(flux-np.roll(flux, 1))/dx
            lapc = (np.roll(c, -1)-2*c+np.roll(c, 1))/dx**2
            c_next = c+dt*(Dc*lapc+alpha*rho-kappa*c)
            rho, c = rho_next, c_next
            t += dt
            nsteps += 1
        density_snapshots.append(rho.copy())
        c_snapshots.append(c.copy())
    dens = np.asarray(density_snapshots)
    c_hist = np.asarray(c_snapshots)
    assert dens.min() > -1e-10, f'Negative density: {dens.min()}'
    assert c_hist.min() > -1e-10, f'Negative chemical field: {c_hist.min()}'
    mass_error = np.max(np.abs(dens.mean(axis=1)-par['rho0']))
    assert mass_error < 1e-9, f'Cell mass not conserved: {mass_error}'
    print(f'[1D KS] chi={chi:.3f}; steps={nsteps}; final peak={dens[-1].max():.3f}; '
          f'mass relative error={mass_error/par["rho0"]:.2e}')
    return x, times, dens, c_hist


def make_collective_figures():
    plot_style()
    low_chi, high_chi = .15, .90
    low = simulate_collective_1d(low_chi)
    high = simulate_collective_1d(high_chi)
    x, t, dens_lo, chem_lo = low
    _, _, dens_hi, chem_hi = high
    # Figure 06: two space-time heat maps with independently scaled colorbars.
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.7), sharey=True,
                             constrained_layout=True)
    for ax, dens, chi in [(axes[0], dens_lo, low_chi),
                           (axes[1], dens_hi, high_chi)]:
        im = ax.imshow(dens/PARAMS['rho0'], origin='lower', aspect='auto',
                       extent=[0, PARAMS['L'], t[0], t[-1]],
                       cmap='viridis', interpolation='nearest')
        gamma = chi*PARAMS['alpha']*PARAMS['rho0']/(PARAMS['D']*PARAMS['kappa'])
        ax.set(xlabel='Position x', title=rf'$\chi = {chi:.2f}$, $\Gamma={gamma:.1f}$')
        cb = fig.colorbar(im, ax=ax, fraction=.046, pad=.035)
        cb.set_label(r'$\rho / \rho_0$')
    axes[0].set_ylabel('Time t')
    fig.suptitle('Collective feedback: spreading versus self-amplifying clusters',
                 fontsize=14, fontweight='semibold')
    finish(fig, '06_collective_1d_space_time')

    # Figure 07: show actual density and chemical, not just a heatmap.
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.4, 6.2), sharex=True)
    for idx, mark, lab in [(0, ':', 'Initial'),
                            (int(.45*(len(t)-1)), '--', f't = {t[int(.45*(len(t)-1))]:g}'),
                            (-1, '-', f't = {t[-1]:g}')]:
        ax1.plot(x, dens_hi[idx], mark, lw=2, label=lab)
    ax1.set(ylabel=r'Cell density $\rho$',
            title='Feedback makes an initially small fluctuation grow')
    ax1.legend(frameon=False, ncol=3)
    ax2.plot(x, dens_hi[-1], label=r'Cells $\rho(x)$', lw=2)
    ax2.plot(x, chem_hi[-1], label='Attractant $c(x)$', lw=2, ls='--')
    ax2.set(xlabel='Position x', ylabel='Field value',
            title='Cells and attractant become spatially correlated')
    ax2.legend(frameon=False)
    fig.tight_layout(h_pad=1.2)
    finish(fig, '07_collective_1d_snapshots')

    # Figure 08: fast-chemical linear growth rate, an optional teacher add-on.
    q = np.linspace(.00001, 2.4, 550)
    fig, ax = plt.subplots(figsize=(8.7, 4.3))
    p = PARAMS
    for chi, lab in [(low_chi, r'$\Gamma=0.5$ (stable)'),
                     (high_chi, r'$\Gamma=3.0$ (unstable)')]:
        sigma = q*q * (chi*p['alpha']*p['rho0'] /
                       (p['kappa']+p['Dc']*q*q) - p['D'])
        ax.plot(q, sigma, label=lab, lw=2.15)
    ax.axhline(0, color='black', lw=.9)
    q_min = 2*np.pi/p['L']
    ax.axvline(q_min, color='gray', lw=1, ls=':',
               label=r'First periodic mode $2\pi/L$')
    ax.set(xlabel='Wave number q', ylabel=r'Approximate growth rate $\sigma(q)$',
           title='Optional: when can a density fluctuation grow?', xlim=(0, 2.4))
    ax.legend(frameon=False, loc='lower left')
    fig.tight_layout()
    finish(fig, '08_optional_growth_rate')


if __name__ == '__main__':
    make_collective_figures()
