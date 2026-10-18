"""Figure 04: chemotaxis in an imposed chemical profile.

The prescribed c(x) does NOT depend on rho. This is directed accumulation,
not spontaneous collective pattern formation.
"""
import numpy as np
import matplotlib.pyplot as plt
from common import plot_style, finish, advdiff_flux_1d


def simulate_fixed_field():
    L, nx, D, chi = 20., 256, .25, .65
    dx = L/nx
    x = np.arange(nx)*dx
    # A truly periodic Gaussian peak (nearest periodic images suffice here).
    width = 2.0
    c = sum(np.exp(-0.5*((x-(L/2+m*L))/width)**2) for m in (-1, 0, 1))
    rho = np.ones(nx)
    t, dt_nominal = 0., .0025
    times = [0., 2., 8., 25.]
    snapshots = [rho.copy()]
    for target in times[1:]:
        while t < target-1e-10:
            flux, u = advdiff_flux_1d(rho, c, dx, D, chi)
            maxout = 2*D/dx**2+2*np.max(np.abs(u))/dx
            dt = min(dt_nominal, target-t, .8/maxout)
            rho = rho-dt*(flux-np.roll(flux, 1))/dx
            t += dt
        snapshots.append(rho.copy())
    assert np.min(rho) >= -1e-12, 'Negative density: check CFL'
    assert abs(rho.mean()-1) < 1e-10, 'Mass conservation failed'
    return x, c, times, snapshots, D, chi


def make_fixed_field_figure():
    plot_style()
    x, c, times, snapshots, D, chi = simulate_fixed_field()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.0, 6.3), sharex=True,
                                    gridspec_kw={'height_ratios':[.9, 1.5]})
    ax1.fill_between(x, 0, c, color='tab:green', alpha=.18)
    ax1.plot(x, c, color='tab:green', lw=2.2)
    ax1.set(ylabel='Chemical $c(x)$', title='Chemotaxis toward a fixed attractant peak')
    ax1.annotate('Field imposed externally', xy=(10., 1.), xytext=(13., .75),
                 arrowprops={'arrowstyle':'->', 'lw':1.1}, fontsize=10)
    for tm, rr in zip(times, snapshots):
        ax2.plot(x, rr, lw=2, label=rf'$t={tm:g}$')
    stationary = np.exp(chi/D*(c-c.max()))
    stationary /= stationary.mean()
    ax2.plot(x, stationary, '--', color='black', alpha=.7, lw=1.5,
             label=r'Zero-flux shape $\propto e^{\chi c/D}$')
    ax2.set(xlabel='Position $x$', ylabel='Cell density $\\rho(x,t)$')
    ax2.legend(frameon=False, ncol=3)
    ax2.set_xlim(x[0], x[-1])
    fig.tight_layout(h_pad=1.2)
    finish(fig, '04_chemotaxis_in_fixed_field')
    print(f'[FIXED] mean cell density={snapshots[-1].mean():.7f}; min={snapshots[-1].min():.5f}')


if __name__ == '__main__':
    make_fixed_field_figure()
