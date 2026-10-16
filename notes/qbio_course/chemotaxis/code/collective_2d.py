"""Asset 09 (article Fig. 18): minimal Keller--Segel on a 2D torus.

A finite-volume update of cell transport uses conservative face fluxes
and upwinding for chemotactic drift. The chemical field evolves via a
reaction-diffusion equation. This is a toy model of pattern formation;
no quantitative comparison to E. coli experiments is implied.
"""
import numpy as np
import matplotlib.pyplot as plt
from common import plot_style, finish


def initial_state_2d(n, L, Dc, alpha, kappa):
    rng = np.random.default_rng(19)
    dx = L/n
    noise = rng.normal(size=(n, n))
    k = 2*np.pi*np.fft.fftfreq(n, dx)
    ky, kx = np.meshgrid(k, k, indexing='ij')
    filt = np.exp(-(kx*kx+ky*ky)/(2*1.35**2))
    smooth = np.fft.ifft2(np.fft.fft2(noise)*filt).real
    smooth -= smooth.mean()
    smooth /= smooth.std()
    rho = 1 + .055*smooth
    c = np.fft.ifft2(alpha*np.fft.fft2(rho)/(kappa+Dc*(kx*kx+ky*ky))).real
    return rho, c


def simulate_2d(n=96, L=24., D=.3, chi=.9, Dc=1., alpha=1., kappa=1.,
                final_time=19.5, nominal_dt=.004):
    dx = L/n
    rho, c = initial_state_2d(n, L, Dc, alpha, kappa)
    saved = {0: rho.copy()}
    t, next_frames = 0., [13., final_time]
    c_limit = .78/(4*Dc/dx**2+kappa)
    steps = 0
    for target in next_frames:
        while t < target - 1e-10:
            xr = np.roll(rho, -1, axis=1)
            yd = np.roll(rho, -1, axis=0)
            ux = chi*(np.roll(c, -1, axis=1)-c)/dx
            uy = chi*(np.roll(c, -1, axis=0)-c)/dx
            Jx = -D*(xr-rho)/dx + np.maximum(ux, 0)*rho + np.minimum(ux, 0)*xr
            Jy = -D*(yd-rho)/dx + np.maximum(uy, 0)*rho + np.minimum(uy, 0)*yd
            drift_bound = 2*(np.max(np.abs(ux))+np.max(np.abs(uy)))/dx
            rho_limit = .80/(4*D/dx**2+drift_bound+1e-12)
            dt = min(nominal_dt, rho_limit, c_limit, target-t)
            rho_next = rho - dt/dx*(Jx-np.roll(Jx, 1, axis=1)
                                    +Jy-np.roll(Jy, 1, axis=0))
            lapc = (np.roll(c, 1, 0)+np.roll(c, -1, 0)
                    +np.roll(c, 1, 1)+np.roll(c, -1, 1)-4*c)/dx**2
            c_next = c+dt*(Dc*lapc+alpha*rho-kappa*c)
            rho, c = rho_next, c_next
            t += dt
            steps += 1
        saved[target] = rho.copy()
    assert min(z.min() for z in saved.values()) >= -1e-10
    assert abs(rho.mean()-1) < 1e-10
    print(f'[2D KS] {steps} steps, final min={rho.min():.3f}, '
          f'max={rho.max():.3f}, mean={rho.mean():.8f}')
    return saved, L


def make_2d_figure():
    plot_style()
    saved, L = simulate_2d()
    fig, axes = plt.subplots(1, 3, figsize=(11.1, 3.75), constrained_layout=True)
    times = list(saved)
    for ax, tm in zip(axes, times):
        arr = saved[tm]
        im = ax.imshow(arr, origin='lower', cmap='viridis',
                       extent=(0, L, 0, L), interpolation='nearest')
        ax.set(title=f't = {tm:g}', xlabel='x', ylabel='y', aspect='equal')
        cb = fig.colorbar(im, ax=ax, fraction=.046, pad=.03)
        cb.set_label(r'$\rho(x,y)$')
    fig.suptitle('Two-dimensional chemotaxis: early density aggregation',
                 fontsize=14, fontweight='semibold')
    finish(fig, '09_optional_2d_patterns', dpi=180)


if __name__ == '__main__':
    make_2d_figure()
