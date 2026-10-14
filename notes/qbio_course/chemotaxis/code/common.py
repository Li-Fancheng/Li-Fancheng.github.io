"""Shared plotting and numerical helpers for the chemotaxis lecture.

Only NumPy and Matplotlib are required. Figures deliberately use English
labels to avoid cross-platform issues with CJK fonts in scientific plots.
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # work on machines without a graphical display
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = ROOT / 'figures'
FIG_DIR.mkdir(parents=True, exist_ok=True)


def plot_style():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans',
        'font.size': 11,
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'axes.titleweight': 'semibold',
        'figure.titlesize': 14,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': False,
        'savefig.facecolor': 'white',
        'figure.facecolor': 'white',
        'svg.fonttype': 'none',
    })


def finish(fig, name, dpi=190, svg=False):
    fig.savefig(FIG_DIR / f'{name}.png', dpi=dpi, bbox_inches='tight', pad_inches=.16)
    if svg:
        fig.savefig(FIG_DIR / f'{name}.svg', bbox_inches='tight', pad_inches=.16)
    plt.close(fig)


def normal_pdf(x, mu, var):
    return np.exp(-0.5 * (x-mu)**2 / var) / np.sqrt(2*np.pi*var)


def advdiff_flux_1d(rho, c, dx, D, chi):
    """J_{i+1/2}: centered diffusion + positivity-preserving upwind drift."""
    rho_right = np.roll(rho, -1)
    u_face = chi * (np.roll(c, -1)-c) / dx
    return (-D*(rho_right-rho)/dx
            + np.maximum(u_face, 0)*rho
            + np.minimum(u_face, 0)*rho_right), u_face


def periodic_helmholtz_1d(rho, L, Dc, alpha, kappa):
    """Solve (kappa - Dc*d_xx)c = alpha*rho on a periodic grid."""
    dx = L/len(rho)
    k = 2*np.pi*np.fft.fftfreq(len(rho), d=dx)
    return np.fft.ifft(alpha*np.fft.fft(rho)/(kappa+Dc*k*k)).real
