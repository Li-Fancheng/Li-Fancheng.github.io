"""Teaching illustration of temporal comparison; no experimental fit or data.

Run from the chemotaxis directory: python code/chemosensing.py
The paper's original figure is kept separately and is never overwritten here.
"""
import numpy as np
import matplotlib.pyplot as plt
from common import plot_style, finish


def make_chemosensing_figure():
    plot_style()
    gain, tau_fast, tau_memory = 1.0, 0.5, 3.0
    factor = gain * tau_memory / (tau_memory - tau_fast)
    time = np.linspace(0, 18, 1801)
    elapsed = np.maximum(time - 2.0, 0.0)
    stimulus = (time >= 2.0).astype(float)
    memory = stimulus * (1 - np.exp(-elapsed / tau_memory))
    response = factor * (np.exp(-elapsed / tau_memory) - np.exp(-elapsed / tau_fast))
    kernel = factor * (np.exp(-time / tau_fast) / tau_fast
                       - np.exp(-time / tau_memory) / tau_memory)

    fig, axes = plt.subplots(2, 1, figsize=(6.2, 7.6), layout='constrained')
    ax = axes[0]
    ax.plot(time, stimulus, '--', color='#64748b', label=r'Input $q$')
    ax.plot(time, memory, color='#c78021', label=r'Memory $m$')
    ax.plot(time, response, color='#087e8b', label=r'Response $y$')
    ax.set(xlim=(0, 18), ylim=(-0.04, 1.12), xlabel='Time (s)', ylabel='Normalized level',
           title='A  Input, memory, and recovery')
    ax.legend(loc='center right', bbox_to_anchor=(1, 0.42), fontsize=16)
    ax.axvline(2, color='#cbd5e1', linewidth=1, zorder=0)

    ax = axes[1]
    ax.axhline(0, color='#64748b', linewidth=1)
    ax.fill_between(time, 0, kernel, where=kernel>=0, color='#5bb3bd', alpha=.35)
    ax.fill_between(time, 0, kernel, where=kernel<0, color='#e8ae65', alpha=.55)
    ax.plot(time, kernel, color='#243747')
    ax.set(xlim=(0, 18), ylim=(-.45, 2.18), xlabel='Time before the present (s)',
           ylabel=r'Weight $R(t)$ (s$^{-1}$)', title='B  Positive and negative weights')
    ax.text(1.3, 1.3, 'Recent', color='#087e8b')
    ax.text(5.4, -.32, 'Older', color='#986017')
    ax.text(7.3, .9, r'$\int_0^\infty R(t)\,dt=0$', fontsize=18,
            bbox=dict(facecolor='white', edgecolor='none', alpha=.95))
    for ax in axes:
        ax.tick_params(labelsize=16)
        ax.xaxis.label.set_fontsize(16)
        ax.yaxis.label.set_fontsize(16)
    finish(fig, 'sensing_01_temporal_filter')


if __name__ == '__main__':
    make_chemosensing_figure()
