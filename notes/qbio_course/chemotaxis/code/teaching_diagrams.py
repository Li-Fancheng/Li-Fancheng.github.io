"""Reproducible diagrams A--E for the step-by-step derivations.

Run: python code/teaching_diagrams.py
Figure labels are in English.
SVG text is converted to paths so readers do not need that font installed.
"""
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
import numpy as np

from common import plot_style, finish


BLUE = '#2563a6'
GREEN = '#218575'
ORANGE = '#cc7040'
INK = '#26384a'
GRAY = '#738294'
def setup():
    plot_style()
    plt.rcParams.update({
        'font.family': 'DejaVu Sans',
        'font.size': 12,
        'axes.unicode_minus': False,
        'axes.titleweight': 'normal',
        'text.color': INK,
        'axes.labelcolor': INK,
        'axes.titlecolor': INK,
        'svg.fonttype': 'path',
    })


def arrow(ax, start, end, color=BLUE, width=2, scale=16, style='-|>'):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle=style,
                                mutation_scale=scale, linewidth=width,
                                color=color, shrinkA=4, shrinkB=4))


def diagram_axes(ax, xlim, ylim):
    ax.set(xlim=xlim, ylim=ylim)
    ax.axis('off')


def node(ax, xy, label, color=BLUE, radius=.22):
    ax.add_patch(Circle(xy, radius, facecolor='white', edgecolor=color,
                        linewidth=2, zorder=3))
    ax.text(*xy, label, ha='center', va='center', fontsize=13, color=color,
            zorder=4)


def make_walk_bookkeeping():
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.6))
    ax = axes[0]
    diagram_axes(ax, (-3.15, 2.65), (-.85, 3.5))
    levels = [(2.6, [0], ['1']), (1.35, [-1, 1], [r'$1/2$', r'$1/2$']),
              (.1, [-2, 0, 2], [r'$1/4$', r'$1/2$', r'$1/4$'])]
    for (y, xs, _), (yn, xn, _) in zip(levels[:-1], levels[1:]):
        for x in xs:
            for dest in (x-1, x+1):
                arrow(ax, (x, y-.22), (dest, yn+.22), GRAY, 1.5)
    for n, (y, xs, labels) in enumerate(levels):
        ax.text(-3.0, y, f'Step {n}', va='center', fontsize=11)
        for x, label in zip(xs, labels):
            node(ax, (x, y), label, GREEN if n == 2 and x == 0 else BLUE,
                 radius=.29)
    for x, label in [(-2, r'$-2a$'), (0, '$0$'), (2, '$2a$')]:
        ax.text(x, -.49, label, ha='center', fontsize=12)
    ax.set_title('Cell motion without a gradient', pad=14)

    ax = axes[1]
    diagram_axes(ax, (-2.8, 2.8), (-.85, 3.5))
    for x, lab in [(-1.7, r'$P_{i-1}^{n}$'), (1.7, r'$P_{i+1}^{n}$')]:
        ax.text(x, 2.7, lab, ha='center', va='center', fontsize=19, color=BLUE)
        arrow(ax, (x, 2.25), (0, 1.12), BLUE, 2)
        ax.text(x*.80, 1.55, r'$\times\,1/2$', ha='center', fontsize=13)
    node(ax, (0, .8), r'$i$', GREEN, radius=.32)
    ax.text(0, .1, r'$P_i^{n+1}=\frac{1}{2} P_{i-1}^{n}+\frac{1}{2} P_{i+1}^{n}$',
            ha='center', fontsize=17)
    ax.text(0, -.55, 'Add the two incoming contributions', ha='center', fontsize=12)
    ax.set_title('Incoming probabilities', pad=14)
    fig.subplots_adjust(left=.03, right=.98, bottom=.09, top=.84, wspace=.18)
    finish(fig, 'A_walk_probability_bookkeeping', svg=True)


def make_flux_balance():
    fig, ax = plt.subplots(figsize=(10.8, 4.6))
    diagram_axes(ax, (0, 10.8), (-.4, 4.3))
    for x0, count, color in [(1.0, 12, BLUE), (6.4, 4, ORANGE)]:
        ax.add_patch(FancyBboxPatch((x0, 1.65), 3.4, 1.7,
                                   boxstyle='round,pad=.04,rounding_size=.12',
                                   facecolor=color+'10', edgecolor=color, linewidth=1.4))
        for k in range(count):
            col, row = k % 6, k // 6
            ax.add_patch(Circle((x0+.36+.50*col, 2.9-.55*row), .10,
                                color=color))
        ax.text(x0+1.7, 3.65, f'{count} cells',
                ha='center', fontsize=15, color=color)
    ax.plot([5.4, 5.4], [1.65, 3.9], '--', color=GRAY, linewidth=1.3)
    ax.text(5.4, 4.03, 'Boundary', ha='center', fontsize=11)
    arrow(ax, (3.0, 1.1), (7.8, 1.1), BLUE, width=4)
    ax.text(5.4, 1.35, 'Right: mean 12 × 1/2 = 6', ha='center', color=BLUE)
    arrow(ax, (7.8, .45), (3.0, .45), ORANGE, width=2)
    ax.text(5.4, .02, 'Left: mean 4 × 1/2 = 2', ha='center', color=ORANGE)
    ax.text(5.4, -.34, 'Net rightward: mean 6 − 2 = 4 cells per step',
            ha='center', fontweight='normal', fontsize=14)
    ax.set_title('Random migration produces a net flux', pad=16)
    fig.subplots_adjust(left=.02, right=.98, bottom=.08, top=.86)
    finish(fig, 'B_diffusion_boundary_flux', svg=True)


def make_gradient_bias():
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    cases = [(.30, .30, .5, .5), (.15, .75, .35, .65)]
    for ax, (cleft, cright, pleft, pright) in zip(axes, cases):
        diagram_axes(ax, (-2.8, 2.8), (-.5, 3.5))
        for x, c in [(-1.8, cleft), (1.8, cright)]:
            ax.add_patch(FancyBboxPatch((x-.45, 1.4), .9, 1.25,
                                       boxstyle='round,pad=.03',
                                       facecolor=GREEN+'12', edgecolor=GREEN))
            ax.add_patch(plt.Rectangle((x-.43, 1.41), .86, 1.2*c,
                                       facecolor=GREEN, alpha=.6))
        ax.text(-1.8, 2.92, r'$c(x-a)$', ha='center', fontsize=14)
        ax.text(1.8, 2.92, r'$c(x+a)$', ha='center', fontsize=14)
        node(ax, (0, 1.95), 'cell', BLUE, radius=.34)
        arrow(ax, (-.4, 1.95), (-1.20, 1.95), BLUE, width=5*pleft)
        arrow(ax, (.4, 1.95), (1.20, 1.95), BLUE, width=5*pright)
        ax.text(-1., 1.05, rf'$p_-={pleft:.2f}$', ha='center', fontsize=15)
        ax.text(1., 1.05, rf'$p_+={pright:.2f}$', ha='center', fontsize=15)
        ax.text(0, .42, r'$v=a(p_+-p_-)/\tau$', ha='center', fontsize=15)
        ax.text(0, -.10, 'Zero mean drift' if pleft == pright
                else 'Rightward bias; left jumps still occur', ha='center', fontsize=12)
        ax.set_title('Equal concentrations' if pleft == pright
                     else 'Higher concentration on the right', pad=14)
    fig.subplots_adjust(left=.03, right=.98, bottom=.08, top=.84, wspace=.14)
    finish(fig, 'C_gradient_to_jump_bias', svg=True)


def make_mode_growth():
    fig, axes = plt.subplots(1, 2, figsize=(11.7, 4.6))
    x = np.linspace(0, 20, 500)
    q = 2*np.pi/10
    ax = axes[0]
    ax.plot(x, 1+.10*np.cos(q*x), color=BLUE, linewidth=2.5)
    ax.axhline(1, color=GRAY, linewidth=1, linestyle='--')
    ax.annotate('', xy=(10, 1.14), xytext=(0, 1.14),
                arrowprops={'arrowstyle':'<->', 'color':INK, 'linewidth':1.3})
    ax.text(5, 1.16, r'$\lambda=2\pi/q$', ha='center', fontsize=13)
    ax.annotate('', xy=(10, 1.10), xytext=(10, 1),
                arrowprops={'arrowstyle':'<->', 'color':ORANGE, 'linewidth':1.5})
    ax.text(11.0, 1.115, '$A(t)$', color=ORANGE, fontsize=13)
    ax.text(18.5, 1.005, r'$\rho_0$', va='bottom', color=GRAY)
    ax.set(xlim=(-.8, 20.8), ylim=(.86, 1.22), xlabel='Position x',
           ylabel='Cell density', xticks=[0, 10, 20])
    ax.set_title('Space: wavelength versus amplitude', pad=15)

    ax = axes[1]
    times = np.linspace(0, 12, 200)
    for chi, color in [(.15, BLUE), (.90, ORANGE)]:
        u, w, h = .3*q*q, 1+q*q, chi*q*q
        sigma = -(u+w)/2 + np.sqrt((u-w)**2+4*h)/2
        label = ('Weak feedback' if chi == .15
                 else 'Strong feedback')
        ax.plot(times, .01*np.exp(sigma*times), color=color, linewidth=2.5,
                label=label+rf' ($\chi={chi:.2f}$)')
    ax.set(xlabel='Time t', ylabel='Amplitude A(t)',
           xlim=(0, 12), ylim=(0, .051))
    ax.set_title('Decay or growth of a density mode', pad=15)
    ax.legend(frameon=False, loc='upper left', fontsize=11)
    ax.text(.05, .38, r'$A(t)=A(0)e^{\sigma_+t}$', transform=ax.transAxes, fontsize=14)
    fig.subplots_adjust(left=.075, right=.98, bottom=.17, top=.83, wspace=.32)
    finish(fig, 'D_cosine_mode_and_growth', svg=True)


def make_signal_kernel():
    fig, axes = plt.subplots(1, 2, figsize=(11.7, 4.6))
    x = np.linspace(-5, 5, 601)
    ell = 1.3
    kernel = np.exp(-np.abs(x)/ell)/(2*ell)  # kappa=1, Dc=ell**2
    ax = axes[0]
    ax.plot(x, kernel, color=GREEN, linewidth=2.5)
    ax.fill_between(x, kernel, color=GREEN, alpha=.09)
    ax.axvline(ell, color=GRAY, linewidth=1, linestyle=':')
    ax.scatter([ell], [1/(2*ell*np.e)], color=ORANGE, zorder=5)
    ax.annotate(r'$G(\ell_c)=G(0)/e$', xy=(ell, 1/(2*ell*np.e)),
                xytext=(1.65, .27), fontsize=12,
                arrowprops={'arrowstyle':'->', 'color':GRAY})
    ax.text(0, .405, 'Unit point source', ha='center', fontsize=11)
    ax.set(xlim=(-5, 5), ylim=(0, .47), xlabel='Position relative to source x',
           ylabel='Steady signal G(x)',
           xticks=[-ell, 0, ell], xticklabels=[r'$-\ell_c$', '0', r'$\ell_c$'])
    ax.set_title('One source: signal decays with distance', pad=15)

    ax = axes[1]
    gleft = np.exp(-np.abs(x+1.8)/ell)/(2*ell)
    gright = np.exp(-np.abs(x-1.5)/ell)/(2*ell)
    ax.plot(x, gleft, '--', color=BLUE, linewidth=1.8, label='Left source')
    ax.plot(x, gright, '--', color=ORANGE, linewidth=1.8, label='Right source')
    ax.plot(x, gleft+gright, color=GREEN, linewidth=2.7, label='Sum')
    for pos, color in [(-1.8, BLUE), (1.5, ORANGE)]:
        ax.scatter([pos], [.018], s=55, color=color, zorder=5)
    # Follow the slope of the total signal toward each local maximum.
    arrow(ax, (-3.8, .52), (-2.1, .52), BLUE, 1.8)
    arrow(ax, (3.8, .52), (1.8, .52), ORANGE, 1.8)
    ax.set(xlim=(-5, 5), ylim=(0, .60), xlabel='Position x',
           ylabel='Total signal concentration')
    ax.legend(frameon=False, loc='upper center', fontsize=10)
    ax.set_title('Multiple sources: add their contributions', pad=15)
    fig.subplots_adjust(left=.075, right=.98, bottom=.17, top=.83, wspace=.32)
    finish(fig, 'E_signal_range_and_superposition', svg=True)


def make_teaching_diagrams():
    setup()
    for make in (make_walk_bookkeeping, make_flux_balance, make_gradient_bias,
                 make_mode_growth, make_signal_kernel):
        make()
    print('Teaching diagrams A--E saved as PNG and SVG.')


if __name__ == '__main__':
    make_teaching_diagrams()
