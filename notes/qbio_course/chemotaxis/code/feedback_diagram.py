"""Figure 05: manually drawn, fully reproducible causal feedback diagram."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from common import plot_style, finish


def make_feedback_figure():
    plot_style()
    fig, ax = plt.subplots(figsize=(10.3, 3.4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3.30)
    ax.axis('off')
    spec = [(0.30, 1.45, 2.55, 1.05, 'More cells locally', r'$\delta \rho > 0$'),
            (3.70, 1.45, 2.55, 1.05, 'More attractant', r'$\delta c > 0$'),
            (7.10, 1.45, 2.55, 1.05, 'Inward drift', r'$J_{\rm chem} \propto \rho \nabla c$')]
    colors = ['#e6f3f6', '#e3eee7', '#ecedf9']
    for (xx, yy, ww, hh, title, formula), color in zip(spec, colors):
        patch = FancyBboxPatch((xx, yy), ww, hh,
                               boxstyle='round,pad=0.10,rounding_size=0.14',
                               linewidth=1.0, edgecolor='#56636d', facecolor=color)
        ax.add_patch(patch)
        ax.text(xx+ww/2, yy+.74, title, ha='center', va='center', fontsize=11.5, fontweight='medium')
        ax.text(xx+ww/2, yy+.27, formula, ha='center', va='center', fontsize=13)
    for xstart, xend, label in [(2.99, 3.57, 'secretion'), (6.39, 6.97, 'chemotaxis')]:
        ax.add_patch(FancyArrowPatch((xstart, 1.94), (xend, 1.94),
                                    arrowstyle='-|>', mutation_scale=16,
                                    linewidth=1.5, color='#44515c'))
        ax.text((xstart+xend)/2, 2.70, label, ha='center', fontsize=9.5)
    ax.add_patch(FancyArrowPatch((8.42, 1.25), (1.55, 1.25),
                                connectionstyle='arc3,rad=-.21',
                                arrowstyle='-|>', mutation_scale=17,
                                linewidth=1.6, color='#b24b3b'))
    ax.text(5.0, .18, 'Positive feedback amplifies a density fluctuation',
            color='#a34133', ha='center', va='center', fontsize=11,
            bbox={'fc':'white','ec':'none', 'pad':1.0})
    ax.set_title('Cells change the field they follow', pad=12)
    finish(fig, '05_field_mediated_feedback')


if __name__ == '__main__':
    make_feedback_figure()
