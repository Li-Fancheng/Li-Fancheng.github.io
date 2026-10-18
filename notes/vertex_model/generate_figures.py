"""Reproduce the 25 teaching figures; no fitted tissue data.

Run from any directory with NumPy and Matplotlib installed. Each figure caption
in the corresponding note states its formula, parameters, and limitations.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle, Polygon, FancyArrowPatch, Ellipse, Rectangle, Arc
from matplotlib.text import Text
import json
import numpy as np

ROOT = Path(__file__).resolve().parent
POLICIES = json.loads((ROOT.parents[1] / "data/vertex_figures.json").read_text())
BLUE, ORANGE, GREEN, INK = "#277da8", "#d67837", "#429582", "#263e50"

for font in ("Arial.ttf", "Arial Bold.ttf"):
    path = Path("/System/Library/Fonts/Supplemental") / font
    if path.exists():
        font_manager.fontManager.addfont(path)
available = {f.name for f in font_manager.fontManager.ttflist}
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": [f for f in ("Arial", "Liberation Sans", "DejaVu Sans") if f in available],
    "font.size": 13,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "legend.fontsize": 11,
    "mathtext.fontset": "cm",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.labelcolor": INK,
    "text.color": INK,
    "axes.titlecolor": INK,
    "lines.linewidth": 2.2,
    "svg.fonttype": "path",
    "svg.hashsalt": "vertex-model-notes",
    "savefig.facecolor": "white",
})


def save(fig, slug, name):
    directory = ROOT / slug / "figures"
    directory.mkdir(parents=True, exist_ok=True)
    # Size labels at the intended web width, including dense multi-panel plots.
    width = POLICIES[f"{slug}/figures/{name}"]["width"]
    label_floor = 13 * fig.get_figwidth() * 72 / width
    for label in fig.findobj(Text):
        if label.get_visible() and label.get_text():
            label.set_fontsize(max(label.get_fontsize(), label_floor))
    fig.savefig(directory / f"{name}.svg", bbox_inches="tight", pad_inches=.2,
                metadata={"Date": None})
    fig.savefig(directory / f"{name}.png", dpi=300, bbox_inches="tight", pad_inches=.2)
    plt.close(fig)


def arrow(ax, start, end, color=INK):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=15,
                                linewidth=1.8, color=color))


def model():
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), layout="constrained")
    ax = axes[0]
    r = .8
    theta = np.arange(6)*np.pi/3
    points = np.c_[r*np.cos(theta), r*np.sin(theta)]
    for offset, color in [((1.5*r, np.sqrt(3)*r/2), "#edf2f5"),
                          ((0, np.sqrt(3)*r), "#edf2f5"), ((0, 0), "#d8edf6")]:
        ax.add_patch(Polygon(points + offset, facecolor=color, edgecolor=INK, lw=1.5))
    ax.scatter(*points.T, s=22, color=BLUE, zorder=3)
    ax.scatter(*points[1], s=45, color=ORANGE, zorder=4)
    ax.annotate(r"$\mathbf{r}_v$", xy=points[1], xytext=(1.25, 1.3),
                arrowprops={"arrowstyle": "->", "color": ORANGE}, fontsize=16)
    ax.text(0, 0, r"$A_i,\ P_i$", ha="center", va="center", fontsize=17)
    ax.set(xlim=(-1.1, 2.4), ylim=(-1, 2.4), aspect="equal", title="A  Vertex coordinates")
    ax.axis("off")
    p = np.linspace(3.2, 4.5, 200)
    axes[1].plot(p, .5*(p-3.8)**2, color=BLUE)
    axes[1].axvline(3.8, color=ORANGE, ls="--", lw=1.3)
    axes[1].set(xlabel=r"$p=P/\sqrt{A_0}$", ylabel=r"$E_P/(K_A A_0^2)$",
                title="B  Preferred perimeter")
    n = np.arange(4, 13)
    pn = np.sqrt(4*n*np.tan(np.pi/n))
    axes[2].plot(n, pn, "o-", color=GREEN, markersize=5)
    axes[2].axhline(np.sqrt(4*np.pi), color=INK, ls=":", lw=1.3)
    axes[2].text(8, 3.57, "circle limit", fontsize=11)
    axes[2].set(xlabel="Number of sides, n", ylabel=r"$p_n=P_n/\sqrt{A}$",
                title="C  Polygonal lower bound", xticks=[4, 6, 8, 10, 12])
    save(fig, "foundations", "model")


def response():
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout="constrained")
    gamma = np.linspace(-.4, .4, 200)
    axes[0].plot(gamma, .6*gamma**2, color=BLUE, label="Affine: u = 0")
    axes[0].plot(gamma, .1*gamma**2, color=ORANGE, label=r"Relaxed: $u=\gamma$")
    axes[0].set(xlabel=r"Shear strain, $\gamma$", ylabel="Energy", title="A  Internal relaxation softens shear")
    axes[0].legend()
    t = np.linspace(0, 7, 200)
    axes[1].plot(t, .2+np.exp(-t), color=BLUE)
    axes[1].axhline(.2, color=ORANGE, ls="--", label="Relaxed modulus = 0.2")
    axes[1].set(xlabel=r"Time, $t/\tau$", ylabel=r"$\sigma(t)/\gamma_0$",
                title="B  Fixed-topology stress relaxation", ylim=(0, 1.3))
    axes[1].legend()
    save(fig, "linear_response", "response")


def incompatibility():
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), layout="constrained")
    x = np.linspace(-1.05, 1.05, 300)
    for ax, delta, color, title in zip(axes, [-.35, 0, .35], [BLUE, ORANGE, GREEN],
                                      ["Incompatible target", "Quartic minimum", "Compatible target"]):
        ax.plot(x, .5*(x*x-delta)**2, color=color)
        minima = [0] if delta <= 0 else [-np.sqrt(delta), np.sqrt(delta)]
        ax.scatter(minima, [.5*(m*m-delta)**2 for m in minima], color=color, zorder=3)
        ax.set(xlabel="Coordinate, x", ylabel="Energy", title=title, ylim=(-.02, 1.1))
        ax.text(.5, .9, rf"$\delta={delta}$", transform=ax.transAxes, ha="center")
    save(fig, "rigidity", "incompatibility")


def t1():
    fig = plt.figure(figsize=(10, 6), layout="constrained")
    grid = fig.add_gridspec(2, 3, height_ratios=[1.2, 1])
    outside = np.array([[-.9, .8], [.9, .8], [-.9, -.8], [.9, -.8]])
    titles = ["1  Old L-R contact", "2  Fourfold vertex", "3  New U-D contact"]
    for j, title in enumerate(titles):
        ax = fig.add_subplot(grid[0, j])
        if j == 0:
            ends = np.array([[0, .3], [0, -.3]])
            links = [ends[0], ends[0], ends[1], ends[1]]
        elif j == 1:
            ends = np.array([[0, 0], [0, 0]])
            links = [ends[0]]*4
        else:
            ends = np.array([[-.3, 0], [.3, 0]])
            links = [ends[0], ends[1], ends[0], ends[1]]
        for outer, inner in zip(outside, links):
            ax.plot([outer[0], inner[0]], [outer[1], inner[1]], color=INK, lw=2)
        ax.plot(*ends.T, color=ORANGE, lw=4)
        ax.scatter(*ends.T, s=30, color=ORANGE, zorder=3)
        for label, pos in [("L", (-.7, 0)), ("R", (.7, 0)), ("U", (0, .65)), ("D", (0, -.65))]:
            ax.text(*pos, label, ha="center", va="center", fontsize=16, color=BLUE)
        ax.set(xlim=(-1, 1), ylim=(-1, 1), aspect="equal", title=title)
        ax.axis("off")
    ax = fig.add_subplot(grid[1, :])
    s = np.linspace(-1.3, 1.3, 400)
    for h, color in [(0, BLUE), (.35, ORANGE)]:
        ax.plot(s, (1-s*s)**2-h*s, color=color, label=rf"$h={h}$")
    ax.set(xlabel="Reaction coordinate, s (old < 0; new > 0)", ylabel="Toy energy",
           title=r"A separate barrier illustration: $E(s)=(1-s^2)^2-hs$")
    ax.legend(loc="upper right")
    save(fig, "rearrangements", "t1")


def rheology():
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), layout="constrained")
    t = np.linspace(0, 12, 240)
    for tau, color in [(.5, BLUE), (2, ORANGE), (8, GREEN)]:
        axes[0].plot(t, np.exp(-t/tau), label=rf"$\tau={tau}$", color=color)
    axes[0].set(xlabel="Time, t", ylabel=r"$G(t)/G_0$", title="A  Maxwell stress relaxation")
    axes[0].legend()
    z = np.logspace(-2, 2, 200)
    axes[1].loglog(z, z*z/(1+z*z), color=BLUE, label=r"$G'/G_0$")
    axes[1].loglog(z, z/(1+z*z), color=ORANGE, label=r"$G''/G_0$")
    axes[1].set(xlabel=r"Frequency, $\omega\tau$", ylabel="Modulus ratio", title="B  Storage and loss", ylim=(1e-4, 2))
    axes[1].legend(loc="lower right")
    ax = axes[2]
    ax.set(xlim=(0, 1), ylim=(0, 1), title="C  Active persistence")
    ax.axis("off")
    arrow(ax, (.08, .55), (.92, .55), BLUE)
    ax.text(.5, .72, r"$\tau_p/\tau_{\rm mech}$", ha="center", fontsize=22)
    ax.text(.12, .38, "Short\nmemory\n< 1", ha="left", va="top")
    ax.text(.88, .38, "Persistent\ndrive\n> 1", ha="right", va="top")
    ax.text(.5, .07, "Timescale comparison; no phase boundary", ha="center", fontsize=10)
    save(fig, "rheology_activity", "rheology")


def dimensions():
    fig, axes = plt.subplots(1, 3, figsize=(11, 4), layout="constrained")
    ax = axes[0]
    ax.add_patch(Circle((0, 0), 1, ec=BLUE, fc="#edf6f9", lw=2))
    t = np.linspace(0, 2*np.pi, 300)
    for a in [.35, .7]:
        ax.plot(a*np.sin(t), np.cos(t), color=BLUE, lw=1)
        ax.plot(np.sin(t), a*np.cos(t), color=BLUE, lw=1)
    ax.text(0, -1.35, "Cells live on a surface", ha="center")
    ax.set(title="A  Curved 2D sheet")
    ax = axes[1]
    top = np.array([[-.85, .35], [.25, .35], [.85, .8], [-.25, .8]])
    bottom = top + [0, -1.2]
    ax.add_patch(Polygon(bottom, fc="#dcece6", ec=GREEN, lw=1.7))
    for j in [0, 1, 2]:
        k = (j+1) % 4
        ax.add_patch(Polygon([top[j], top[k], bottom[k], bottom[j]], fc="#edf0f2", ec=INK, lw=1.3))
    ax.add_patch(Polygon(top, fc="#f7ddc8", ec=ORANGE, lw=1.7))
    ax.text(-.65, 1.03, "Apical", color=ORANGE)
    ax.annotate("Lateral", xy=(.6, -.2), xytext=(.75, -.65), fontsize=12,
                ha="center", arrowprops={"arrowstyle": "->", "color": INK})
    ax.text(-.65, -1.18, "Basal", color=GREEN)
    ax.set(title="B  Finite-thickness monolayer")
    ax = axes[2]
    # A projected polyhedron indicates a bulk cell, not a space-filling packing.
    a = np.array([[-.8, -.6], [.3, -.8], [.9, -.05], [.35, .8], [-.65, .6]])
    center = np.array([.0, .05])
    ax.add_patch(Polygon(a, fc="#edf6f9", ec=BLUE, lw=2))
    for point in a:
        ax.plot([point[0], center[0]], [point[1], center[1]], color=BLUE, lw=1.5)
    ax.text(0, -1.35, r"Volume $V_i$, surface area $S_i$", ha="center")
    ax.set(title="C  Bulk 3D tissue")
    for ax in axes:
        ax.set(xlim=(-1.55, 1.55), ylim=(-1.65, 1.5), aspect="equal")
        ax.axis("off")
    save(fig, "geometry_3d", "dimensions")


def inference():
    fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
    ax.set(xlim=(0, 10), ylim=(0, 5))
    ax.axis("off")
    box = {"boxstyle": "round,pad=.6", "facecolor": "#edf6f9", "edgecolor": BLUE}
    ax.text(5, 4.45, "Geometry + cell trajectories", ha="center", va="center", bbox=box, fontsize=15)
    ax.text(5, 3.2, "Cell deformation / rearrangements / boundary motion",
            ha="center", va="center", bbox=box, fontsize=14)
    arrow(ax, (5, 4.1), (5, 3.6))
    hypotheses = [(1.65, "Junction tension"), (5, "Protrusive drive"), (8.35, "External traction")]
    for x, label in hypotheses:
        arrow(ax, (5, 2.8), (x, 2.2))
        ax.text(x, 1.9, label, ha="center", va="center", fontsize=13,
                bbox={**box, "facecolor": "#fbf0e6", "edgecolor": ORANGE})
        arrow(ax, (x, 1.55), (5, .95))
    ax.text(5, .55, "Predict a new perturbation response; compare with data",
            ha="center", va="center", fontsize=14,
            bbox={**box, "facecolor": "#e8f3ed", "edgecolor": GREEN})
    save(fig, "morphogenesis", "inference")


def force_balance():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.1), layout="constrained")
    ends = np.array([[0, 1], [-np.sqrt(3)/2, -.5], [np.sqrt(3)/2, -.5]])
    for j, ax in enumerate(axes):
        for k, end in enumerate(ends):
            removed = j == 1 and k == 0
            ax.plot([0, end[0]], [0, end[1]], color=ORANGE if k == 0 else INK,
                    lw=2, ls="--" if removed else "-")
            if not removed:
                arrow(ax, (0, 0), tuple(.65*end), ORANGE if k == 0 else BLUE)
        ax.scatter([0], [0], s=55, color=INK, zorder=4)
        if j == 0:
            ax.text(.13, .7, r"$T_e$", color=ORANGE, fontsize=18)
            ax.text(0, -1.05, "Balanced tension vectors", ha="center")
        else:
            arrow(ax, (.1, -.05), (.1, -.85), GREEN)
            ax.text(.23, -.82, r"$\zeta\mathbf{v}_{0+}$", color=GREEN, fontsize=17)
            ax.text(0, -1.05, "One tension removed", ha="center")
        ax.set(xlim=(-1.35, 1.35), ylim=(-1.25, 1.25), aspect="equal",
               title=["A  Before cutting", "B  Initial recoil"][j])
        ax.axis("off")
    save(fig, "foundations", "force_balance")


def topology():
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout="constrained")
    for ax, sides, title, subtitle in [
        (axes[0], [5, 6, 7], "A  Periodic plane", r"$\sum_i(6-n_i)=0$"),
        (axes[1], [5], "B  Closed sphere", r"$\sum_i(6-n_i)=12$")]:
        for j, n in enumerate(sides):
            a = np.arange(n)*2*np.pi/n+np.pi/2
            center = np.array([1.4 if len(sides)==1 else j*1.4, 0])
            points = .5*np.c_[np.cos(a), np.sin(a)]+center
            ax.add_patch(Polygon(points, fc={5:"#f7ddc8",6:"#edf2f5",7:"#d8edf6"}[n], ec=INK))
            ax.text(*center, f"{n}", ha="center", va="center", fontsize=16)
            ax.text(center[0], -.8, f"charge {6-n:+d}", ha="center", fontsize=12)
        ax.text(1.4, -1.4, subtitle, ha="center", fontsize=20)
        if ax is axes[1]:
            ax.text(1.4, 1.02, "Twelve pentagons if all others are hexagons", ha="center", fontsize=11)
            ax.text(2.6, 0, r"$\times12$", ha="center", va="center", fontsize=19)
        ax.set(xlim=(-.8, 3.6), ylim=(-1.7, 1.3), aspect="equal", title=title)
        ax.axis("off")
    save(fig, "foundations", "topology")


def affine_hexagon():
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7), layout="constrained")
    radius = np.sqrt(2/(3*np.sqrt(3)))
    angles = np.arange(6)*np.pi/3
    vertices = radius*np.c_[np.cos(angles), np.sin(angles)]
    for gamma, color in [(-.5, BLUE), (0, INK), (.5, ORANGE)]:
        xy = vertices @ np.array([[1, 0], [gamma, 1]])
        xy = np.vstack([xy, xy[0]])
        axes[0].plot(*xy.T, color=color, label=rf"$\gamma={gamma}$")
    axes[0].set(aspect="equal", title="A  Area-preserving shear")
    axes[0].axis("off"); axes[0].legend(loc="upper center", bbox_to_anchor=(.5, .03), ncol=3, fontsize=10)
    gamma = np.linspace(-.5, .5, 201)
    perimeter = []
    for g in gamma:
        xy = vertices @ np.array([[1, 0], [g, 1]])
        perimeter.append(np.linalg.norm(np.roll(xy,-1,axis=0)-xy,axis=1).sum())
    p6 = np.sqrt(8*np.sqrt(3))
    axes[1].plot(gamma, perimeter, color=BLUE)
    axes[1].plot(gamma, p6+3*p6/16*gamma**2, "--", color=ORANGE, label="Quadratic expansion")
    axes[1].set(xlabel=r"$\gamma$", ylabel="Perimeter, P", title="B  Perimeter along the path")
    axes[1].legend(loc="upper center", fontsize=10)
    p0 = np.linspace(3.2, 4.05, 200)
    axes[2].plot(p0, 3*p6*(p6-p0)/8, color=GREEN)
    axes[2].axhline(0,color=INK,lw=1); axes[2].axvline(p6,color=ORANGE,ls="--",lw=1)
    axes[2].set(xlabel=r"Target perimeter, $p_0$", ylabel="Affine energy curvature",
                title="C  Stability of this path")
    save(fig, "linear_response", "affine_hexagon")


def mode_coupling():
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), layout="constrained")
    tau = np.array([.2, 5.]); a = np.array([.6, .4]); lamb = 1/tau; weights = a*lamb
    x = np.arange(2)
    axes[0].bar(x-.16, lamb, width=.3, color=BLUE, label=r"$\lambda_\alpha$")
    axes[0].bar(x+.16, weights, width=.3, color=ORANGE, label=r"$w_\alpha$")
    axes[0].set(yscale="log", xticks=x, xticklabels=["Mode 1","Mode 2"], title="A  Stiffness and coupling")
    axes[0].legend()
    t = np.linspace(0, 15, 400)
    axes[1].plot(t,.2+sum(v*np.exp(-t/r) for v,r in zip(a,tau)),color=BLUE)
    axes[1].axhline(.2,color=INK,lw=1,ls=":")
    axes[1].set(xlabel="Time",ylabel="G(t)",title="B  Two relaxation times")
    w = np.logspace(-2, 2, 300); complex_mod=.2+sum(v*1j*w*r/(1+1j*w*r) for v,r in zip(a,tau))
    axes[2].semilogx(w,complex_mod.real,color=BLUE,label=r"$G'$")
    axes[2].semilogx(w,complex_mod.imag,color=ORANGE,label=r"$G''$")
    axes[2].set(xlabel="Angular frequency",ylabel="Modulus",title="C  Collective response")
    axes[2].legend()
    save(fig, "linear_response", "mode_coupling")


def polygon_bounds():
    fig, axes = plt.subplots(1, 4, figsize=(11, 3), layout="constrained")
    for ax, n, color in zip(axes,[4,5,6,7],[BLUE,ORANGE,GREEN,INK]):
        radius=np.sqrt(2/(n*np.sin(2*np.pi/n)))
        theta=np.arange(n)*2*np.pi/n+np.pi/2
        ax.add_patch(Polygon(radius*np.c_[np.cos(theta),np.sin(theta)],fc=color+"22",ec=color,lw=2))
        ax.text(0,0,"A = 1",ha="center",va="center",fontsize=15)
        ax.text(0,-1.25,rf"$p_{{{n}}}={np.sqrt(4*n*np.tan(np.pi/n)):.3f}$",ha="center",fontsize=18)
        ax.set(title=f"{n} sides",xlim=(-1.1,1.1),ylim=(-1.55,1.1),aspect="equal")
        ax.axis("off")
    save(fig, "rigidity", "polygon_bounds")


def taut_springs():
    fig, axes=plt.subplots(1,2,figsize=(10,3.8),layout="constrained")
    y=.5
    axes[0].scatter([-1,1],[0,0],color=INK,s=55)
    axes[0].scatter([0],[y],color=ORANGE,s=60,zorder=3)
    for end in [(-1,0),(1,0)]:
        axes[0].plot([end[0],0],[0,y],color=BLUE)
        arrow(axes[0],(0,y),(end[0]*.45,y*.55),BLUE)
    axes[0].plot([-1,1],[0,0],color=INK,lw=1,ls=":")
    axes[0].text(.1,.7,"Transverse displacement",fontsize=12)
    axes[0].set(xlim=(-1.3,1.3),ylim=(-.4,1.1),aspect="equal",title="A  Prestress supplies a restoring force")
    axes[0].axis("off")
    y=np.linspace(-1,1,301)
    for l0,color,label in [(.8,BLUE,"Tension"),(1,ORANGE,"Unstressed"),(1.2,GREEN,"Compression")]:
        energy=(np.sqrt(1+y*y)-l0)**2-(1-l0)**2
        axes[1].plot(y,energy,color=color,label=label)
    axes[1].axhline(0,color=INK,lw=1,ls=":")
    axes[1].set(xlabel="Transverse coordinate, y",ylabel="E(y) - E(0)",title="B  Same geometry, different stability")
    axes[1].legend()
    save(fig,"rigidity","taut_springs")


def event_strain():
    fig,axes=plt.subplots(1,3,figsize=(12,3.6),layout="constrained")
    t=np.linspace(0,12,1201);times=np.arange(1,11);increment=.055
    for signs,color,label in [(np.ones(10),BLUE,"Aligned events"),((-1.)**np.arange(10),ORANGE,"Alternating events")]:
        count=np.sum(t[:,None]>=times[None,:],axis=1)
        plastic=np.sum((t[:,None]>=times[None,:])*signs[None,:]*increment,axis=1)
        axes[0].plot(t,count,color=color,label=label)
        axes[1].plot(t,plastic,color=color)
        axes[2].plot(t,.07*t-plastic,color=color)
    for ax in axes:ax.set_xlabel("Time")
    axes[0].set(ylabel="Cumulative event count",title="A  Same count");axes[0].legend(fontsize=10)
    axes[1].set(ylabel=r"$\gamma_{\rm pl}$",title="B  Different signed strain")
    axes[2].set(ylabel=r"$\sigma=\gamma-\gamma_{\rm pl}$",title="C  Different stress release")
    save(fig,"rearrangements","event_strain")


def barrier_statistics():
    fig,axes=plt.subplots(1,3,figsize=(12,3.8),layout="constrained")
    e=np.linspace(0,5,300)
    axes[0].plot(e,4*e*np.exp(-2*e),color=BLUE)
    axes[0].set(xlabel=r"$E/\bar E$",ylabel="Normalized probability density",title="A  A chosen Gamma distribution")
    ratio=np.linspace(0,20,250)
    axes[1].semilogy(ratio,np.exp(-ratio),color=ORANGE,label="Use the mean barrier")
    axes[1].semilogy(ratio,(1+ratio/2)**-2,color=BLUE,label="Average the rates")
    axes[1].set(xlabel=r"$\bar E/X$",ylabel=r"$\langle k\rangle/k_0$",title="B  The low-barrier tail matters")
    axes[1].legend(fontsize=10)
    h=np.linspace(-2,2,200)
    axes[2].plot(h,np.exp(h),color=BLUE,label="Forward")
    axes[2].plot(h,np.exp(-h),color=ORANGE,label="Backward")
    axes[2].set(xlabel=r"$\sigma\Omega^{\ddagger}/X$",ylabel="Rate / zero-stress rate",title="C  Stress biases direction")
    axes[2].legend()
    save(fig,"rearrangements","barrier_statistics")


def protocols():
    fig,axes=plt.subplots(2,3,figsize=(12,5.8),layout="constrained")
    t=np.linspace(0,10,600);w=1.5
    strain=[np.ones_like(t),np.cos(w*t),.2*t]
    stress=[np.exp(-t),w*w/(1+w*w)*np.cos(w*t)-w/(1+w*w)*np.sin(w*t),.2*(1-np.exp(-t))]
    for j,name in enumerate(["Step strain","Oscillatory strain","Constant strain rate"]):
        axes[0,j].plot(t,strain[j],color=BLUE);axes[0,j].set_title(name)
        axes[1,j].plot(t,stress[j],color=ORANGE);axes[1,j].set_xlabel("Time")
    axes[0,0].plot([0,0],[0,1],color=BLUE)
    axes[0,0].set_ylabel("Input strain");axes[1,0].set_ylabel("Stress response")
    save(fig,"rheology_activity","protocols")


def active_msd():
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout="constrained")
    t=np.logspace(-2,3,400)
    for dr,color in [(.2,BLUE),(1,ORANGE),(5,GREEN)]:
        msd=2/dr**2*(dr*t+np.expm1(-dr*t))
        axes[0].loglog(t,msd,color=color,label=rf"$D_r={dr}$")
    axes[0].set(xlabel="Time",ylabel=r"$\langle|\Delta\mathbf{R}|^2\rangle$",title="A  A free persistent particle")
    axes[0].legend()
    lam=np.logspace(-2,2,300)
    for tau,color in [(.1,BLUE),(1,ORANGE),(10,GREEN)]:
        axes[1].loglog(lam,1/(lam*(lam+1/tau)),color=color,label=rf"$\tau_p={tau}$")
    axes[1].set(xlabel=r"Mode stiffness, $\lambda$",ylabel=r"$\langle u^2\rangle$",title="B  Fixed force variance")
    axes[1].legend()
    save(fig,"rheology_activity","active_msd")


def shape_tensor():
    fig,axes=plt.subplots(1,3,figsize=(12,3.8),layout="constrained")
    angles=[np.arange(12)*np.pi/12,np.zeros(12)]
    order=[]
    for ax,theta,title in zip(axes[:2],angles,["A  Distributed directions","B  Aligned directions"]):
        for j,angle in enumerate(theta):
            center=(j%4,(j//4)*.9)
            ax.add_patch(Ellipse(center,.85,.35,angle=angle*180/np.pi,facecolor="#d8edf6",edgecolor=BLUE,lw=1.6))
        q=(.85**2-.35**2)/(2*(.85**2+.35**2))
        order.append(q*np.sqrt(np.mean(np.cos(2*theta))**2+np.mean(np.sin(2*theta))**2))
        ax.set(xlim=(-.6,3.6),ylim=(-.6,2.5),aspect="equal",title=title);ax.axis("off")
    axes[2].bar([0,1],order,color=[ORANGE,BLUE],width=.5)
    axes[2].set(xticks=[0,1],xticklabels=["Distributed","Aligned"],ylabel="Mean tensor amplitude",title="C  Same individual elongation")
    save(fig,"geometry_3d","shape_tensor")


def curvature_thickness():
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout="constrained")
    for radius,color,label in [(.8,BLUE,"Inner"),(1,INK,"Mid-surface"),(1.2,ORANGE,"Outer")]:
        axes[0].add_patch(Circle((0,0),radius,fill=False,ec=color,lw=2,ls=":" if radius==1 else "-"))
    axes[0].plot([.8,1.2],[0,0],color=GREEN,lw=4)
    axes[0].text(1,.17,"h",ha="center",color=GREEN,fontsize=17)
    axes[0].text(0,-1.5,"Finite thickness changes surface areas",ha="center",fontsize=12)
    axes[0].set(xlim=(-1.7,1.7),ylim=(-1.7,1.5),aspect="equal",title="A  A spherical layer");axes[0].axis("off")
    z=np.linspace(0,.8,200)
    axes[1].plot(z,((1+z/2)/(1-z/2))**2,color=BLUE,label="Exact geometry")
    axes[1].plot(z,1+2*z,"--",color=ORANGE,label="Thin-shell expansion")
    axes[1].set(xlabel="Thickness / radius, h/R",ylabel=r"$A_{\rm out}/A_{\rm in}$",title="B  Inner and outer area")
    axes[1].legend()
    save(fig,"geometry_3d","curvature_thickness")


def identifiability():
    fig,axes=plt.subplots(1,3,figsize=(12,3.8),layout="constrained")
    x=np.linspace(-1,1.5,250)
    for scale,color in [(1,BLUE),(3,ORANGE),(8,GREEN)]:
        energy=scale*((x-1)**2+.5*(x+1)**2)
        axes[0].plot(x,energy-energy.min(),color=color,label=f"Scale {scale}")
        force=np.linspace(0,1,100)
        axes[2].plot(force,force/(3*scale),color=color)
    axes[0].set(xlabel="Coordinate, x",ylabel="Energy above minimum",title="A  Identical equilibrium shape");axes[0].legend()
    a=np.linspace(.2,5,150);b=np.linspace(.2,5,150);aa,bb=np.meshgrid(a,b)
    cs=axes[1].contour(aa,bb,(aa-bb)/(aa+bb),levels=[-.5,0,.5],colors=[BLUE,INK,ORANGE])
    axes[1].clabel(cs,inline=True,fontsize=12)
    axes[1].set(xlabel="Parameter a",ylabel="Parameter b",title="B  Contours of the same shape")
    axes[2].set(xlabel="Calibrated external force",ylabel="Additional displacement",title="C  A measurement breaks degeneracy")
    save(fig,"morphogenesis","identifiability")


def plastic_lock():
    fig,axes=plt.subplots(1,3,figsize=(12,3.8),layout="constrained")
    t=np.linspace(0,55,11001);dt=t[1]-t[0]
    drive=np.where(t<35,.6+.5*np.sin(2*np.pi*t/7),0)
    axes[0].plot(t,drive,color=INK)
    axes[0].set(xlabel="Time",ylabel="External force",title="A  Cyclic loading, then unloading")
    for beta,color,label in [(0,ORANGE,"Fixed reference"),(.25,BLUE,"Plastic reference")]:
        length=np.ones_like(t);rest=np.ones_like(t)
        for j in range(len(t)-1):
            velocity=drive[j]-(length[j]-rest[j])
            reference_rate=beta*max(velocity,0) if length[j]>rest[j] else 0
            length[j+1]=length[j]+dt*velocity
            rest[j+1]=rest[j]+dt*reference_rate
        axes[1].plot(t,length,color=color,label=label)
        axes[2].plot(t,rest,color=color)
    axes[1].set(xlabel="Time",ylabel="Length, L",title="B  Residual elongation");axes[1].legend(fontsize=10)
    axes[2].set(xlabel="Time",ylabel=r"Reference length, $\lambda$",title="C  Memory in a reference state")
    for ax in axes:ax.axvline(35,color=GREEN,ls=":",lw=1.2)
    save(fig,"morphogenesis","plastic_lock")


def contact_adhesion():
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.3), layout="constrained")
    ax = axes[0]
    half_angle = np.deg2rad(50)
    free = np.array([np.sin(half_angle), np.cos(half_angle)])
    contact = np.array([0, -2*np.cos(half_angle)])
    for direction in (free, free*np.array([-1, 1])):
        ax.plot([0, 1.5*direction[0]], [0, 1.5*direction[1]], color="#cbd8e0", lw=4)
        arrow(ax, (0, 0), direction, BLUE)
    ax.plot([0, 0], [0, -1.6], color="#cbd8e0", lw=4)
    arrow(ax, (0, 0), contact, ORANGE)
    ax.add_patch(Arc((0, 0), 1.1, 1.1, theta1=40, theta2=140, color=INK, lw=1.5))
    ax.text(0, .70, r"$\theta=100^\circ$", ha="center")
    ax.text(1.15, .85, r"$\gamma_{\rm cf}=1$", ha="left")
    ax.text(.13, -1.10, r"$\gamma_{\rm cc}=2\cos(\theta/2)$", ha="left")
    ax.text(0, 1.2, r"$\alpha=\gamma_{\rm cc}/(2\gamma_{\rm cf})=\cos(\theta/2)$", ha="center")
    ax.scatter([0], [0], color=INK, s=40, zorder=5)
    ax.set(xlim=(-1.6, 2.3), ylim=(-1.7, 1.5), aspect="equal", title="A  Local force balance")
    ax.axis("off")
    for patch in ax.patches:
        if isinstance(patch, FancyArrowPatch):
            patch.set_zorder(3)
    length = np.linspace(0, 5, 100)
    for work, color in [(0, INK), (.8, BLUE), (1.6, ORANGE)]:
        axes[1].plot(length, 10-work*length, color=color, label=f"W = {work:g}")
    axes[1].set(xlabel=r"Contact length, $L_{\rm c}$", ylabel=r"Surface energy, $E_{\rm surf}$",
                title="B  Fixed total perimeter = 10")
    axes[1].legend()
    axes[1].text(.04, .1, r"$E_{\rm surf}=\gamma_{\rm cf}\sum_iP_i-WL_{\rm c}$", transform=axes[1].transAxes)
    save(fig, "foundations", "contact_adhesion")


def contact_rigidity():
    points = np.array([[-1.2, .2], [-.5, 1.1], [0, 0], [.8, -.9], [1.35, .35]])
    base = [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)]
    fig, axes = plt.subplots(1, 3, figsize=(12.3, 4), layout="constrained")
    spectra = []
    for ax, edges, title in [(axes[0], base, "A  Connected, one hinge"),
                              (axes[1], base+[(1, 4)], "B  One added independent brace")]:
        jac = np.zeros((len(edges), 2*len(points)))
        for e, (i, j) in enumerate(edges):
            direction = (points[i]-points[j])/np.linalg.norm(points[i]-points[j])
            jac[e, 2*i:2*i+2] = direction
            jac[e, 2*j:2*j+2] = -direction
            ax.plot(*points[[i, j]].T, color=ORANGE if (i, j)==(1, 4) else BLUE, lw=2.5)
        ax.scatter(*points.T, s=60, color=INK, zorder=4)
        vals = np.linalg.svd(jac, compute_uv=False)
        spectra.append(np.sort(np.r_[vals, np.zeros(10-len(vals))]))
        rank = np.linalg.matrix_rank(jac)
        ax.text(.5, -.14, f"rank J = {rank}; zero modes = {10-rank}", ha="center", transform=ax.transAxes)
        ax.set(xlim=(-1.4, 1.6), ylim=(-1.2, 1.4), aspect="equal", title=title)
        ax.axis("off")
    angle = np.deg2rad(24)
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    moved = points.copy()
    moved[3:] = points[3:]@rotation.T
    for i, j in [(2, 3), (3, 4), (4, 2)]:
        axes[0].plot(*moved[[i, j]].T, color=ORANGE, ls="--", lw=1.8)
    axes[0].text(.5, -.04, "Dashed: finite hinge rotation", ha="center", transform=axes[0].transAxes, fontsize=10)
    for vals, offset, color, label in [(spectra[0], -.15, BLUE, "Hinge"), (spectra[1], .15, ORANGE, "Braced")]:
        axes[2].bar(np.arange(1, 11)+offset, vals, width=.28, color=color, label=label)
    axes[2].set(xlabel="Coordinate-mode index", ylabel="Singular value of J", title="C  Independent constraints")
    axes[2].set_xticks([1, 4, 7, 10])
    axes[2].legend()
    save(fig, "rigidity", "contact_rigidity")


def activity_bias():
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), layout="constrained")
    biases = [-.6, 0, .6]
    activity = 4
    x = np.arange(3)
    axes[0].bar(x-.17, activity*(1+np.array(biases))/2, width=.34, color=BLUE, label="Positive events")
    axes[0].bar(x+.17, activity*(1-np.array(biases))/2, width=.34, color=ORANGE, label="Negative events")
    axes[0].set_xticks(x, [str(b) for b in biases])
    axes[0].set(xlabel="Directional bias, b", ylabel="Event rate", title=r"A  Same activity, $\Gamma=4$")
    axes[0].set_ylim(0, 4.8)
    axes[0].legend(fontsize=10)
    gamma = np.linspace(0, 8, 100)
    for bias, color in [(0, INK), (.25, BLUE), (.75, ORANGE)]:
        axes[1].plot(gamma, gamma*bias, color=color, label=f"b = {bias:g}")
    axes[1].set(xlabel=r"Total event activity, $\Gamma$", ylabel="Net plastic strain rate", title="B  Activity needs a directional bias")
    axes[1].legend()
    bias = np.linspace(-1, 1, 150)
    axes[2].plot(bias, activity*bias, color=GREEN)
    axes[2].axhline(0, color=INK, lw=1)
    axes[2].axvline(0, color=INK, lw=1, ls=":")
    axes[2].set(xlabel="Directional bias, b", ylabel="Net plastic strain rate", title=r"C  $\dot\varepsilon^{pl}=\varepsilon_e\Gamma b$")
    save(fig, "rheology_activity", "activity_bias")


def porosity_feedback():
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.3), layout="constrained")
    ax = axes[0]
    labels = [("Nodal–Lefty\nsignal", (.15, .82)), ("Contact\ngeometry", (.76, .82)),
              ("Fluid\nchannels", (.76, .25)), ("Transport\n+ reaction", (.15, .25))]
    for label, pos in labels:
        ax.text(*pos, label, ha="center", va="center", fontsize=11,
                bbox={"boxstyle":"round,pad=.5", "fc":"#eaf2f5", "ec":BLUE})
    for a, b in [((.36, .82), (.58, .82)), ((.76, .66), (.76, .41)),
                 ((.59, .25), (.37, .25)), ((.15, .41), (.15, .66))]:
        arrow(ax, a, b)
    ax.set(xlim=(-.15, 1.08), ylim=(0, 1), title="A  Geometry mediates the feedback")
    ax.axis("off")
    t = np.linspace(0, 1, 120)
    eps = .2-.1*t
    axes[1].plot(t, eps/.2, color=BLUE, label=r"Porosity, $\epsilon/\epsilon_0$")
    axes[1].plot(t, .2/eps, color=ORANGE, label=r"Concentration, $c/c_0$")
    axes[1].set(xlabel="Prescribed closure progress", ylabel="Ratio to initial value", title=r"B  Retained amount: $c=n/\epsilon$")
    axes[1].legend(fontsize=10)
    dp = (eps/.2)**2
    axes[2].plot(t, dp, color=BLUE, label=r"$D_{\rm p}/D_0=(\epsilon/\epsilon_0)^2$")
    axes[2].plot(t, np.sqrt(dp), color=GREEN, label=r"$\ell/\ell_0=\sqrt{D_{\rm p}/D_0}$")
    axes[2].set(xlabel="Prescribed closure progress", ylabel="Ratio to initial value", title="C  Illustrative transport law")
    axes[2].legend(fontsize=10)
    save(fig, "geometry_3d", "porosity_feedback")


if __name__ == "__main__":
    for draw in (model, response, incompatibility, t1, rheology, dimensions, inference,
                 force_balance, topology, affine_hexagon, mode_coupling, polygon_bounds,
                 taut_springs, event_strain, barrier_statistics, protocols, active_msd,
                 shape_tensor, curvature_thickness, identifiability, plastic_lock,
                 contact_adhesion, contact_rigidity, activity_bias, porosity_feedback):
        draw()
    print("Wrote 25 SVG teaching figures and 25 PNG exports.")
