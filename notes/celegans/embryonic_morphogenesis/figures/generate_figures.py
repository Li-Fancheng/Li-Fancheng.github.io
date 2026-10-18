"""Original explanatory figures; cartoons and model curves are not embryo data.

Run with Python, NumPy and Matplotlib. SVG glyphs are outlined for portable
display; PNG companions support local visual inspection and static export.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Polygon, FancyBboxPatch, Rectangle
from matplotlib import font_manager
import numpy as np

OUT = Path(__file__).resolve().parent
for name in ("Arial.ttf", "Arial Bold.ttf"):
    p = Path("/System/Library/Fonts/Supplemental") / name
    if p.exists():
        font_manager.fontManager.addfont(p)
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 15, "axes.titlesize": 17, "axes.labelsize": 14,
    "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12,
    "mathtext.fontset": "cm", "svg.fonttype": "path",
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "text.color": "#233849", "axes.labelcolor": "#233849",
})
INK, MUTED = "#233849", "#657887"
GREEN, TEAL, PINK, GOLD, PURPLE = "#2c9176", "#288fa1", "#ce6791", "#d4a542", "#8164ac"


def save(fig, name):
    fig.savefig(OUT / (name + ".svg"), bbox_inches="tight", pad_inches=.2)
    fig.savefig(OUT / (name + ".png"), dpi=180, bbox_inches="tight", pad_inches=.2)
    plt.close(fig)


def diagram(ax, title, xlim=(0, 10), ylim=(0, 7)):
    ax.set(xlim=xlim, ylim=ylim, aspect="equal")
    ax.axis("off")
    ax.set_title(title, loc="left", fontweight="bold", pad=14)


def arrow(ax, a, b, color=INK, lw=2):
    ax.annotate("", xy=b, xytext=a,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw))


def cell(ax, x, y, w, h, color, label=""):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.16",
                      facecolor=color, edgecolor="white", lw=2)
    ax.add_patch(p)
    if label:
        ax.text(x+w/2, y+h/2, label, ha="center", va="center", fontsize=13, color="white")


def overview():
    fig, ax = plt.subplots(figsize=(12, 6.5), layout="constrained")
    ax.set(xlim=(0, 10), ylim=(-.5, 6.5)); ax.axis("off")
    labels = ["Cleavage + polarity", "Gastrulation", "Dorsal intercalation",
              "Ventral enclosure", "Elongation", "Twitching + hatching"]
    starts, ends = [0, 1.5, 4.3, 4.8, 6.1, 7.1], [3.0, 4.5, 6.4, 6.7, 9.1, 9.8]
    colors = [MUTED, GREEN, TEAL, PINK, GOLD, PURPLE]
    for i, (label, x0, x1, c) in enumerate(zip(labels, starts, ends, colors)):
        y = 5-i
        ax.add_patch(FancyBboxPatch((x0, y), x1-x0, .56,
                     boxstyle="round,pad=0.02,rounding_size=0.12",
                     facecolor=c, alpha=.9))
        ax.text(x0+.12, y+.28, label, fontsize=14, va="center", color="white")
    arrow(ax, (0, -.22), (9.9, -.22))
    ax.text(5, -.6, "Developmental order only: bar lengths are not durations",
            ha="center", fontsize=13)
    ax.set_title("Overlapping processes build internal organization and external shape",
                 loc="left", fontweight="bold", pad=18)
    save(fig, "01_developmental_map")


def gastrulation():
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.9), layout="constrained")
    for k, ax in enumerate(axes):
        diagram(ax, ["A  Surface exposure", "B  Coupled constriction", "C  Internalization"][k])
        ax.plot([.3, 9.7], [5.1, 5.1], ls="--", color=MUTED, lw=1.2)
        ax.text(5, 6.25, "Exterior", ha="center", fontsize=13)
        ax.text(5, .4, "Embryo interior", ha="center", fontsize=13)
        w = [1.75, 1.0, 0.0][k]
        top = [5.1, 4.5, 2.5][k]
        for center, label in [(3.9, "Ea"), (6.1, "Ep")]:
            if k < 2:
                points = [(center-w/2, top), (center+w/2, top),
                          (center+1.05, 1.8), (center-1.05, 1.8)]
                ax.add_patch(Polygon(points, facecolor=GREEN, edgecolor="white", lw=2))
                ax.plot([center-w/2, center+w/2], [top, top], color=GOLD, lw=6)
            else:
                ax.add_patch(Ellipse((center, 2.5), 2.0, 2.2, facecolor=GREEN, edgecolor="white", lw=2))
            ax.text(center, 2.9 if k < 2 else 2.5, label, color="white", ha="center", fontsize=14)
        if k == 0:
            cell(ax, .5, 2.1, 2.4, 3.0, TEAL, "Neighbor")
            cell(ax, 7.1, 2.1, 2.4, 3.0, TEAL, "Neighbor")
        elif k == 1:
            cell(ax, .5, 2.5, 2.9, 2.6, TEAL)
            cell(ax, 6.6, 2.5, 2.9, 2.6, TEAL)
            arrow(ax, (2, 5.7), (3.5, 5.7), GOLD)
            arrow(ax, (8, 5.7), (6.5, 5.7), GOLD)
            arrow(ax, (5, 3.7), (5, 2.1), INK)
        else:
            cell(ax, .5, 4.2, 4.5, .9, TEAL)
            cell(ax, 5, 4.2, 4.5, .9, TEAL)
            ax.text(5, 5.65, "New neighbor contact", ha="center", fontsize=13)
    fig.supxlabel("Conceptual section through Ea/Ep; cell shapes and sizes are schematic", fontsize=13)
    save(fig, "02_gastrulation")


def intercalation():
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 5.3), layout="constrained")
    for ax in axes[:2]:
        diagram(ax, "A  Two rows" if ax is axes[0] else "B  Interdigitated row",
                xlim=(0, 6), ylim=(0, 8))
        ax.plot([3, 3], [.5, 7.5], ls=":", color=MUTED)
        ax.text(.15, 7.5, "A", fontsize=13); ax.text(.15, .25, "P", fontsize=13)
    for i in range(3):
        cell(axes[0], .7, 1.1+i*1.9, 2.1, 1.7, TEAL, "L"+str(3-i))
        cell(axes[0], 3.2, 1.1+i*1.9, 2.1, 1.7, "#75b9c4", "R"+str(3-i))
        arrow(axes[0], (2.2, 2.0+i*1.9), (3.3, 2.0+i*1.9), INK)
        arrow(axes[0], (3.8, 2.0+i*1.9), (2.8, 2.0+i*1.9), INK)
    for i, name in enumerate(["R3", "L3", "R2", "L2", "R1", "L1"]):
        cell(axes[1], 1, .8+i*1.05, 4, .95, TEAL if name[0]=="L" else "#75b9c4", name)
    ax=axes[2]; diagram(ax, "C  Track the contact graph", xlim=(0, 6), ylim=(0, 8))
    pos={"i":(1.1,5.9),"j":(4.9,5.9),"k":(1.1,2.3),"l":(4.9,2.3)}
    for a,b in [("i","j"),("j","l"),("l","k"),("k","i")]:
        ax.plot([pos[a][0],pos[b][0]],[pos[a][1],pos[b][1]],color="#c8d3db",lw=3)
    ax.plot([1.1,4.9],[5.9,2.3],ls="--",color=PINK,lw=3)
    ax.plot([4.9,1.1],[5.9,2.3],color=GREEN,lw=3)
    for label,(x,y) in pos.items():
        ax.add_patch(Ellipse((x,y),1,1,facecolor=INK)); ax.text(x,y,label,ha="center",va="center",color="white")
    ax.text(3,7.3,"Abstract local graph",ha="center",fontsize=13)
    ax.text(3,1.0,"Dashed: lost contact\nGreen: gained contact",ha="center",fontsize=13)
    fig.supxlabel("Dorsal view: A–P is vertical. The graph is illustrative, not a literal T1 sequence.", fontsize=13)
    save(fig, "03_intercalation")


def enclosure():
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 5.3), layout="constrained")
    for k, ax in enumerate(axes):
        diagram(ax, ["A  Leading-cell migration", "B  New contacts + pocket", "C  Continuous epidermis"][k],
                xlim=(0, 6), ylim=(0, 8))
        ax.add_patch(Ellipse((3,4.0),5.5,7.1,facecolor="#edf1f4",edgecolor=MUTED,lw=1.5))
        ax.plot([3,3],[.5,7.5],ls=":",color=MUTED)
        gap = [1.6,.75,0][k]
        for i in range(5):
            y=1.3+i*1.05
            g = gap if i < 3 else (gap if k==0 else 0)
            cell(ax,.8,y,2.2-g/2,.92,PINK)
            cell(ax,3+g/2,y,2.2-g/2,.92,PINK)
        if k == 0:
            arrow(ax,(1.8,5.9),(2.95,5.9),INK); arrow(ax,(4.2,5.9),(3.05,5.9),INK)
            ax.text(3,6.9,"Leading pairs",ha="center",fontsize=13)
        if k == 1:
            ax.add_patch(Ellipse((3,3),.85,3.2,fill=False,edgecolor=GOLD,lw=4))
            ax.text(3,6.9,"First stable contacts",ha="center",fontsize=13)
            ax.text(3,.15,"Posterior pocket remains",ha="center",fontsize=13)
        if k == 2:
            ax.plot([3,3],[1.3,6.4],color=INK,lw=2)
            ax.text(3,6.9,"Ventral junctions",ha="center",fontsize=13)
        ax.text(.05,7.5,"A",fontsize=12); ax.text(.05,.25,"P",fontsize=12)
    fig.supxlabel("Ventral surface schematic: underlying neuroblasts are the gray substrate", fontsize=13)
    save(fig, "04_enclosure")


def elongation():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6), layout="constrained")
    ax=axes[0]; diagram(ax,"A  Different physical roles",xlim=(0,9),ylim=(0,8))
    for y,h,c,label in [(5.1,1.8,TEAL,"Dorsal epidermis"),(3.2,1.7,GOLD,"Seam cells"),
                         (1.2,1.8,PINK,"Ventral epidermis")]:
        cell(ax,.7,y,7.4,h,c)
        ax.text(8.8,y+h/2,label,rotation=90,va="center",ha="center",fontsize=12)
    for y0,y1 in [(1.3,2.9),(5.2,6.8)]:
        for x in np.linspace(1.2,7.7,15):
            ax.plot([x,x],[y0,y1],color=INK,lw=1.3)
    for x in [2,4.5,7]:
        arrow(ax,(x,4.8),(x,4.2),INK); arrow(ax,(x,3.3),(x,3.9),INK)
    arrow(ax,(1,7.5),(7.8,7.5),INK)
    ax.text(4.4,7.75,"Anterior–posterior axis",ha="center",fontsize=13)
    ax.text(4.4,.35,"Unrolled epidermal surface (schematic)",ha="center",fontsize=12)
    ax=axes[1]; diagram(ax,"B  Constant-volume geometry",xlim=(0,9),ylim=(0,8))
    ax.add_patch(Ellipse((3.1,5.9),4.5,2.5,facecolor="#dbe9ed",edgecolor=INK,lw=2))
    ax.add_patch(Ellipse((4.5,2.3),7.5,1.55,facecolor="#dbe9ed",edgecolor=INK,lw=2))
    arrow(ax,(4.5,4.5),(4.5,3.5),GREEN)
    ax.text(6.7,5.8,r"$L_0,\ R_0$",ha="center")
    ax.text(4.5,.6,r"$R/R_0=\lambda^{-1/2}$  if  $V=\pi R^2L$",ha="center")
    ax.text(4.5,7.7,"Shorter circumference; longer body",ha="center",fontsize=13)
    save(fig,"05_elongation_mechanics")


def ratchet():
    fig, axes=plt.subplots(1,2,figsize=(12,4.8),layout="constrained")
    t=np.linspace(0,4,800)
    elastic=.06*np.sin(2*np.pi*t)
    plastic=.035*t
    axes[0].plot(t,elastic,color=MUTED,lw=2,label="Elastic response")
    axes[0].plot(t,elastic+plastic,color=PURPLE,lw=2.5,label="With evolving reference shape")
    axes[0].plot(t,plastic,ls="--",color=GREEN,label="Retained component")
    axes[0].set(title="A  Cycles with or without accumulation",
                xlabel="Time / contraction period",ylabel="Axial strain (model)")
    axes[0].legend(loc="upper left",fontsize=11)
    e=np.linspace(0,1,200)
    axes[1].plot(e,e,color=MUTED,lw=2,label="Load")
    axes[1].plot([1,.35],[1,0],color=PURPLE,lw=2.5,label="Unload after remodeling")
    axes[1].plot([.35,.35],[0,-.1],color=GREEN,ls="--")
    axes[1].annotate("Residual strain",xy=(.35,0),xytext=(.45,.22),
                    arrowprops=dict(arrowstyle="->",color=GREEN),fontsize=13)
    axes[1].set(title="B  A changed reference configuration",
                xlabel="Strain (normalized)",ylabel="Stress (normalized)",
                xlim=(-.05,1.1),ylim=(-.13,1.15))
    axes[1].legend(loc="upper left")
    for ax in axes: ax.grid(alpha=.15)
    fig.supxlabel("Illustrative model curves; no experimental values or fit",fontsize=12)
    save(fig,"06_viscoplastic_lock")


def measurements():
    fig,axes=plt.subplots(1,3,figsize=(13.5,4.8),layout="constrained")
    ax=axes[0]
    theta=np.linspace(-.9*np.pi,.9*np.pi,250); x=np.cos(theta); y=np.sin(theta)
    ax.plot(x,y,color=TEAL,lw=8,solid_capstyle="round")
    ax.plot(x,y,color="white",lw=1.8,ls="--")
    ax.plot([x[0],x[-1]],[y[0],y[-1]],color=PINK,lw=2)
    ax.set(aspect="equal",title="A  Body length vs. endpoints")
    ax.text(0,0,"Centerline: L\nEndpoint distance: D",ha="center",va="center",fontsize=13)
    ax.axis("off")
    u=np.linspace(0,1,200)
    axes[1].plot(u,1-u,color=TEAL,lw=2.5,label="Constant area-loss rate")
    axes[1].plot(u,(1-u)**2,color=PINK,lw=2.5,label="Constant radial speed")
    axes[1].set(title="B  Closure laws can differ",xlabel="Normalized time",
                ylabel=r"Open area $A/A_0$",xlim=(0,1),ylim=(0,1.05))
    axes[1].legend(fontsize=10)
    lam=np.linspace(1,4,200)
    axes[2].plot(lam,lam**(-.5),color=TEAL,lw=2.5,label=r"$R/R_0=\lambda^{-1/2}$")
    axes[2].plot(lam,lam**.5,color=GOLD,lw=2.5,label=r"$S/S_0\simeq\lambda^{1/2}$")
    axes[2].set(title="C  Fixed volume is not fixed area",xlabel=r"Length ratio $\lambda=L/L_0$",
                ylabel="Geometric ratio")
    axes[2].legend(fontsize=11)
    for ax in axes[1:]: ax.grid(alpha=.15)
    fig.supxlabel("A: schematic centerline. B–C: idealized geometric predictions, not embryo data.",fontsize=12)
    save(fig,"07_quantitative_geometry")


if __name__=="__main__":
    for make in (overview,gastrulation,intercalation,enclosure,elongation,ratchet,measurements):
        make()
