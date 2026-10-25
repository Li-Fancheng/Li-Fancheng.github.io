"""Plot only the saved numerical results; never replace a failed solve with a curve."""

from pathlib import Path
import json
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib import font_manager

HERE = Path(__file__).resolve().parent
ROOT = (
    HERE.parent
    if (HERE.parent / "rigidity" / "body.md").exists()
    else HERE / "generated"
)
BLUE, ORANGE, GREEN, INK = "#277da8", "#d67837", "#429582", "#263e50"
for filename in ["Arial.ttf", "Arial Bold.ttf"]:
    font = Path("/System/Library/Fonts/Supplemental") / filename
    if font.exists():
        font_manager.fontManager.addfont(font)
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 13,
        "axes.titlesize": 14,
        "axes.labelsize": 13,
        "legend.fontsize": 10.5,
        "mathtext.fontset": "cm",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "lines.linewidth": 2.2,
        "svg.fonttype": "path",
        "svg.hashsalt": "vertex-unified-network-20261025",
        "savefig.facecolor": "white",
    }
)


def save(fig, slug, name):
    dest = ROOT / slug / "figures"
    dest.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        dest / (name + ".svg"),
        bbox_inches="tight",
        pad_inches=0.2,
        metadata={"Date": None},
    )
    fig.savefig(dest / (name + ".png"), bbox_inches="tight", pad_inches=0.2, dpi=300)
    plt.close(fig)


def polygon_points(cell, x, L=4.0):
    p = x[cell]
    jumps = -np.rint((np.roll(p, -1, axis=0) - p) / L)
    images = np.vstack([np.zeros(2), np.cumsum(jumps[:-1], axis=0)])
    p = p + L * images
    return p - L * np.floor(p.mean(axis=0) / L)


def draw_mesh(ax, data, x, title, labels=False):
    for ci, cell in enumerate(data["preparation"]["cells"]):
        p = polygon_points(cell, x)
        fill = {5: "#f3e3d6", 6: "#e6edf2", 7: "#dcece5"}.get(len(cell), "#ebe6ed")
        for dx in [-4, 0, 4]:
            for dy in [-4, 0, 4]:
                q = p + [dx, dy]
                if (
                    q[:, 0].max() < 0
                    or q[:, 0].min() > 4
                    or q[:, 1].max() < 0
                    or q[:, 1].min() > 4
                ):
                    continue
                ax.add_patch(Polygon(q, facecolor=fill, edgecolor=INK, lw=1.1))
        if labels:
            center = p.mean(axis=0) % 4
            ax.text(*center, str(ci), ha="center", va="center", fontsize=10)
    u, v = 2, 3
    d = x[v] - x[u]
    d -= np.rint(d / 4) * 4
    line = np.array([x[u], x[u] + d])
    line -= 4 * np.floor(line.mean(axis=0) / 4)
    for dx in [-4, 0, 4]:
        for dy in [-4, 0, 4]:
            q = line + [dx, dy]
            ax.plot(q[:, 0], q[:, 1], color=ORANGE, lw=4, solid_capstyle="round")
    ax.set(xlim=(0, 4), ylim=(0, 4), aspect="equal", title=title)
    ax.set_xticks([])
    ax.set_yticks([])


def network(data):
    fig, grid = plt.subplots(2, 2, figsize=(10.7, 8.4), layout="constrained")
    axes = grid.ravel()
    for ax, p0, title in zip(
        axes[:2],
        [3.65, 3.84],
        ["A  Same topology, p0 = 3.65", "B  Same topology, p0 = 3.84"],
    ):
        draw_mesh(ax, data, np.array(data["states"][str(p0)]), title)
    p = np.array([r["p0"] for r in data["static"]])
    for key, label, color in [
        ("area_rms", r"Area residual RMS", BLUE),
        ("perimeter_rms", r"Perimeter residual RMS", ORANGE),
    ]:
        axes[3].semilogy(
            p,
            np.maximum([r[key] for r in data["static"]], 1e-12),
            "o-",
            color=color,
            label=label,
            markersize=4,
        )
    x = np.array(data["states"]["3.75"])
    row = next(r for r in data["static"] if r["p0"] == 3.75)
    mode = np.array(row["softest_mode"])
    draw_mesh(axes[2], data, x, "C  Softest mode, p0 = 3.75")
    axes[2].quiver(
        *(x % 4).T,
        *mode.T,
        color=BLUE,
        angles="xy",
        scale_units="xy",
        scale=1.2,
        width=0.005,
    )
    axes[3].set(
        xlabel=r"Target perimeter, $p_0$",
        ylabel="Residual",
        title="D  Incompatible to compatible",
    )
    axes[3].legend()
    save(fig, "rigidity", "unified_network")


def stiffness(data):
    rows = data["static"]
    p = np.array([r["p0"] for r in rows])
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.1), layout="constrained")
    for key, label, color in [
        ("h_min", r"$H$", BLUE),
        ("gram_min", r"$J^{\mathrm{T}}J$", GREEN),
        ("s_min", r"$S$", ORANGE),
    ]:
        y = np.array([r[key] for r in rows])
        y[abs(y) < 1e-9] = 0
        axes[0].plot(p, y, "o-", color=color, label=label, markersize=4)
    axes[0].axhline(0, color=INK, lw=0.8, alpha=0.4)
    axes[0].set_yscale("symlog", linthresh=1e-5)
    axes[0].set(
        xlabel=r"$p_0$",
        ylabel="Lowest non-translational eigenvalue",
        title="A  Three different spectra",
    )
    axes[0].legend()
    picked = [rows[i] for i in [0, 2, 4]]
    loc = np.arange(3)
    gram = np.array([r["mode_gram"] for r in picked])
    prestress = np.array([r["mode_stress"] for r in picked])
    axes[1].bar(
        loc,
        gram,
        width=0.6,
        color=GREEN,
        label=r"$u_{\min}^{\mathrm{T}}J^{\mathrm{T}}Ju_{\min}$",
    )
    axes[1].bar(
        loc,
        prestress,
        bottom=gram,
        width=0.6,
        color=ORANGE,
        label=r"$u_{\min}^{\mathrm{T}}Su_{\min}$",
    )
    axes[1].set(
        xticks=loc,
        xticklabels=[str(r["p0"]) for r in picked],
        xlabel=r"$p_0$",
        ylabel="Stiffness along the same mode",
        title="B  What restores the soft mode?",
    )
    axes[1].legend()
    axes[2].plot(
        p, [r["G_affine"] for r in rows], "o-", color=INK, label="Affine", markersize=4
    )
    axes[2].plot(
        p,
        [r["G_relaxed"] for r in rows],
        "o-",
        color=BLUE,
        label="Relaxed (Hessian)",
        markersize=4,
    )
    axes[2].plot(
        p,
        [r["G_finite_difference"] for r in rows],
        "x",
        color=ORANGE,
        ms=8,
        label="Relaxed (finite difference)",
    )
    axes[2].set(
        xlabel=r"$p_0$", ylabel="Shear modulus", title="C  Internal relaxation matters"
    )
    axes[2].legend()
    save(fig, "rigidity", "unified_hessian")


def t1(data):
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.1), layout="constrained")
    draw_mesh(
        axes[0],
        data,
        np.array(data["states"]["3.75"]),
        "A  Selected edge (2, 3)",
        labels=True,
    )
    colors = [INK, BLUE, ORANGE, GREEN, "#916fa3"]
    for row, color in zip(data["paths"], colors):
        s = np.array(row["s"])
        energy = np.array(row["energy"]) - row["initial_energy"]
        axes[1].plot(s, energy, color=color, label=rf"$p_0={row['p0']:g}$")
    axes[1].axvline(0, color=INK, lw=1, ls=":")
    axes[1].set(
        xlabel="Signed edge length (old > 0)",
        ylabel=r"$E(s)-E_{\rm initial}$",
        title="B  Constrained relaxed paths",
    )
    axes[1].legend()
    p = np.array([r["p0"] for r in data["paths"]])
    barrier = np.array([r["path_barrier"] for r in data["paths"]])
    axes[2].semilogy(p, np.maximum(barrier, 1e-8), "o-", color=BLUE, markersize=5)
    low = barrier < 1e-8
    axes[2].plot(p[low], np.full(low.sum(), 1e-8), "v", color=ORANGE, ms=7)
    axes[2].axhline(1e-8, color=INK, ls=":", lw=1)
    axes[2].set(
        xlabel=r"$p_0$",
        ylabel="Barrier along this path",
        title="C  One edge, one preparation",
    )
    axes[2].text(
        0.98,
        0.15,
        r"Below $10^{-8}$: unresolved",
        transform=axes[2].transAxes,
        ha="right",
        fontsize=10.5,
    )
    save(fig, "rearrangements", "unified_t1")


def dynamics(data):
    d = data["dynamics"]
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2), layout="constrained")
    for key, color, ls, label in [
        ("passive_fixed", INK, "-", "Passive, fixed topology"),
        ("passive_t1", ORANGE, "--", "Passive, T1 allowed"),
    ]:
        row = d[key]
        axes[0].plot(row["time"], row["response"], color=color, ls=ls, label=label)
    active = np.array([r["response"] for r in d["active"]])
    t = np.array(d["active"][0]["time"])
    mean = active.mean(axis=0)
    sem = active.std(axis=0, ddof=1) / np.sqrt(len(active))
    axes[0].plot(t, mean, color=BLUE, label="Active + T1 (mean)")
    axes[0].fill_between(t, mean - sem, mean + sem, color=BLUE, alpha=0.18, linewidth=0)
    axes[0].set(
        xlabel="Time",
        ylabel=r"$R_{\gamma}(t)$, $\gamma=0.15$",
        title="A  Symmetric finite-step response",
    )
    axes[0].legend(loc="lower left", fontsize=9.5)
    events = np.array([r["events"] for r in d["active"]]) / 16
    mean = events.mean(axis=0)
    sem = events.std(axis=0, ddof=1) / np.sqrt(len(events))
    axes[1].plot(t, mean, color=BLUE, label="Active + T1")
    axes[1].fill_between(t, mean - sem, mean + sem, color=BLUE, alpha=0.18, linewidth=0)
    axes[1].plot(
        d["passive_t1"]["time"],
        np.array(d["passive_t1"]["events"]) / 16,
        color=ORANGE,
        ls="--",
        label="Passive, T1 allowed",
    )
    axes[1].set(
        xlabel="Time", ylabel="Cumulative T1 / cell", title="B  Events in the same runs"
    )
    axes[1].legend()
    row = d["passive_fixed"]
    axes[2].plot(row["time"], row["response"], color=INK, label=r"$\gamma=0.15$")
    for row, color in zip(d["amplitude_checks"], [GREEN, BLUE]):
        axes[2].plot(
            row["time"],
            row["response"],
            color=color,
            label=rf"$\gamma={row['gamma']:g}$",
        )
    geq = next(r["G_relaxed"] for r in data["static"] if r["p0"] == 3.75)
    axes[2].axhline(geq, color=ORANGE, ls=":", label=r"$G_{\rm eq}^{\rm(top)}$")
    axes[2].set(
        xlabel="Time",
        ylabel=r"Passive $R_{\gamma}(t)$",
        title="C  A finite step is not a modulus",
        ylim=(0.03, 0.11),
    )
    axes[2].legend()
    save(fig, "rheology_activity", "unified_dynamics")


if __name__ == "__main__":
    data = json.loads((HERE / "unified_vertex_data.json").read_text())
    network(data)
    stiffness(data)
    t1(data)
    dynamics(data)
    print("Plotted four figures from the saved periodic-network calculations.")
