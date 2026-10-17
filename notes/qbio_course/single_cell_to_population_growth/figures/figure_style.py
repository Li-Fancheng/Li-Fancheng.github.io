"""Shared typography and export settings for the lecture's generated figures.

STIX is bundled with Matplotlib, so text and mathematics do not depend on
system fonts or a TeX installation. SVG glyph outlines preserve that typography
in browsers; 300-dpi PNGs remain available for raster workflows.
"""

from pathlib import Path

import matplotlib as mpl


def configure_style(*, line_width=2.4, text_color="black"):
    mpl.rcParams.update({
        "font.family": "STIXGeneral",
        "font.size": 12,
        "mathtext.fontset": "stix",
        "mathtext.default": "it",
        "text.usetex": False,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "text.color": text_color,
        "axes.labelcolor": text_color,
        "legend.frameon": False,
        "lines.linewidth": line_width,
        "svg.fonttype": "path",
        "svg.hashsalt": "cell-cycle-lecture",
    })


def save_figure(fig, directory: Path, name: str):
    fig.savefig(directory / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(directory / f"{name}.svg", bbox_inches="tight",
                metadata={"Date": None})
