"""Typography and vector export shared by the quantitative biology notes.

Text uses Arial where available, with portable sans-serif fallbacks. Mathematics
uses Computer Modern, matching the LaTeX-style equations on the website. SVG
glyph outlines preserve both fonts without requiring fonts on readers' devices.
The display widths in data/qbio_figures.json also determine the smallest label
size, so exporting a large canvas does not create tiny lettering on the page.
"""

import json
from pathlib import Path
import re

import matplotlib as mpl
from matplotlib import font_manager
from matplotlib.text import Text
from matplotlib.ticker import AutoLocator

COURSE = Path(__file__).resolve().parent
POLICIES = json.loads((COURSE.parents[1] / "data/qbio_figures.json").read_text())


def configure_style(*, line_width=2.3, text_color="#243747"):
    # macOS does not always include Supplemental fonts in its font cache.
    for directory in (Path("/System/Library/Fonts/Supplemental"),):
        for name in ("Arial.ttf", "Arial Bold.ttf", "Arial Italic.ttf", "Arial Bold Italic.ttf"):
            candidate = directory / name
            if candidate.exists():
                font_manager.fontManager.addfont(candidate)
    installed = {font.name for font in font_manager.fontManager.ttflist}
    sans = [name for name in ("Arial", "Liberation Sans", "DejaVu Sans") if name in installed]
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": sans,
        "font.size": 15,
        "font.weight": "normal",
        "mathtext.fontset": "cm",
        "mathtext.default": "it",
        "text.usetex": False,
        "axes.titlesize": 17,
        "axes.titleweight": "bold",
        "axes.labelsize": 15,
        "axes.labelweight": "normal",
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 14,
        "figure.titlesize": 18,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "text.color": text_color,
        "axes.labelcolor": text_color,
        "axes.titlecolor": text_color,
        "legend.frameon": False,
        "lines.linewidth": line_width,
        "svg.fonttype": "path",
        "svg.hashsalt": "qbio-course-figures",
    })


def _wrap(text, renderer, max_width):
    """Wrap prose without splitting a $...$ expression or an explicit line."""
    if not text.get_text():
        return
    original = getattr(text, "_qbio_original", text.get_text())
    text._qbio_original = original
    lines = []
    for paragraph in original.split("\n"):
        current = ""
        for word in re.findall(r"\$[^$]*\$|[^\s$]+", paragraph):
            candidate = f"{current} {word}".strip()
            width, _, _ = renderer.get_text_width_height_descent(
                candidate, text.get_fontproperties(), ismath="$" in candidate)
            if current and width > max_width:
                lines.append(current)
                current = word
            else:
                current = candidate
        lines.append(current)
    text.set_text("\n".join(lines))


def save_figure(fig, directory, name, *, min_label_pixels=13):
    directory = Path(directory)
    key = str((directory / name).relative_to(COURSE))
    policy = POLICIES[key]
    width = policy["width"]
    family = mpl.rcParams["font.sans-serif"][0]

    # Increase explicit small labels as well as defaults. Keep larger formulae
    # and annotations large; these sizes are measured at the intended web width.
    floor = min_label_pixels * fig.get_figwidth() * 72 / width
    for text in fig.findobj(Text):
        text.set_fontfamily(family)
        text.set_math_fontfamily("cm")
        if text.get_fontweight() == "medium":
            text.set_fontweight("normal")
        elif text.get_fontweight() == "semibold":
            text.set_fontweight("bold")
        text.set_fontsize(max(text.get_fontsize(), floor))
    for ax in fig.axes:
        if ax.get_aspect() != "auto":
            ax.set_anchor("N")
        legend = ax.get_legend()
        if legend is not None:
            legend.set_frame_on(True)
            legend.get_frame().set_facecolor("white")
            legend.get_frame().set_edgecolor("none")
            legend.get_frame().set_alpha(.94)
        for axis in (ax.xaxis, ax.yaxis):
            locator = axis.get_major_locator()
            if isinstance(locator, AutoLocator):
                locator.set_params(nbins=5)

    for _ in range(2):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        for ax in fig.axes:
            for title in (ax.title, ax._left_title, ax._right_title):
                _wrap(title, renderer, ax.get_window_extent(renderer).width * .98)
            if ax.get_xlabel():
                _wrap(ax.xaxis.label, renderer, ax.get_window_extent(renderer).width * 1.05)
        if fig._suptitle:
            _wrap(fig._suptitle, renderer, fig.bbox.width * .94)
        if fig.get_layout_engine() is None and all(ax.axison for ax in fig.axes):
            fig.tight_layout(pad=1.2, w_pad=2.2, h_pad=2.2)

    # Give each panel in a row the same title depth. A short title should not
    # sit one line below its wrapped neighbour.
    rows = {}
    for ax in fig.axes:
        spec = ax.get_subplotspec()
        if spec is not None and ax.get_label() != "<colorbar>":
            rows.setdefault((id(spec.get_gridspec()), spec.rowspan.start), []).append(ax)
    for axes in rows.values():
        titles = [title for ax in axes for title in (ax.title, ax._left_title, ax._right_title)
                  if title.get_text()]
        depth = max((title.get_text().count("\n") for title in titles), default=0)
        for title in titles:
            title.set_text(title.get_text() + "\n" * (depth - title.get_text().count("\n")))

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    # Include export padding in the screen-size check, rather than assuming
    # that the original canvas width is the final SVG viewBox width.
    export_width = (fig.get_tightbbox(renderer).width + .28) * 72
    actual_min = min(t.get_fontsize() * width / export_width
                     for t in fig.findobj(Text) if t.get_visible() and t.get_text())
    if actual_min < min_label_pixels - .4:
        factor = min_label_pixels / actual_min
        for text in fig.findobj(Text):
            text.set_fontsize(text.get_fontsize() * factor)
        fig.canvas.draw()

    svg = directory / f"{name}.svg"
    fig.savefig(svg, bbox_inches="tight", pad_inches=.14, metadata={"Date": None})
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
    fig.savefig(directory / f"{name}.png", dpi=300, bbox_inches="tight", pad_inches=.14)
