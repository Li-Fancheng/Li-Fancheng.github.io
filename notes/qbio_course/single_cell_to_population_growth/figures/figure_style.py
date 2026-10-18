"""Use the course-wide Arial / Computer Modern plotting and export settings."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from plot_style import configure_style, save_figure
