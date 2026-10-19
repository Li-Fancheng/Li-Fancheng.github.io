"""Regenerate every figure from the source code.

From the directory containing README.md:
    python code/make_all_figures.py
"""
from random_walk import make_random_walk_figures
from fixed_field import make_fixed_field_figure
from feedback_diagram import make_feedback_figure
from collective_1d import make_collective_figures
from collective_2d import make_2d_figure
from teaching_diagrams import make_teaching_diagrams
from chemosensing import make_chemosensing_figure


def main():
    make_random_walk_figures()
    make_fixed_field_figure()
    make_feedback_figure()
    make_collective_figures()
    make_2d_figure()
    make_teaching_diagrams()
    make_chemosensing_figure()
    print('Done. All 15 model figures are in figures/; paper originals are preserved.')


if __name__ == '__main__':
    main()
