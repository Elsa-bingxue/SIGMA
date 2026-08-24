"""Regenerate versioned SIGMA manuscript figures from a reference project."""
from __future__ import annotations

import argparse

from sigma_spatial import reproduce_manuscript_figures


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("project_root")
    parser.add_argument("--figures", nargs="*", default=None)
    args = parser.parse_args()
    for output in reproduce_manuscript_figures(args.project_root, figures=args.figures):
        print(output)


if __name__ == "__main__":
    main()
