"""Run one SIGMA dataset from a portable JSON configuration."""
from __future__ import annotations

import argparse

from sigma_spatial import run_analysis_config


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="path to a SIGMA JSON configuration")
    args = parser.parse_args()
    result = run_analysis_config(args.config)
    print(result.output_dir)


if __name__ == "__main__":
    main()
