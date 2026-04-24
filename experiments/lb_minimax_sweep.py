"""Sweep x over valid range and find min sup chord over 3-robot-on-boundary
configurations with |E| = 3x.

Compare against candidate formulas:
    - 2 cos(x/4)  (midpoint-of-unexplored chord in 3-fold symmetric)
    - 2 sin(pi/3 + x/2)  (endpoint chord in 3-fold symmetric)
    - 2 sin(3x/2)  (CGK's chord across unexplored arc)
"""
from __future__ import annotations

import math

import numpy as np

import sys
sys.path.insert(0, ".")
from experiments.lb_minimax_v2 import min_refined_chord_over_boundary_configs

TWO_PI = 2.0 * math.pi


def main():
    print("Sweep of min sup chord over 3-robot boundary configurations with |E| = 3x")
    print()
    print("  x        | min sup chord | 2 cos(x/4) | 2 sin(pi/3+x/2) | 2 sin(3x/2) |")
    print("  ---------+---------------+-------------+-----------------+--------------")
    for x in np.linspace(0.9, 2.0, 12):
        chord_min, info = min_refined_chord_over_boundary_configs(x, n_restarts=12)
        c_sym_mid = 2 * math.cos(x / 4)
        c_sym_end = 2 * math.sin(math.pi / 3 + x / 2)
        c_cgk = 2 * math.sin(3 * x / 2) if x <= 2 * math.pi / 3 else float("nan")
        print(f"  {x:.4f}   | {chord_min:.7f}     | {c_sym_mid:.7f}   | {c_sym_end:.7f}       | {c_cgk:.7f}     |")


if __name__ == "__main__":
    main()
