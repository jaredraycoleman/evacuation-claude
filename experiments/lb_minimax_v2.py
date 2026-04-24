"""Improved numerical minimax: for 3 robots on boundary with |E| = 3x,
find the configuration minimizing sup_{p in E^c} max_j ||R_j - p||.

Uses an ANALYTIC sup (not grid sup): since robot positions and scan arcs
are arcs on the circle, max chord is piecewise-smooth in beta and the
maxima are at arc boundaries OR at interior critical points (which for
chord = 2|sin((r - beta)/2)| are at beta = r + pi, i.e., antipodal to
some robot -- if antipodal is in unexplored).

For a given config, the sup is therefore max over a FINITE set:
    - Unexplored arc endpoints
    - Antipodes of robots that land in unexplored
"""
from __future__ import annotations

import math
from itertools import product

import numpy as np
from scipy.optimize import differential_evolution

TWO_PI = 2.0 * math.pi


def chord_angle(a: float, b: float) -> float:
    """Chord length for boundary angles a, b (shorter arc)."""
    d = (b - a) % TWO_PI
    if d > math.pi:
        d = TWO_PI - d
    return 2.0 * math.sin(d / 2.0)


def in_arc_halfopen(theta: float, a: float, b_plus_len: float) -> bool:
    """Is theta in the half-open arc (a, b_plus_len) going CCW on S^1?

    Here a in [0, 2 pi) and arc has length b_plus_len - a > 0, possibly
    wrapping past 2 pi.
    """
    t = (theta - a) % TWO_PI
    return 0 < t < (b_plus_len - a)


def sup_chord_in_unexplored(
    robot_angles: list[float],
    scan_arcs: list[tuple[float, float]],  # (lo, hi) with hi = lo + L, hi may > 2pi
) -> float:
    """Compute sup over unexplored of max_j chord(R_j, beta) analytically.

    Unexplored endpoints: all scan_arc endpoints.
    Interior critical points: antipodes of robots (beta = R_j + pi mod 2pi),
    when these fall in unexplored.
    """
    # Collect candidate beta values: scan-arc endpoints + robot antipodes
    candidates: list[float] = []
    for lo, hi in scan_arcs:
        candidates.append(lo % TWO_PI)
        candidates.append(hi % TWO_PI)
    for r in robot_angles:
        candidates.append((r + math.pi) % TWO_PI)

    # Filter: only keep those in UNEXPLORED (outside all scan arcs).
    def in_explored(theta: float) -> bool:
        for lo, hi in scan_arcs:
            if in_arc_halfopen(theta, lo % TWO_PI, hi):
                return True
        return False

    # For endpoints of arcs (on the boundary of explored), include them
    # (they're limits of unexplored points; chord there is valid as sup).
    best_chord = 0.0
    for beta in candidates:
        if in_explored(beta):
            continue  # strictly in explored: adversary can't pick
        c = max(chord_angle(r, beta) for r in robot_angles)
        if c > best_chord:
            best_chord = c
    return best_chord


def min_refined_chord_over_boundary_configs(x: float, n_restarts: int = 40) -> tuple[float, dict]:
    """
    For 3 robots on boundary with 3 disjoint scan arcs, each of length L = x
    (so |E| = 3x), minimize sup chord over configurations.

    Parameterization: (alpha_1, alpha_2, alpha_3, s_1, s_2, s_3).
    Constraint: arcs disjoint. We enforce by normalizing.
    """
    L = x

    def arcs_from_params(params):
        alphas = params[0:3]
        ss = np.sign(params[3:6])
        ss[ss == 0] = 1.0
        arcs = []
        robot_angles = []
        for j in range(3):
            if ss[j] > 0:
                lo = alphas[j] % TWO_PI
                hi = lo + L
                ra = lo + L
            else:
                lo = (alphas[j] - L) % TWO_PI
                hi = lo + L
                ra = lo  # robot at scan end (= lower bound when CW)
            arcs.append((lo, hi))
            robot_angles.append(ra % TWO_PI)
        return robot_angles, arcs

    def disjoint(arcs) -> float:
        """Return overlap amount (0 if disjoint)."""
        total = 0.0
        # Normalize: convert each arc to a set on [0, 2 pi) possibly wrap
        sets: list[list[tuple[float, float]]] = []
        for lo, hi in arcs:
            if hi <= TWO_PI:
                sets.append([(lo, hi)])
            else:
                sets.append([(lo, TWO_PI), (0.0, hi - TWO_PI)])
        # Compute pairwise overlap
        for i in range(3):
            for j in range(i + 1, 3):
                for (a1, b1) in sets[i]:
                    for (a2, b2) in sets[j]:
                        total += max(0.0, min(b1, b2) - max(a1, a2))
        return total

    def objective(params):
        robot_angles, arcs = arcs_from_params(params)
        overlap = disjoint(arcs)
        chord = sup_chord_in_unexplored(robot_angles, arcs)
        return chord + 100.0 * overlap

    # Optimize
    bounds = [(0.0, TWO_PI)] * 3 + [(-1.0, 1.0)] * 3
    best = float("inf")
    best_info = None
    for trial in range(n_restarts):
        res = differential_evolution(
            objective, bounds, seed=trial, maxiter=300, popsize=40, tol=1e-10,
            init="sobol", polish=True, workers=1,
        )
        if res.fun < best:
            best = res.fun
            ra, arcs = arcs_from_params(res.x)
            best_info = {"chord": res.fun, "robot_angles": ra, "arcs": arcs, "params": res.x}
    return best, best_info


def main():
    x_cgk = (2.0 / 3.0) * math.acos(-1.0 / 3.0)
    print(f"x_CGK = {x_cgk:.8f}")
    print()

    # Symbolic prediction for 3-fold symmetric: 2 sin(pi/3 + x/2)
    pred_sym = 2 * math.sin(math.pi / 3 + x_cgk / 2)
    print(f"Prediction (3-fold symmetric): {pred_sym:.8f}")
    print()

    # Run numerical minimax with analytic sup
    print("Running numerical minimax with analytic sup (40 restarts)...")
    chord_min, info = min_refined_chord_over_boundary_configs(x_cgk, n_restarts=40)
    print(f"Min sup chord found: {chord_min:.8f}")
    print(f"Refined LB = {1 + x_cgk + chord_min:.8f}")
    print()
    print("Best configuration:")
    print(f"  robot angles: {info['robot_angles']}")
    print(f"  scan arcs:    {info['arcs']}")
    print()

    # Check a few specific configurations
    print("Sanity checks on specific configurations:")

    # 3-fold symmetric CCW
    ra = [L * j_factor for j_factor, L in zip([0, 0, 0], [0, 0, 0])]  # placeholder
    alphas = [0, TWO_PI / 3, 2 * TWO_PI / 3]
    ra = [(a + x_cgk) % TWO_PI for a in alphas]
    arcs = [(a, a + x_cgk) for a in alphas]
    c = sup_chord_in_unexplored(ra, arcs)
    print(f"  3-fold symmetric CCW, L=x: sup chord = {c:.8f}")

    # Single-arc [0, 3x]
    arcs = [(0, x_cgk), (x_cgk, 2 * x_cgk), (2 * x_cgk, 3 * x_cgk)]
    ra = [x_cgk, 2 * x_cgk, 3 * x_cgk]
    c = sup_chord_in_unexplored(ra, arcs)
    print(f"  Single-arc [0, 3x]:        sup chord = {c:.8f}")

    # 2-adjacent + 1-opposite
    arcs = [(0, x_cgk), (x_cgk, 2 * x_cgk), (math.pi, math.pi + x_cgk)]
    ra = [x_cgk, 2 * x_cgk, (math.pi + x_cgk) % TWO_PI]
    c = sup_chord_in_unexplored(ra, arcs)
    print(f"  2-adjacent + 1-opposite:    sup chord = {c:.8f}")


if __name__ == "__main__":
    main()
