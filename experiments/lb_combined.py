"""Combined lower bound via max(state-tight-CGK, refined).

For an algorithm with state (R_j, E) at time 1+x:

    state_tight_CGK(state, x) = 1 + x + 2 sin((2 pi - |E|) / 2)  if |E| >= pi
                              = 1 + x + 2                         if |E| < pi
    refined(state, x)          = (1+x) + sup_{p in E^c} max_j ||R_j - p||

    combined(state, x) = max(state_tight_CGK, refined)

Universal LB: inf_{valid state at 1+x} combined(state, x), maximised over x.

The CONSERVATIVE CGK published at 4.159 uses state_tight_CGK with |E| = 3x
(max coverage). state_tight_CGK with |E| < 3x is strictly larger, so the
min-state always has |E| = 3x at x = x_{CGK}^* = 1.274. Question: for any
state with |E| = 3x, is refined > 4.159?

If yes: universal LB > 4.159.
If no: find a "bad" alg state where both bounds are <= 4.159.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import differential_evolution

TWO_PI = 2.0 * math.pi


def state_tight_cgk(E_len: float, x: float) -> float:
    """Alg-specific CGK LB: 1 + x + chord across unexplored set."""
    if E_len >= math.pi:
        # shorter arc = 2 pi - |E| <= pi; chord = 2 sin(shorter/2)
        u = TWO_PI - E_len
        return 1.0 + x + 2.0 * math.sin(u / 2.0)
    else:
        # unexplored length > pi; exist antipodal unexplored pts
        return 1.0 + x + 2.0


# --------------------------------------------------------------------------
# Parameterise: 3 robots on boundary, with 3 disjoint scan arcs of total
# length 3x. Each robot is at end of scan arc.
#
# Parameters (6): arc start angles alpha_j in [0, 2 pi), scan directions
# s_j in {-1, +1} (encoded as sign of a float).
#
# For disjoint arcs totaling 3x: arcs are determined by starts and
# lengths. We'll assume each robot scans the SAME length L = x, and
# the 3 starts are at alpha_1, alpha_2, alpha_3. Constraint: arcs
# disjoint.
#
# Robot positions: R_j = e(alpha_j + s_j x).
# Explored: union of [alpha_j, alpha_j + s_j x] (absolute).
# --------------------------------------------------------------------------


def arc_of(alpha: float, s: float, L: float) -> tuple[float, float]:
    """Normalize arc to (lo, hi) with lo in [0, 2pi), hi = lo + L, possibly > 2pi."""
    if s > 0:
        lo = alpha % TWO_PI
    else:
        lo = (alpha - L) % TWO_PI
    return (lo, lo + L)


def in_arc(theta: float, arc: tuple[float, float]) -> bool:
    lo, hi = arc
    theta = theta % TWO_PI
    if hi <= TWO_PI:
        return lo <= theta <= hi
    return theta >= lo or theta <= hi - TWO_PI


def arcs_overlap(a1: tuple[float, float], a2: tuple[float, float]) -> bool:
    """Check if two arcs on circle overlap (strict interior overlap)."""
    # Sample test
    samples1 = np.linspace(a1[0], a1[1], 50)[1:-1]  # interior
    for t in samples1:
        if in_arc(t % TWO_PI, a2):
            return True
    samples2 = np.linspace(a2[0], a2[1], 50)[1:-1]
    for t in samples2:
        if in_arc(t % TWO_PI, a1):
            return True
    return False


def refined_on_boundary_state(alphas: np.ndarray, ss: np.ndarray, L: float,
                              grid: np.ndarray) -> float:
    """
    Compute (1+x) + max over unexplored p of max_j ||R_j - p||,
    for a state with 3 robots on boundary, scans length L each, disjoint.

    Returns sup over boundary unexplored points, as a CHORD value (not LB).
    """
    # Robot angles
    robot_angles = [alphas[j] + ss[j] * L for j in range(3)]
    robots = np.array([[math.cos(a), math.sin(a)] for a in robot_angles])

    # Explored mask on grid
    covered = np.zeros_like(grid, dtype=bool)
    for j in range(3):
        a_start = alphas[j]
        a_end = alphas[j] + ss[j] * L
        lo = min(a_start, a_end) % TWO_PI
        hi = (min(a_start, a_end) + L) % (TWO_PI + 0.0)  # preserve
        if ss[j] > 0:
            lo_ = a_start % TWO_PI
            hi_ = lo_ + L
        else:
            lo_ = (a_start - L) % TWO_PI
            hi_ = lo_ + L
        if hi_ <= TWO_PI:
            covered |= (grid >= lo_) & (grid <= hi_)
        else:
            covered |= (grid >= lo_) & (grid <= TWO_PI)
            covered |= (grid >= 0.0) & (grid <= hi_ - TWO_PI)

    unexp_mask = ~covered
    if not unexp_mask.any():
        return 0.0

    unexp_thetas = grid[unexp_mask]
    ps = np.stack([np.cos(unexp_thetas), np.sin(unexp_thetas)], axis=1)
    diffs = ps[None, :, :] - robots[:, None, :]
    dists = np.linalg.norm(diffs, axis=2)
    max_chord = dists.max(axis=0)
    return float(max_chord.max())


def min_refined_over_states(x: float, L: float, n_grid: int = 720,
                            n_restarts: int = 10) -> tuple[float, np.ndarray]:
    """
    Minimize refined chord over alg states with all-on-boundary, scan L each,
    disjoint arcs.

    For disjoint arcs of length L: 3L <= 2 pi, so L <= 2 pi / 3.
    We optimize over (alpha_1, alpha_2, alpha_3, s_1, s_2, s_3).
    """
    grid = np.linspace(0.0, TWO_PI, n_grid, endpoint=False)

    bounds = [(0.0, TWO_PI)] * 3 + [(-1.0, 1.0)] * 3

    def objective(params):
        alphas = params[0:3]
        ss = np.sign(params[3:6])
        ss[ss == 0] = 1.0
        # Check arcs are disjoint (approximately)
        arcs = [arc_of(alphas[j], ss[j], L) for j in range(3)]
        # Penalize overlap strongly by adding a large value
        penalty = 0.0
        for i in range(3):
            for jj in range(i + 1, 3):
                if arcs_overlap(arcs[i], arcs[jj]):
                    penalty += 10.0
        chord = refined_on_boundary_state(alphas, ss, L, grid)
        return chord + penalty

    best = float("inf")
    best_params = None
    for t in range(n_restarts):
        res = differential_evolution(
            objective, bounds, seed=t, maxiter=100, popsize=25, tol=1e-8,
            init="sobol", workers=1, polish=True,
        )
        if res.fun < best:
            best, best_params = res.fun, res.x
    return best, best_params


def main():
    x_cgk = (2.0 / 3.0) * math.acos(-1.0 / 3.0)
    cgk_universal = 1 + x_cgk + 2 * math.sin(3 * x_cgk / 2)
    print(f"CGK+14 universal LB (|E| = 3x tight): {cgk_universal:.6f} at x = {x_cgk:.5f}")
    print()

    # Scan: for each x, set L = x (full scan), find min of refined over
    # boundary configs; this is the "min refined when CGK is tight".
    print("For algs with |E| = 3x (all-on-boundary, full scan):")
    print("  x       |  L = x  | min refined chord | refined LB   | CGK LB      | combined min")
    print("  --------+---------+-------------------+--------------+-------------+--------------")

    best_combined = 0.0
    best_x = None
    for x in [1.0, 1.1, 1.2, 1.274, 1.35, 1.4, 1.486, 1.55, 1.65, 1.8]:
        L = x
        if L > TWO_PI / 3 - 1e-6:
            L = TWO_PI / 3 - 1e-6
        chord, params = min_refined_over_states(x, L, n_grid=720, n_restarts=4)
        refined_lb = 1 + x + chord
        cgk_at_x = 1 + x + 2 * math.sin(3 * x / 2)
        combined = max(cgk_at_x, refined_lb)
        print(f"  {x:.3f}   | {L:.4f}  |  {chord:.5f}          | {refined_lb:.6f}    | {cgk_at_x:.6f}   | {combined:.6f}")
        if combined > best_combined:
            best_combined = combined
            best_x = x

    print()
    print(f"Best combined LB (over x, all-on-boundary full-scan algs):")
    print(f"  combined = {best_combined:.6f} at x = {best_x}")
    print(f"  CGK universal = {cgk_universal:.6f}")
    print(f"  improvement   = {best_combined - cgk_universal:+.5f}")


if __name__ == "__main__":
    main()
