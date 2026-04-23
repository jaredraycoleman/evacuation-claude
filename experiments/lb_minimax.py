"""Numerical minimax for the refined lower bound.

For any 3-robot wireless algorithm:
    evac >= (1+x) + max_j ||R_j(1+x) - p||   for any unexplored p.

To get a universal LB:
    LB_ref = inf_{alg} sup_{x} sup_{p in E^c(1+x)} [(1+x) + max_j ||R_j - p||].

We parametrize algorithm STATES at time 1+x, not full trajectories.
A state is consistent if it can be achieved by some valid unit-speed
trajectory from origin by time 1+x. Relaxed consistency: each robot j
has scanned an arc of length L_j <= x on boundary and then (optionally)
moved inward by r_j in [0, x - L_j] radially from scan end.

Parameters of a state at time 1+x (fix x):
    - 3 scan-start angles alpha_j in [0, 2 pi)
    - 3 scan directions s_j in {+1, -1}
    - 3 scan lengths L_j in [0, x]
    - 3 inward displacements r_j in [0, max(0, x - L_j)]

Then:
    scan_arc_j = [alpha_j, alpha_j + s_j L_j]  (as interval on S^1)
    E = union of scan arcs
    R_j = (1 - r_j) * e(alpha_j + s_j L_j)    (position at 1+x)

This is a 12-parameter algorithm-state space (for each x).

For each state, evaluate sup_{p in E^c} (1+x) + max_j ||R_j - p||,
then minimize over states (for fixed x), then maximize over x.

If the result exceeds 4.159, the CGK LB can be tightened.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import differential_evolution, minimize

TWO_PI = 2.0 * math.pi


def e(theta: float) -> np.ndarray:
    return np.array([math.cos(theta), math.sin(theta)])


def wrap(a: float) -> float:
    return a % TWO_PI


def arc_as_interval(alpha: float, s: float, L: float) -> tuple[float, float]:
    """Return normalized [a, b] on [0, 2 pi), possibly wrap-around.

    If the scan arc wraps, we represent it as (a, b) with b > 2 pi,
    or split into two arcs. For simplicity, return (a, b) with a in
    [0, 2 pi) and b = a + L (possibly > 2 pi). Caller handles wrap.
    """
    if s > 0:
        a = wrap(alpha)
        return (a, a + L)
    else:
        a = wrap(alpha - L)
        return (a, a + L)


def covered_mask(alphas: np.ndarray, ss: np.ndarray, Ls: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Return boolean array of which grid points are in the union of scan arcs."""
    covered = np.zeros_like(grid, dtype=bool)
    for j in range(3):
        a, b = arc_as_interval(alphas[j], ss[j], Ls[j])
        # Mark grid points in [a, b] mod 2 pi.
        lo, hi = a, b
        if hi <= TWO_PI:
            covered |= (grid >= lo) & (grid <= hi)
        else:
            covered |= (grid >= lo) & (grid <= TWO_PI)
            covered |= (grid >= 0.0) & (grid <= hi - TWO_PI)
    return covered


def state_refined_lb(x: float, params: np.ndarray, grid: np.ndarray) -> float:
    """
    Compute sup_{p in E^c} (1+x) + max_j ||R_j - p||  for a given state.

    params: 12 values
      params[0:3]   = alphas (in [0, 2 pi))
      params[3:6]   = signs encoded as float in [-1, 1], sign taken
      params[6:9]   = Ls in [0, x]
      params[9:12]  = rs in [0, max(0, x - L_j)]
    """
    alphas = params[0:3]
    ss = np.sign(params[3:6])
    ss[ss == 0] = 1.0
    Ls = np.clip(params[6:9], 0.0, x)
    rs_raw = params[9:12]
    max_r = np.maximum(0.0, x - Ls)
    rs = np.clip(rs_raw, 0.0, max_r)

    # Robot positions
    robots = np.array([(1.0 - rs[j]) * e(alphas[j] + ss[j] * Ls[j]) for j in range(3)])

    # Explored mask on grid
    covered = covered_mask(alphas, ss, Ls, grid)
    unexp_mask = ~covered
    if not unexp_mask.any():
        return 1.0 + x  # no unexplored - this state isn't possible for relevant x

    # For each unexplored grid point, compute max_j chord
    unexp_thetas = grid[unexp_mask]
    ps = np.stack([np.cos(unexp_thetas), np.sin(unexp_thetas)], axis=1)  # shape (n, 2)
    # Distances |R_j - p| for each j, p
    diffs = ps[None, :, :] - robots[:, None, :]  # (3, n, 2)
    dists = np.linalg.norm(diffs, axis=2)        # (3, n)
    max_chord = dists.max(axis=0)                # (n,)

    return 1.0 + x + float(max_chord.max())


def minimize_state_lb(x: float, n_grid: int = 720, n_restarts: int = 10,
                      seed: int = 0) -> tuple[float, np.ndarray]:
    """
    At fixed x, find algorithm state that minimizes the refined LB.
    Uses differential_evolution + local refinement.

    Returns (best_lb, best_params).
    """
    grid = np.linspace(0.0, TWO_PI, n_grid, endpoint=False)

    bounds = (
        [(0.0, TWO_PI)] * 3 +  # alphas
        [(-1.0, 1.0)] * 3 +    # signs
        [(0.0, x)] * 3 +        # Ls
        [(0.0, x)] * 3          # rs (clamped internally)
    )

    def objective(params):
        return state_refined_lb(x, params, grid)

    best_lb = float("inf")
    best_params = None

    for trial in range(n_restarts):
        result = differential_evolution(
            objective, bounds,
            seed=seed + trial,
            maxiter=150,
            popsize=30,
            tol=1e-8,
            polish=True,
            init="sobol",
            workers=1,
        )
        if result.fun < best_lb:
            best_lb = result.fun
            best_params = result.x.copy()

    return best_lb, best_params


def main():
    x_cgk = (2.0 / 3.0) * math.acos(-1.0 / 3.0)
    cgk_lb = 1.0 + x_cgk + 2.0 * math.sin(3 * x_cgk / 2)
    print(f"CGK+14: x = {x_cgk:.6f}, LB = {cgk_lb:.6f}")
    print()

    # Quick scan - high grid res, moderate restarts
    xs_to_probe = [1.0, 1.1, 1.2, 1.274, 1.35, 1.4, 1.486, 1.55, 1.65, 1.8, 1.95]

    print("  x       | min refined LB  | CGK LB      | diff ")
    print("  --------+-----------------+-------------+--------")
    overall_max = -float("inf")
    overall_x = None
    for x in xs_to_probe:
        lb, params = minimize_state_lb(x, n_grid=720, n_restarts=6)
        cgk = 1 + x + 2 * math.sin(3 * x / 2) if math.pi / 3 <= x <= TWO_PI / 3 else float("nan")
        diff = lb - cgk if not math.isnan(cgk) else float("nan")
        print(f"  {x:.3f}   | {lb:.7f}        | {cgk:.7f}   | {diff:+.5f}")
        if lb > overall_max:
            overall_max = lb
            overall_x = x

    print()
    print(f"Universal refined LB (numerical, over probed x):")
    print(f"  max_x min_alg sup_p LB = {overall_max:.7f} at x = {overall_x}")
    print(f"  CGK LB                 = {cgk_lb:.7f}")
    print(f"  improvement            = {overall_max - cgk_lb:+.5f}")


if __name__ == "__main__":
    main()
