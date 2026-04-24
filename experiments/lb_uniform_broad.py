"""Broad numerical minimax for the uniform LB.

Searches over a parametric family that includes:
    - 2-coincident with asymmetric slack (3 params: tau_1, tau_2, tau_3)
    - 3-fold symmetric with inward motion
    - Off-center clusters
    - Fully-distinct robots with inward motion

Parametrization: at time 1 + x^*, each robot j has
    - scan arc angle alpha_j in [0, 2 pi) (scan start)
    - scan direction s_j in {-1, +1}
    - scan length L_j in [0, x^*]
    - inward distance tau_j in [0, x^* - L_j] (move radially inward
      from scan end)
so position R_j = (1 - tau_j) e(alpha_j + s_j L_j), on or inside disk.

Computes combined = max(state-tight CGK, refined) and minimises it.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import differential_evolution

TWO_PI = 2.0 * math.pi
X_STAR = (2.0 / 3.0) * math.acos(-1.0 / 3.0)


def state_tight_cgk(e: float) -> float:
    if e >= math.pi:
        return 1.0 + X_STAR + 2.0 * math.sin((TWO_PI - e) / 2.0)
    return 1.0 + X_STAR + 2.0


def combined_for_params(params: np.ndarray) -> float:
    alphas = params[0:3]
    ss = np.sign(params[3:6])
    ss[ss == 0] = 1.0
    Ls = np.clip(params[6:9], 0.0, X_STAR)
    max_tau = np.maximum(0.0, X_STAR - Ls)
    taus_raw = params[9:12]
    taus = np.clip(taus_raw, 0.0, np.minimum(max_tau, 1.0))

    # Robot positions
    robots = np.array([
        (1.0 - taus[j]) * np.array([math.cos(alphas[j] + ss[j] * Ls[j]),
                                    math.sin(alphas[j] + ss[j] * Ls[j])])
        for j in range(3)
    ])

    # Scan arcs (as intervals, possibly wrapping): normalize each to (lo, lo+L)
    arcs = []
    for j in range(3):
        if ss[j] > 0:
            lo = alphas[j] % TWO_PI
        else:
            lo = (alphas[j] - Ls[j]) % TWO_PI
        arcs.append((lo, lo + Ls[j]))

    # |E| = union length. Compute disjoint total minus overlaps.
    # Use a grid to compute both explored mask and sup chord.
    n_grid = 1000
    grid = np.linspace(0.0, TWO_PI, n_grid, endpoint=False)
    covered = np.zeros_like(grid, dtype=bool)
    for (lo, hi) in arcs:
        lo = lo % TWO_PI
        if hi <= TWO_PI:
            covered |= (grid >= lo) & (grid <= hi)
        else:
            covered |= (grid >= lo) & (grid <= TWO_PI)
            covered |= (grid >= 0.0) & (grid <= hi - TWO_PI)

    e_len = covered.sum() * (TWO_PI / n_grid)
    cgk = state_tight_cgk(e_len)

    # Sup chord over unexplored
    unexp_mask = ~covered
    if unexp_mask.sum() == 0:
        refined = 1.0 + X_STAR  # vacuous: no unexplored
    else:
        unexp_thetas = grid[unexp_mask]
        ps = np.stack([np.cos(unexp_thetas), np.sin(unexp_thetas)], axis=1)
        diffs = ps[None, :, :] - robots[:, None, :]
        dists = np.linalg.norm(diffs, axis=2)
        max_chord = dists.max(axis=0).max()
        refined = 1.0 + X_STAR + float(max_chord)

    return max(cgk, refined)


def main():
    print(f"x^* = {X_STAR:.8f}")
    print()

    bounds = (
        [(0.0, TWO_PI)] * 3       # alphas
        + [(-1.0, 1.0)] * 3       # signs
        + [(0.0, X_STAR)] * 3     # L_j
        + [(0.0, 1.0)] * 3        # taus (clipped internally)
    )

    best = float("inf")
    best_x = None
    for trial in range(12):
        res = differential_evolution(
            combined_for_params, bounds, seed=trial,
            maxiter=400, popsize=50, tol=1e-10,
            init="sobol", workers=1, polish=True,
        )
        if res.fun < best:
            best = res.fun
            best_x = res.x
        print(f"  trial {trial}: {res.fun:.6f}")

    print()
    print(f"Minimum combined LB (over broad alg states): {best:.8f}")
    print(f"CGK LB:                                       {1 + X_STAR + 2*math.sin(3*X_STAR/2):.8f}")
    print(f"Scan-tight refined:                           {1 + X_STAR + 2*math.cos(X_STAR/4):.8f}")
    print(f"Balanced 2-coincident opt:                    4.16702689  (from lb_uniform.py)")
    print()
    if best_x is not None:
        alphas = best_x[0:3]
        ss = np.sign(best_x[3:6])
        Ls = np.clip(best_x[6:9], 0.0, X_STAR)
        taus = np.clip(best_x[9:12], 0.0, np.maximum(0.0, X_STAR - Ls))
        print(f"Best state params:")
        print(f"  alphas = {alphas}")
        print(f"  signs  = {ss}")
        print(f"  L_j    = {Ls}")
        print(f"  tau_j  = {taus}")
        print(f"  r_j    = {1 - taus}")
        print(f"  robot angles = {alphas + ss * Ls}")


if __name__ == "__main__":
    main()
