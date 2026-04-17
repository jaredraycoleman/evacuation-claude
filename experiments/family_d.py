"""Family D: 3-fold symmetric, each robot has a two-phase trajectory.

All three robots do the same thing rotated by 2*pi/3. Per robot:
    approach alpha_i -> scan of length a (CCW) -> return to origin
    -> approach alpha_i + 2*pi/3 + delta -> scan of length b (direction sb)
with alpha_i = 2*pi*i/3. For full 3-fold symmetry in deployment, the
second-phase approach angle is offset by 2*pi/3 + delta from the first.

Coverage: by symmetry, the three arcs of each type are rotations of each
other, so if a single robot's two arcs together cover an angular measure
of a+b in its "sector" and the sectors tile, total coverage is 3(a+b).
We need 3(a+b) >= 2*pi, and the specific placement of the second arcs
must not leave gaps.

Free parameters (by symmetry):
    a         first-scan length (CCW)
    b         second-scan length
    sb in {+,-} second-scan direction
    delta     offset of second approach within its sector
    sa in {+,-} first-scan direction (we fix sa = +1 WLOG by reflection)

Total continuous: 3 (a, b, delta). Discrete: sb. 2 combos to try.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import differential_evolution, minimize

from experiments.wireless_k3_opt2 import (
    _robot_phases_two,
    worst_case_time,
)
from experiments.family_c import _arc_union_gap

TWO_PI = 2.0 * np.pi
EVAC_A3 = 4.2185169939


def robots_from(a, b, delta, sb):
    """Build the 3 robots with 3-fold symmetry. Robot i starts by
    approaching 2*pi*i/3, scans CCW by a, returns to origin, then
    approaches 2*pi*i/3 + 2*pi/3 + delta (i.e. into the next sector
    plus an offset delta), then scans by b in direction sb."""
    robots = []
    for i in range(3):
        alpha = (TWO_PI * i / 3.0) % TWO_PI
        beta = (alpha + TWO_PI / 3.0 + delta) % TWO_PI
        r = _robot_phases_two(a=alpha, La=a, sa=+1, b=beta, Lb=b, sb=sb)
        robots.append(r)
    return robots


def coverage_gap_exact(a, b, delta, sb):
    arcs = []
    for i in range(3):
        alpha = (TWO_PI * i / 3.0) % TWO_PI
        beta = (alpha + TWO_PI / 3.0 + delta) % TWO_PI
        arcs.append((alpha, a))                              # first scan CCW
        if sb > 0:
            arcs.append((beta, b))                           # second scan CCW
        else:
            arcs.append(((beta - b) % TWO_PI, b))            # second scan CW
    return _arc_union_gap(arcs)


def objective(x, sb: int = +1):
    a, b, delta = x
    if a < 0 or b < 0 or a > TWO_PI or b > TWO_PI:
        return 1000.0
    gap = coverage_gap_exact(a, b, delta, sb)
    if gap > 1e-9:
        return 500.0 + 100.0 * gap
    robots = robots_from(a, b, delta, sb)
    evac, _ = worst_case_time(robots, n_theta=4001)
    if not np.isfinite(evac):
        return 1000.0
    return evac


def de_search(sb: int = +1, *, maxiter=400, popsize=40):
    bounds = [
        (0.0, TWO_PI),     # a
        (0.0, TWO_PI),     # b
        (-TWO_PI/3, TWO_PI/3),   # delta
    ]
    return differential_evolution(
        objective, bounds=bounds, args=(sb,),
        maxiter=maxiter, popsize=popsize, seed=0, tol=1e-8,
        polish=True, init="sobol", updating="deferred", workers=-1,
    )


def local_search(a0, b0, delta0, sb=+1, *, maxiter=20000):
    x0 = np.array([a0, b0, delta0])
    return minimize(
        objective, x0, args=(sb,), method="Nelder-Mead",
        options={"xatol": 1e-9, "fatol": 1e-10, "maxiter": maxiter, "adaptive": True},
    )


def report(res, sb):
    a, b, delta = res.x
    robots = robots_from(a, b, delta, sb)
    evac_fine, theta_fine = worst_case_time(robots, n_theta=40001)
    gap = coverage_gap_exact(a, b, delta, sb)
    print(f"  sb = {sb:+d}:  a = {a:.6f}  b = {b:.6f}  delta = {delta:.6f}")
    print(f"    coverage gap = {gap:.3e}")
    print(f"    objective    = {res.fun:.10f}")
    print(f"    fine-grid WC = {evac_fine:.10f}  (theta = {theta_fine:.6f})")
    delta_vs_A3 = EVAC_A3 - evac_fine
    marker = "better" if delta_vs_A3 > 1e-6 else ("tie" if abs(delta_vs_A3) < 1e-6 else "worse")
    print(f"    vs A_3       = {delta_vs_A3:+.10f}  ({marker})")


if __name__ == "__main__":
    print("=== DE over Family D (3-fold symmetric, two phases) ===")
    for sb in (+1, -1):
        print(f"\n[sb = {sb:+d}]")
        res = de_search(sb=sb, maxiter=500, popsize=40)
        report(res, sb)

    # Also try local searches from a few natural seeds:
    #   - "no second scan": a = 2*pi/3, b = 0, delta = 0  (degenerate A_3 symmetric)
    #   - "half-and-half":  a = pi/3, b = pi/3, delta = 0
    print()
    print("=== Local searches from natural seeds ===")
    seeds = [
        ("b=0 (symmetric A_3)", (2*np.pi/3, 0.0, 0.0)),
        ("a=b=pi/3 split",      (np.pi/3,   np.pi/3, 0.0)),
        ("a=b=pi/3 staggered",  (np.pi/3,   np.pi/3, np.pi/6)),
        ("a=1.21, b=0.88",      (1.2159,    0.8817, 0.0)),
    ]
    for name, (a0, b0, d0) in seeds:
        for sb in (+1, -1):
            print(f"\n[{name}, sb={sb:+d}]")
            res = local_search(a0, b0, d0, sb=sb)
            report(res, sb)
