"""Family C: A_3 extended with r_3 scanning after redeployment.

In A_3 (paper) and our balanced variant, r_3 has a two-phase trajectory
    approach -> scan of length y -> return-to-origin -> approach -> park
with no scan at the redeployment end. Here we release that "park" into a
genuine scan phase of length L_b in some direction s_b, and re-balance.

Free parameters:
    y         deployment offset (B angle = 2 pi - y)
    La        r_3 first-scan length (no longer tied to y)
    b         r_3 redeployment angle
    Lb        r_3 second-scan length (NEW; was 0 in A_3)
    sb        r_3 second-scan direction in {+1, -1}  (discrete)
    dL        L_1 - L_2 asymmetry  (was 0 at A_3 optimum)
    shift     r_1/r_2 arc shift; fixed at 0 for now (gauge)

Coverage constraint: union of four arcs (r_1, r_2, r_3 first, r_3 second)
covers [0, 2 pi).  Here we enforce coverage numerically rather than
algebraically (the arc structure is no longer a tidy partition once r_3
has two disjoint scans).

Objective: worst-case evacuation time over exit angle theta.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import differential_evolution, minimize

from experiments.wireless_k3_opt2 import (
    _robot_phases_simple,
    _robot_phases_two,
    worst_case_time,
    _hit_times_phases,
)

TWO_PI = 2.0 * np.pi
Y_OPT_A3 = 1.2158578321292429   # balanced A_3 y_opt
EVAC_A3  = 4.2185169939          # balanced A_3 worst-case evac


def robots_from(y, La, L1, L2, b, Lb, sb):
    """Build the 3-robot trajectory set. r_3 scans La then redeploys and
    scans Lb in direction sb. Gauge: r_1 at angle 0 scanning CCW."""
    aB = TWO_PI - y
    r1 = _robot_phases_simple(a=0.0, L=L1, s=+1)
    r2 = _robot_phases_simple(a=aB, L=L2, s=-1)
    r3 = _robot_phases_two(a=aB, La=La, sa=+1, b=b, Lb=Lb, sb=sb)
    return [r1, r2, r3]


def _arc_union_gap(arcs):
    """Exact coverage gap for a union of CCW arcs on the circle.

    `arcs` is a list of (start, length) pairs. An arc is the CCW sweep
    from angle `start` for `length` radians (mod 2*pi). Returns the total
    uncovered measure on [0, 2*pi).
    """
    # Normalize each arc to (s, e) with s in [0, 2*pi), e = s + length,
    # and possibly split into two pieces if it wraps past 2*pi.
    segments = []
    for start, length in arcs:
        if length <= 0:
            continue
        s = start % TWO_PI
        e = s + min(length, TWO_PI)
        if e <= TWO_PI:
            segments.append((s, e))
        else:
            segments.append((s, TWO_PI))
            segments.append((0.0, e - TWO_PI))
    if not segments:
        return TWO_PI

    # Merge overlapping / touching segments.
    segments.sort()
    merged = [segments[0]]
    for s, e in segments[1:]:
        last_s, last_e = merged[-1]
        if s <= last_e + 1e-12:
            merged[-1] = (last_s, max(last_e, e))
        else:
            merged.append((s, e))

    covered = sum(e - s for s, e in merged)
    return max(0.0, TWO_PI - covered)


def coverage_gap_exact(y, La, L1, L2, b, Lb, sb):
    """Arc union gap for Family C."""
    aB = (TWO_PI - y) % TWO_PI
    arcs = [
        (0.0,               L1),              # r_1 scans [0, L_1] CCW
        (aB - L2, L2),                        # r_2 scans [aB-L2, aB] (CCW-equivalent)
        (aB,                La),              # r_3 first scan CCW by La
    ]
    if sb > 0:
        arcs.append((b, Lb))                  # r_3 second scan CCW by Lb
    else:
        arcs.append((b - Lb, Lb))             # r_3 second scan CW -> CCW-equivalent
    # Normalize starts mod 2*pi.
    arcs = [(s % TWO_PI, ell) for s, ell in arcs]
    return _arc_union_gap(arcs)


def objective(x, sb: int = +1):
    y, La, L1, L2, b, Lb = x
    if y <= 0 or La < 0 or L1 <= 0 or L2 <= 0 or Lb < 0:
        return 1000.0
    if Lb > TWO_PI or La > TWO_PI or L1 > TWO_PI or L2 > TWO_PI:
        return 1000.0
    gap = coverage_gap_exact(y, La, L1, L2, b, Lb, sb)
    if gap > 1e-9:
        return 500.0 + 100.0 * gap
    robots = robots_from(y, La, L1, L2, b, Lb, sb)
    evac, _ = worst_case_time(robots, n_theta=4001)
    if not np.isfinite(evac):
        return 1000.0
    return evac


def search_from_a3(sb: int = +1, *, maxiter: int = 2000):
    """Local search starting at the A_3 optimum plus a small nudge on Lb."""
    y0 = Y_OPT_A3
    L0 = np.pi - y0 / 2.0
    # Seed: A_3 optimum with a token Lb = 0.1 to escape the boundary.
    x0 = np.array([y0, y0, L0, L0, np.pi - y0 / 2.0, 0.1])
    res = minimize(
        objective, x0, args=(sb,), method="Nelder-Mead",
        options={"xatol": 1e-9, "fatol": 1e-10, "maxiter": maxiter, "adaptive": True},
    )
    return res


def de_search(sb: int = +1, *, maxiter: int = 300, popsize: int = 30):
    bounds = [
        (0.3, 1.8),     # y
        (0.0, TWO_PI),  # La
        (0.1, TWO_PI),  # L1
        (0.1, TWO_PI),  # L2
        (0.0, TWO_PI),  # b
        (0.0, TWO_PI),  # Lb
    ]
    return differential_evolution(
        objective, bounds=bounds, args=(sb,),
        maxiter=maxiter, popsize=popsize, seed=0, tol=1e-8,
        polish=True, init="sobol", updating="deferred", workers=-1,
    )


def report(res, sb):
    y, La, L1, L2, b, Lb = res.x
    robots = robots_from(y, La, L1, L2, b, Lb, sb)
    evac_fine, theta_fine = worst_case_time(robots, n_theta=40001)
    gap = coverage_gap_exact(y, La, L1, L2, b, Lb, sb)
    print(f"  sb = {sb:+d}")
    print(f"  y = {y:.8f}  La = {La:.8f}  Lb = {Lb:.8f}")
    print(f"  L1 = {L1:.8f}  L2 = {L2:.8f}  b = {b:.8f}")
    print(f"  coverage gap (exact) = {gap:.3e}")
    print(f"  objective value      = {res.fun:.10f}")
    print(f"  fine-grid worst-case = {evac_fine:.10f}  at theta = {theta_fine:.6f}")
    delta = EVAC_A3 - evac_fine
    print(f"  improvement over A_3 = {delta:+.10f}  ({'better' if delta > 0 else 'worse'})")


if __name__ == "__main__":
    print("=== Local search from A_3 seed with Lb active ===")
    for sb in (+1, -1):
        print(f"\n[sb = {sb:+d}]")
        res = search_from_a3(sb=sb, maxiter=20000)
        report(res, sb)

    print()
    print("=== Differential evolution (wider search) ===")
    for sb in (+1, -1):
        print(f"\n[sb = {sb:+d}]")
        res = de_search(sb=sb, maxiter=300, popsize=30)
        report(res, sb)
