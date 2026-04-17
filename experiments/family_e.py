"""Family E: r_1 simple, r_2 and r_3 both have two-phase trajectories.

In A_3 only r_3 does the return-to-origin + redeploy trick. Family D
imposes full 3-fold symmetry (all robots redeploy), which collapses to
the naive 4.826 bound. Family E sits in between: r_2 and r_3 redeploy,
r_1 simply scans.

Free continuous parameters (gauge: r_1 at angle 0, scan CCW):
    y            angular offset of r_2/r_3's first target
    L1           r_1 scan length
    a2, L2a, b2, L2b   r_2 first+second phase parameters
    a3, L3a, b3, L3b   r_3 first+second phase parameters
Discrete (scan directions, mostly ±): s_{2a}, s_{2b}, s_{3a}, s_{3b}.

Coverage constraint: union of up to 5 arcs (r_1 + r_2 x2 + r_3 x2) covers
[0, 2 pi).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import differential_evolution, minimize

from experiments.wireless_k3_opt2 import (
    _robot_phases_simple,
    _robot_phases_two,
    worst_case_time,
)
from experiments.family_c import _arc_union_gap

TWO_PI = 2.0 * np.pi
EVAC_A3 = 4.2185169939


def robots_from(L1, a2, L2a, b2, L2b, s2a, s2b, a3, L3a, b3, L3b, s3a, s3b):
    r1 = _robot_phases_simple(a=0.0, L=L1, s=+1)
    r2 = _robot_phases_two(a=a2, La=L2a, sa=s2a, b=b2, Lb=L2b, sb=s2b)
    r3 = _robot_phases_two(a=a3, La=L3a, sa=s3a, b=b3, Lb=L3b, sb=s3b)
    return [r1, r2, r3]


def coverage_gap(L1, a2, L2a, b2, L2b, s2a, s2b, a3, L3a, b3, L3b, s3a, s3b):
    arcs = [(0.0, L1)]
    for (a, L, s) in [(a2, L2a, s2a), (b2, L2b, s2b),
                      (a3, L3a, s3a), (b3, L3b, s3b)]:
        if s > 0:
            arcs.append((a % TWO_PI, L))
        else:
            arcs.append(((a - L) % TWO_PI, L))
    return _arc_union_gap(arcs)


def objective(x, s2a, s2b, s3a, s3b):
    L1, a2, L2a, b2, L2b, a3, L3a, b3, L3b = x
    if any(v < 0 for v in (L1, L2a, L2b, L3a, L3b)):
        return 1000.0
    if any(v > TWO_PI for v in (L1, L2a, L2b, L3a, L3b)):
        return 1000.0
    gap = coverage_gap(L1, a2, L2a, b2, L2b, s2a, s2b, a3, L3a, b3, L3b, s3a, s3b)
    if gap > 1e-9:
        return 500.0 + 100.0 * gap
    robots = robots_from(L1, a2, L2a, b2, L2b, s2a, s2b, a3, L3a, b3, L3b, s3a, s3b)
    evac, _ = worst_case_time(robots, n_theta=4001)
    if not np.isfinite(evac):
        return 1000.0
    return evac


def seed_as_A3_plus_r2_redeploy():
    """Seed Family E from the balanced A_3 optimum, with r_2 given a
    degenerate 2nd phase (L2b = 0 at an arbitrary redeploy angle)."""
    y = 1.2158578321
    L = np.pi - y / 2.0
    return {
        "L1":  L,
        "a2":  TWO_PI - y, "L2a": L,        "b2":  np.pi - y/2, "L2b": 0.0,
        "s2a": -1,         "s2b": +1,
        "a3":  TWO_PI - y, "L3a": y,        "b3":  np.pi - y/2, "L3b": 0.0,
        "s3a": +1,         "s3b": +1,
    }


def local_search(seed, maxiter=30000):
    s2a = seed["s2a"]; s2b = seed["s2b"]; s3a = seed["s3a"]; s3b = seed["s3b"]
    x0 = np.array([seed["L1"], seed["a2"], seed["L2a"], seed["b2"], seed["L2b"],
                   seed["a3"], seed["L3a"], seed["b3"], seed["L3b"]])
    return minimize(
        objective, x0, args=(s2a, s2b, s3a, s3b), method="Nelder-Mead",
        options={"xatol": 1e-9, "fatol": 1e-10, "maxiter": maxiter, "adaptive": True},
    ), (s2a, s2b, s3a, s3b)


def report(res, discrete):
    s2a, s2b, s3a, s3b = discrete
    L1, a2, L2a, b2, L2b, a3, L3a, b3, L3b = res.x
    robots = robots_from(L1, a2, L2a, b2, L2b, s2a, s2b,
                         a3, L3a, b3, L3b, s3a, s3b)
    gap = coverage_gap(L1, a2, L2a, b2, L2b, s2a, s2b,
                       a3, L3a, b3, L3b, s3a, s3b)
    evac_fine, theta_fine = worst_case_time(robots, n_theta=40001)
    print(f"  s = (s2a={s2a:+d}, s2b={s2b:+d}, s3a={s3a:+d}, s3b={s3b:+d})")
    print(f"  r_1: a=0  L={L1:.6f}  s=+1")
    print(f"  r_2: a={a2:.4f} La={L2a:.4f}  b={b2:.4f} Lb={L2b:.4f}")
    print(f"  r_3: a={a3:.4f} La={L3a:.4f}  b={b3:.4f} Lb={L3b:.4f}")
    print(f"  coverage gap = {gap:.3e}")
    print(f"  objective    = {res.fun:.10f}")
    print(f"  fine-grid WC = {evac_fine:.10f}  at theta={theta_fine:.6f}")
    delta = EVAC_A3 - evac_fine
    marker = "BETTER" if delta > 1e-6 else ("tie" if abs(delta) < 1e-6 else "worse")
    print(f"  vs A_3       = {delta:+.10f}  ({marker})")


if __name__ == "__main__":
    print("=== Family E: r_1 simple, r_2 and r_3 both two-phase ===")
    print()
    seed = seed_as_A3_plus_r2_redeploy()
    print("[Seed: A_3 optimum + r_2 with degenerate (L2b=0) redeploy]")
    res, discrete = local_search(seed)
    report(res, discrete)

    print()
    print("[Seed variations: L2b active to break the degeneracy]")
    for L2b0 in (0.05, 0.2, 0.5):
        seed["L2b"] = L2b0
        res, discrete = local_search(seed)
        print(f"\n-- L2b seed = {L2b0} --")
        report(res, discrete)
