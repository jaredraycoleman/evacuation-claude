"""Local search over the full A3 parameter family.

A3 has hidden rigid choices that the paper takes as given:
    - L_1 = L_2 = pi - y/2  (tight coverage, symmetric)
    - r_3 first-phase scan length La = y
    - r_3 redeploy angle b = pi - y/2 (antipodal to midpoint of AB)
    - r_3 second-phase scan length Lb = 0

We relax to let the optimizer pick these around the balanced y seed and see
if the 4.21852 bound tightens further.

Free params (all continuous):
    y            in [0.8, 1.6]   (B offset from A)
    L1, L2       in [y/2, pi]    (r_1, r_2 scan lengths; coverage constraint)
    La           in [0, pi]      (r_3 first-scan length)
    b_offset     in [-0.5, 0.5]  (redeploy angle = pi - y/2 + b_offset)

Coverage: arcs [0, L1] ∪ [2pi-y-L2, 2pi-y] ∪ [2pi-y, 2pi-y+La] must cover
[0, 2pi). We parametrize so equality holds at the tight end and add a
penalty for gaps.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from experiments.wireless_k3_opt2 import (
    _robot_phases_simple,
    _robot_phases_two,
    worst_case_time,
    _coverage_ok,
)

TWO_PI = 2.0 * np.pi


def robots_from(y, dL, b_offset):
    """A3 structure with hard tight coverage: L_a = y, L_1 + L_2 = 2 pi - y.
    Asymmetry parameter dL = L_1 - L_2.
    """
    L_mid = (TWO_PI - y) / 2.0     # = pi - y/2
    L1 = L_mid + dL / 2.0
    L2 = L_mid - dL / 2.0
    La = y
    aB = TWO_PI - y
    b = np.pi - y / 2.0 + b_offset
    r1 = _robot_phases_simple(a=0.0, L=L1, s=+1)
    r2 = _robot_phases_simple(a=aB, L=L2, s=-1)
    r3 = _robot_phases_two(a=aB, La=La, sa=+1, b=b, Lb=0.0, sb=+1)
    return [r1, r2, r3]


def objective(x):
    y, dL, b_offset = x
    if y <= 0:
        return 1000.0
    L_mid = (TWO_PI - y) / 2.0
    if abs(dL) > 2 * L_mid - 0.01:
        return 1000.0
    robots = robots_from(y, dL, b_offset)
    evac, _ = worst_case_time(robots, n_theta=16001)
    if not np.isfinite(evac):
        return 1000.0
    return evac


if __name__ == "__main__":
    # Seed at the balanced-y A3.
    y_seed = 1.2158578321
    x0 = np.array([y_seed, 0.0, 0.0])

    print(f"Seed (balanced y): y={y_seed:.6f}, dL=0, b_offset=0")
    print(f"Seed objective: {objective(x0):.10f}")
    print()

    res = minimize(
        objective, x0, method="Nelder-Mead",
        options={"xatol": 1e-9, "fatol": 1e-10, "maxiter": 20000, "adaptive": True},
    )
    y_o, dL_o, b_o = res.x
    print("After local search:")
    print(f"  y          = {y_o:.10f}")
    print(f"  dL         = {dL_o:.10f}     (L1 = pi-y/2 + dL/2, L2 = pi-y/2 - dL/2)")
    print(f"  b_offset   = {b_o:.10f}     (redeploy = pi - y/2 + b_offset)")
    print(f"  opt evac   = {res.fun:.10f}")

    robots = robots_from(*res.x)
    evac_fine, theta_fine = worst_case_time(robots, n_theta=200001)
    cov = _coverage_ok(robots, n_grid=8001)
    print(f"  fine-grid  = {evac_fine:.10f}  (theta = {theta_fine:.6f})")
    print(f"  coverage OK = {cov}")

    print()
    paper = 4*np.pi/9 + (2*np.sqrt(3)+5)/3 + 1/600
    balanced_y = 1.0 + 2.0*np.pi/3.0 - y_seed/2.0 + np.sqrt(3.0)
    print(f"  paper UB       = {paper:.10f}")
    print(f"  balanced-y UB  = {balanced_y:.10f}  (current best 1-parameter)")
    print(f"  local-search   = {evac_fine:.10f}")
    print(f"  improvement vs paper  = {paper - evac_fine:.10f}")
    print(f"  improvement vs y-only = {balanced_y - evac_fine:.10f}")
