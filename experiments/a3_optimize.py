"""Optimize the A3 deployment parameter y to balance the two worst-case regimes.

Hand-analysis of Czyzowicz et al.'s A3 (three robots, r_3 with return-to-
center + redeploy) identifies two dominant worst-case angles:

  Case 1 (exit in r_1's arc; bottleneck r_2 diametrically opposite):
    theta* = 2 pi/3 - y/2
    evac_1(y) = 1 + 2 pi/3 - y/2 + sqrt 3            # decreasing in y

  Case 2 (exit at antipodal meeting point of r_1, r_2 arcs; bottleneck r_3
          still transiting outward along the exit direction):
    theta** = pi - y/2
    evac_2(y) = 3 + y                                # increasing in y

Balance: evac_1 = evac_2  =>  y* = 4 pi/9 + 2 sqrt 3/3 - 4/3.
The paper chooses y_paper = 4 pi/9 + 2 sqrt 3/3 - 401/300 = y* - 1/300.

Expected improvement at y*: UB = 3 + y* = 5/3 + 4 pi/9 + 2 sqrt 3/3,
which is 1/600 smaller than the paper's UB of 4 pi/9 + (2 sqrt 3 + 5)/3
+ 1/600.

This script plugs y* into our simulator and confirms the new worst case.
It also sweeps y near y* to show the tradeoff curve.
"""
from __future__ import annotations

import numpy as np

from experiments.wireless_k3_opt2 import (
    _robot_phases_simple,
    _robot_phases_two,
    worst_case_time,
    _coverage_ok,
)

TWO_PI = 2.0 * np.pi


def a3_robots_at(y: float):
    L = np.pi - y / 2.0
    aB = TWO_PI - y
    a_redeploy = np.pi - y / 2.0
    r1 = _robot_phases_simple(a=0.0, L=L, s=+1)
    r2 = _robot_phases_simple(a=aB, L=L, s=-1)
    r3 = _robot_phases_two(a=aB, La=y, sa=+1, b=a_redeploy, Lb=0.0, sb=+1)
    return [r1, r2, r3]


def evac_case1(y: float) -> float:
    return 1.0 + 2.0 * np.pi / 3.0 - y / 2.0 + np.sqrt(3.0)


def evac_case2(y: float) -> float:
    return 3.0 + y


if __name__ == "__main__":
    y_paper = 4.0 * np.pi / 9.0 + 2.0 * np.sqrt(3.0) / 3.0 - 401.0 / 300.0
    y_star = 4.0 * np.pi / 9.0 + 2.0 * np.sqrt(3.0) / 3.0 - 4.0 / 3.0

    print(f"y_paper      = {y_paper:.10f}")
    print(f"y*  (balance) = {y_star:.10f}")
    print(f"y_paper - y*  = {y_paper - y_star:.10f}   (expected -1/300 = {-1/300:.10f})")
    print()

    for label, y in [("paper", y_paper), ("y*", y_star)]:
        robots = a3_robots_at(y)
        ok = _coverage_ok(robots, n_grid=8001)
        evac_sim, theta = worst_case_time(robots, n_theta=200001)
        print(f"[{label}]  y={y:.10f}  coverage_ok={ok}")
        print(f"    evac_case1({y:.4f}) = {evac_case1(y):.10f}")
        print(f"    evac_case2({y:.4f}) = {evac_case2(y):.10f}")
        print(f"    simulator max       = {evac_sim:.10f}  at theta={theta:.6f}")
        print()

    print("Paper UB  (4pi/9 + (2sqrt3+5)/3 + 1/600) = "
          f"{4*np.pi/9 + (2*np.sqrt(3)+5)/3 + 1/600:.10f}")
    print("Balanced  (3 + y*) = 5/3 + 4pi/9 + 2sqrt3/3 = "
          f"{5.0/3.0 + 4*np.pi/9 + 2*np.sqrt(3)/3:.10f}")
    print(f"Improvement:     = 1/600 = {1.0/600.0:.10f}")

    # Small sweep around y* to show the tradeoff curve.
    print()
    print("Worst-case evac as y varies (simulator):")
    for dy in np.linspace(-0.02, 0.02, 9):
        y = y_star + dy
        robots = a3_robots_at(y)
        evac_sim, _ = worst_case_time(robots, n_theta=40001)
        print(f"  y = y* + {dy:+.4f}  ->  evac = {evac_sim:.10f}  "
              f"(case1={evac_case1(y):.6f}, case2={evac_case2(y):.6f})")

    # Golden-section search for best y using the simulator (finds the
    # true 3-way balance).
    print()
    print("Scalar search for best y:")
    from scipy.optimize import minimize_scalar

    def obj(y):
        evac, _ = worst_case_time(a3_robots_at(float(y)), n_theta=40001)
        return evac

    res = minimize_scalar(obj, bracket=(y_paper - 0.05, y_paper, y_paper + 0.05),
                          method="brent", tol=1e-10)
    y_opt = float(res.x)
    evac_opt, theta_opt = worst_case_time(a3_robots_at(y_opt), n_theta=200001)
    print(f"  optimal y (scalar)     = {y_opt:.10f}")
    print(f"  optimal evac           = {evac_opt:.10f}")
    print(f"  worst theta            = {theta_opt:.6f}")
    print(f"  improvement over paper = {4*np.pi/9 + (2*np.sqrt(3)+5)/3 + 1/600 - evac_opt:.10f}")
