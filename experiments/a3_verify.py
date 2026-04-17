"""Verify the simulator reproduces Czyzowicz et al.'s A3 upper bound at k=3.

A3 setup (Theorem 6):
    y = 4 pi/9 + 2 sqrt(3)/3 - 401/300                           # ~ 1.2143
    A at angle 0; B at angle 2 pi - y (i.e. y CW of A).
    r_1 -> A, scans CCW with length L = pi - y/2                 # ~ 2.5344
    r_2 -> B, scans CW  with length L = pi - y/2
    r_3 -> B, scans CCW for time y, returns to origin (time 1),
           then walks out to boundary angle pi - y/2 (time 1); parks.

The three scan arcs have lengths (L, L, y) with L + L + y = 2 pi: tight
coverage, no overlap.

Worst case (from hand-analysis):
    exit theta = pi - y/2: found simultaneously by r_1 and r_2 at
    t_found = 1 + L = 1 + pi - y/2. At that moment r_3 has been
    transiting out for time t_found - (2 + y) = pi - 3y/2 - 1 along
    its new radial direction, which points exactly at theta. Its chord
    to the exit is 1 - (pi - 3y/2 - 1) = 2 - pi + 3y/2.
    Evac time = (1 + pi - y/2) + (2 - pi + 3y/2) = 3 + y.

The paper states UB ~ 4 pi/9 + (2 sqrt 3 + 5)/3 + 1/600 ~ 4.2193, which
is 3 + y + 1/200 = 3 + y + 0.005 — an extra 0.005 of proof slop on top of
the actual worst case, so the real simulator value should be ~ 3 + y =
4.2143, not 4.2193.
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


def a3_robots():
    y = 4.0 * np.pi / 9.0 + 2.0 * np.sqrt(3.0) / 3.0 - 401.0 / 300.0
    L = np.pi - y / 2.0
    aA = 0.0
    aB = TWO_PI - y
    a_redeploy = np.pi - y / 2.0

    r1 = _robot_phases_simple(a=aA, L=L, s=+1)
    r2 = _robot_phases_simple(a=aB, L=L, s=-1)
    r3 = _robot_phases_two(a=aB, La=y, sa=+1, b=a_redeploy, Lb=0.0, sb=+1)
    return [r1, r2, r3], y, L, a_redeploy


if __name__ == "__main__":
    robots, y, L, a_redeploy = a3_robots()
    print("A3 parameters:")
    print(f"  y          = {y:.10f}")
    print(f"  L = pi-y/2 = {L:.10f}")
    print(f"  redeploy   = {a_redeploy:.10f} rad = pi - y/2")
    print()
    print(f"Coverage OK?  {_coverage_ok(robots, n_grid=8001)}")

    # Fine-grid worst-case sweep.
    for n_theta in (2001, 8001, 40001, 200001):
        evac, theta = worst_case_time(robots, n_theta=n_theta)
        print(f"  n_theta={n_theta:>7d}:  worst={evac:.10f}  theta={theta:.6f}")

    print()
    print("Reference values:")
    print(f"  3 + y                            = {3.0 + y:.10f}")
    print(f"  Paper UB 4pi/9+(2sqrt3+5)/3+1/600 = "
          f"{4*np.pi/9 + (2*np.sqrt(3)+5)/3 + 1/600:.10f}")
    print(f"  Naive baseline 1 + 2pi/3 + sqrt 3 = {1.0 + 2*np.pi/3 + np.sqrt(3):.10f}")
