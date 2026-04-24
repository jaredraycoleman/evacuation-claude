"""Numerical minimax for the uniform LB at x = x_CGK^*.

For any algorithm state at time 1 + x^* = 2.274, let
    e = |E(2.274)| in [0, 3x^*],
    R_j = robot positions in closed disk,
where the state is reachable by some valid trajectory.

Bounds:
    state-tight CGK(e) = 1 + x^* + 2 sin((2pi - e)/2)   for e in [pi, 3x^*],
                       = 1 + x^* + 2                     for e in [0, pi).
    refined(state)     = 1 + x^* + sup_{p in E^c} max_j ||R_j - p||.

combined(state) = max(state-tight CGK, refined).
universal LB = inf_{valid state} combined(state).

We explore the "2-coincident with balanced slack" parametrization:
    tau in [0, x^*]: slack budget per robot.
    R_1 = (1 - tau) e(d), R_2 = R_3 = (1 - tau) e(0).
    Scan arcs of length L = x^* - tau each.
    d = pi + L/2 (optimal for this subfamily).
    |E| = 3L = 3(x^* - tau) = 3x^* - 3 tau.

With delta := 3 tau, |E| = 3x^* - delta.

Also explore other parametrizations (3-fold symmetric with slack,
fully asymmetric slack) to check whether 2-coincident is the true
minimizer.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar

TWO_PI = 2.0 * math.pi
X_STAR = (2.0 / 3.0) * math.acos(-1.0 / 3.0)  # ~ 1.2737554908
ONE_PLUS_XSTAR = 1.0 + X_STAR
CGK_LB = 1.0 + X_STAR + 2.0 * math.sin(3.0 * X_STAR / 2.0)
SCAN_TIGHT_REFINED = 1.0 + X_STAR + 2.0 * math.cos(X_STAR / 4.0)


def state_tight_cgk(e: float) -> float:
    if e >= math.pi:
        return 1.0 + X_STAR + 2.0 * math.sin((TWO_PI - e) / 2.0)
    return 1.0 + X_STAR + 2.0


def chord_point_to_boundary(R: np.ndarray, beta: float) -> float:
    p = np.array([math.cos(beta), math.sin(beta)])
    return float(np.linalg.norm(R - p))


def sup_chord_balanced_2coincident(tau: float) -> float:
    """Sup chord for 2-coincident config with balanced inward by tau each.

    R_1 = (1 - tau) e(d), R_2 = R_3 = ((1 - tau), 0).
    Scan arcs of length L = x^* - tau each.
    d = pi + L/2.
    """
    r = 1.0 - tau
    L = X_STAR - tau
    d = math.pi + L / 2.0

    R1 = np.array([r * math.cos(d), r * math.sin(d)])
    R2 = np.array([r, 0.0])

    # Unexplored arcs: (L, d - L) and (d, 2*pi - L).
    betas = np.concatenate([
        np.linspace(L + 1e-9, d - L - 1e-9, 5000),
        np.linspace(d + 1e-9, TWO_PI - L - 1e-9, 5000),
    ])

    max_c = 0.0
    for beta in betas:
        p = np.array([math.cos(beta), math.sin(beta)])
        c = max(np.linalg.norm(R1 - p), np.linalg.norm(R2 - p))
        if c > max_c:
            max_c = c
    return max_c


def combined_balanced(tau: float) -> float:
    L = X_STAR - tau
    e = 3 * L  # |E| = 3 L for 3 disjoint arcs of length L each
    cgk = state_tight_cgk(e)
    sup_c = sup_chord_balanced_2coincident(tau)
    refined = 1.0 + X_STAR + sup_c
    return max(cgk, refined)


def main():
    print(f"x^* = {X_STAR:.10f}")
    print(f"CGK LB:           {CGK_LB:.8f}")
    print(f"Scan-tight refined (tau=0): {SCAN_TIGHT_REFINED:.8f}")
    print()
    print("Sweep of balanced 2-coincident config (tau = delta/3):")
    print("  tau     | delta = 3 tau | |E|     | CGK(e) | refined | combined")
    print("  --------+---------------+---------+--------+---------+----------")
    for tau in np.linspace(0.0, 0.05, 26):
        delta = 3 * tau
        e = 3 * X_STAR - delta
        cgk = state_tight_cgk(e)
        sc = sup_chord_balanced_2coincident(tau)
        refined = ONE_PLUS_XSTAR + sc
        combined = max(cgk, refined)
        print(f"  {tau:.5f} | {delta:.5f}       | {e:.4f}  | {cgk:.5f} | {refined:.5f} | {combined:.5f}")

    # Find minimum of combined over tau
    res = minimize_scalar(combined_balanced, bounds=(1e-6, 0.1), method="bounded",
                          options={"xatol": 1e-9})
    print()
    print(f"Min combined over tau (balanced 2-coincident): {res.fun:.8f}")
    print(f"At tau = {res.x:.6f}, delta = {3*res.x:.6f}")
    print(f"  = 4.15937 + {res.fun - CGK_LB:.6f}")
    print(f"  = 4.17306 - {SCAN_TIGHT_REFINED - res.fun:.6f}")

    # Also try: pure inward (no scan reduction). i.e., reduce L_j and move inward,
    # but WITHOUT balancing: one robot stays on boundary, others inward.
    # This shows whether asymmetric slack can be smaller.

    print()
    print("For reference, scan-tight gives combined = max(CGK(3x*), 4.17306) =",
          f"{max(CGK_LB, SCAN_TIGHT_REFINED):.6f}")


if __name__ == "__main__":
    main()
