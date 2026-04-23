"""Refined lower bound on 3-robot wireless evacuation time.

Central inequality (see notes/lb-refinement.md for derivation). For any
3-robot wireless algorithm, any x > 0, and any boundary exit p that is
unexplored at time 1+x:

    evac(p) >= (1+x) + max_j || R_j(1+x) - p ||,

where R_j(1+x) is robot j's position at time 1+x.

Strategy to get a universal LB:
    For each x, over all algorithm states (R_1, R_2, R_3, E) consistent
    with the algorithm having traveled <= 1+x units,
        LB(x) = inf_alg  sup_{p in E^c}  (1+x) + max_j || R_j - p ||.
    Then LB = sup_x LB(x).

Compare to CGK's (1+x) + 2 sin(xk/2), attained by the "two-point"
adversary where the bound is the chord across the unexplored arc.

This script explores algorithmic families in a 3-fold symmetric
parameterization and computes numerical LBs.

Takeaway we want: for some x range, the refined LB beats the CGK
closed-form bound of 4.15937 at its optimum x = (2/3) arccos(-1/3).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.optimize import minimize, minimize_scalar

TWO_PI = 2.0 * math.pi
K = 3


def chord_between_angles(alpha: float, beta: float) -> float:
    """Chord length between boundary points e(alpha), e(beta) on unit circle."""
    d = (beta - alpha) % TWO_PI
    if d > math.pi:
        d = TWO_PI - d
    return 2.0 * math.sin(d / 2.0)


def chord_point_to_boundary(point: np.ndarray, beta: float) -> float:
    """Chord from arbitrary 2D point to e(beta) on unit circle."""
    bx, by = math.cos(beta), math.sin(beta)
    return math.hypot(point[0] - bx, point[1] - by)


# --------------------------------------------------------------------------
# CGK's bound for reference
# --------------------------------------------------------------------------

def cgk_bound(x: float, k: int = K) -> float:
    """CGK+14 Lemma 6: for pi/k <= x <= 2 pi / k, evac >= 1 + x + 2 sin(xk/2)."""
    return 1.0 + x + 2.0 * math.sin(x * k / 2.0)


def cgk_optimum(k: int = K) -> tuple[float, float]:
    """Return (x_opt, LB_opt) from CGK+14 Thm 7."""
    # Stationarity: 1 + k cos(xk/2) = 0, so cos(xk/2) = -1/k.
    x_opt = (2.0 / k) * math.acos(-1.0 / k)
    lb = cgk_bound(x_opt, k)
    return x_opt, lb


# --------------------------------------------------------------------------
# Refined LB for algorithm family A:
#
#   3-fold symmetric, each robot j goes to angle phi_j := phi_0 + 2 pi (j-1)/3,
#   arriving at time 1; then scans CCW by arclength L, so scan arc
#   = [phi_j, phi_j + L], ending at R_j = e(phi_j + L) at time 1 + L;
#   then waits there.
#
# At any time t >= 1 + L, state is:
#   R_j = e(phi_j + L), E = union of three scan arcs [phi_j, phi_j + L].
# Unexplored = three arcs of length 2 pi/3 - L each.
# --------------------------------------------------------------------------

def family_A_lb(x: float, L: float, phi0: float = 0.0) -> float:
    """
    Refined LB for family A at time 1+x with scan length L.

    For valid state we need L <= x (so robot has completed scan by 1+x)
    and L < 2 pi / 3 (so unexplored is nonempty).
    """
    if L >= TWO_PI / 3.0 - 1e-12:
        return 1.0 + x  # vacuous; no unexplored
    robot_angles = [phi0 + 2.0 * math.pi * j / 3.0 + L for j in range(3)]
    # Unexplored arcs:
    unexplored_starts = [phi0 + 2.0 * math.pi * j / 3.0 + L for j in range(3)]
    unexplored_ends = [phi0 + 2.0 * math.pi * (j + 1) / 3.0 for j in range(3)]
    # Sample boundary points on each unexplored arc, compute max_j chord,
    # and find sup over the arc.
    best = 0.0
    best_beta = None
    for (a, b) in zip(unexplored_starts, unexplored_ends):
        for frac in np.linspace(0.0, 1.0, 401):
            beta = a + (b - a) * frac
            mx = max(chord_between_angles(ra, beta) for ra in robot_angles)
            if mx > best:
                best, best_beta = mx, beta
    return 1.0 + x + best


# --------------------------------------------------------------------------
# Closed-form max-chord for family A: by symmetry and the computations in
# notes/lb-refinement.md, the adversary maximum over an unexplored arc
# midpoint gives max angular distance pi - L/2 so chord = 2 cos(L/4).
# But we also need to check the arc ENDPOINTS, which are the ends of the
# scan arcs where a robot sits. At those endpoints, one robot is at
# chord 0; the OPPOSITE robot is at angular distance 4 pi/3 - L, chord
# 2 sin((4 pi/3 - L)/2) (as long as < pi). So the sup is
# max(2 cos(L/4), 2 sin((4 pi/3 - L)/2)).
# --------------------------------------------------------------------------

def family_A_closed_form(L: float) -> float:
    """Closed-form max chord for family A as a function of scan length L."""
    # Midpoint-of-unexplored-arc chord (interior critical):
    c_mid = 2.0 * math.cos(L / 4.0)
    # Endpoint chord (opposite-robot from a scan endpoint).
    # At beta = phi_1 + L (right endpoint of scan 1 = left endpoint of unexplored arc 1):
    # distance to R_3 = e(phi_3 + L) = e(phi_1 + 4 pi/3 + L) is 4 pi/3 (mod 2 pi).
    # Actually no, let me reconsider: at beta = phi_1 + 2 pi/3 - epsilon (right
    # edge of unexplored 1, = phi_2 - epsilon just before scan 2), distance to
    # R_1 = phi_1 + L is 2 pi/3 - L. Distance to R_2 = phi_2 + L - (phi_2 - eps)
    # = L + eps ~ L. Distance to R_3 is 2 pi/3 + L (from the other side, minimum
    # is 4 pi/3 - L, hmm depending on L).
    # Let's just compute numerically for the endpoint.
    beta = TWO_PI / 3.0 - 1e-9  # right edge of unexplored arc 1, assuming phi0 = 0
    robot_angles = [2.0 * math.pi * j / 3.0 + L for j in range(3)]
    c_edge = max(chord_between_angles(ra, beta) for ra in robot_angles)
    return max(c_mid, c_edge)


# --------------------------------------------------------------------------
# Scan over x for family A: for each x, find best L (algorithm minimizes)
# subject to L <= x and L < 2 pi / 3. Report LB(x).
# --------------------------------------------------------------------------

def family_A_min_lb_at_x(x: float, n_L: int = 100) -> tuple[float, float, float]:
    """Return (best_L, best_max_chord, LB = 1 + x + max_chord) for family A."""
    L_max = min(x, TWO_PI / 3.0 - 1e-6)
    L_min = 1e-6
    Ls = np.linspace(L_min, L_max, n_L)
    best_chord = math.inf
    best_L = None
    for L in Ls:
        c = family_A_closed_form(L)
        if c < best_chord:
            best_chord, best_L = c, L
    # Refine with scipy.
    res = minimize_scalar(
        family_A_closed_form, bounds=(L_min, L_max), method="bounded",
        options={"xatol": 1e-10},
    )
    if res.fun < best_chord:
        best_chord, best_L = res.fun, res.x
    return best_L, best_chord, 1.0 + x + best_chord


# --------------------------------------------------------------------------
# Main driver
# --------------------------------------------------------------------------

def main() -> None:
    x_cgk, lb_cgk = cgk_optimum()
    print(f"CGK+14 optimum: x = {x_cgk:.6f}, LB = {lb_cgk:.6f}")
    print()

    print("Family A (3-fold symmetric, equal scans of length L, robots at scan ends):")
    print("  x       | best L    | max chord | refined LB  | CGK LB     | diff")
    print("  --------+-----------+-----------+-------------+------------+---------")
    best_refined = 0.0
    best_x = None
    best_L = None
    for x in np.linspace(0.6, 2.0, 29):
        L_opt, chord, lb = family_A_min_lb_at_x(x)
        cgk = cgk_bound(x) if math.pi / K <= x <= TWO_PI / K else float("nan")
        diff = lb - cgk if not math.isnan(cgk) else float("nan")
        print(f"  {x:.3f}   | {L_opt:.5f}   | {chord:.5f}   | {lb:.7f}    | {cgk:.7f}   | {diff:+.4f}")
        if lb > best_refined:
            best_refined, best_x, best_L = lb, x, L_opt

    print()
    print(f"Best refined LB (over x, family A): {best_refined:.7f} at x = {best_x:.4f}, L = {best_L:.4f}")
    print(f"Improvement over CGK:               {best_refined - lb_cgk:+.5f}")
    print()
    print("NOTE: family A is one specific algorithm class. The TRUE LB is the")
    print("      infimum of this quantity over ALL algorithms (all valid")
    print("      trajectories). For a valid proof we need to either widen the")
    print("      family or argue that family A is close to worst case.")


if __name__ == "__main__":
    main()
