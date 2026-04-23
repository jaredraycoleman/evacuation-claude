"""Evaluate the refined LB (evac >= (1+x) + max_j ||R_j(1+x) - p||)
on the SPECIFIC algorithm A_3(y_opt), for a range of x and adversary
exit angles.

This is a sanity check: for A_3 we KNOW evac = 4.21852 (UB = proved
worst-case). So the refined LB applied to A_3 should be <= 4.21852 for
every (x, p). The question is: is it > 4.15937 (CGK's LB) at some (x, p)?

If yes, then the refined LB inequality is non-trivially strong on A_3
itself. Whether it also gives a UNIVERSAL LB > 4.15937 is the next
question - that requires showing every algorithm state at its 1+x has
max_j ||R_j - p*|| sufficiently large for some unexplored p*.
"""
from __future__ import annotations

import math

import numpy as np

TWO_PI = 2.0 * math.pi

# Constants from paper / analysis/a3_balance.wls
Y_OPT = 1.2158578321292429094342098320811533017778
L_OPT = math.pi - Y_OPT / 2.0  # = pi - y_opt/2
T_STAR = 4.218516993897451331118770347651630972185


# --------------------------------------------------------------------------
# A_3(y) trajectories as a function of (y, t).
# --------------------------------------------------------------------------

def e(theta: float) -> np.ndarray:
    return np.array([math.cos(theta), math.sin(theta)])


def a3_position(j: int, t: float, y: float) -> np.ndarray:
    """Position of robot j in {1, 2, 3} at time t, for algorithm A_3(y)."""
    L = math.pi - y / 2.0
    if t <= 0:
        return np.array([0.0, 0.0])
    if t <= 1:
        # Radial approach 0 -> start point
        if j == 1:
            end = e(0.0)
        elif j == 2:
            end = e(TWO_PI - y)
        else:  # j == 3
            end = e(TWO_PI - y)
        return t * end

    elapsed = t - 1.0  # time spent past arrival at boundary
    if j == 1:
        # Scan CCW from e(0) on [1, 1+L]; park at e(L) after.
        if elapsed <= L:
            return e(elapsed)
        return e(L)
    if j == 2:
        # Scan CW from e(2 pi - y) on [1, 1+L]; park at e(pi - y/2) after.
        if elapsed <= L:
            return e((TWO_PI - y) - elapsed)
        return e(TWO_PI - y - L)  # = e(L)
    # j == 3
    #   [1, 1+y]: scan CCW from e(2 pi - y) to e(2 pi) = e(0)
    #   [1+y, 2+y]: straight line e(0) -> origin
    #   [2+y, 3+y]: straight line origin -> e(L)
    #   [3+y, inf): parked at e(L)
    if elapsed <= y:
        return e((TWO_PI - y) + elapsed)
    if elapsed <= 1 + y:
        # Boundary-to-origin along real axis.
        frac = (elapsed - y) / 1.0
        return e(0.0) * (1.0 - frac)
    if elapsed <= 2 + y:
        frac = elapsed - (1 + y)
        return frac * e(L)
    return e(L)


def a3_explored(t: float, y: float) -> list[tuple[float, float]]:
    """Return the set of arcs on [0, 2 pi) explored by time t under A_3(y).

    Each arc is (a, b) with a < b in [0, 2 pi), meaning e(theta) for theta in
    [a, b] has been visited by some robot. Arcs are normalized so a,b < 2 pi
    and may wrap (represented as (a, b+2pi) with b < a; we keep them simple).
    """
    L = math.pi - y / 2.0
    if t <= 1:
        return []
    elapsed = t - 1.0
    arcs: list[tuple[float, float]] = []
    # Robot 1: CCW from 0
    a1 = min(elapsed, L)
    if a1 > 0:
        arcs.append((0.0, a1))
    # Robot 2: CW from 2pi - y (scanning CW means angle decreases).
    a2 = min(elapsed, L)
    if a2 > 0:
        arcs.append((TWO_PI - y - a2, TWO_PI - y))
    # Robot 3: CCW from 2pi - y for up to time y.
    a3 = min(elapsed, y)
    if a3 > 0:
        # Scan from (2pi - y) to (2pi - y + a3). For a3 = y, endpoint = 2 pi.
        start = TWO_PI - y
        end = start + a3
        # Wrap: if end > 2 pi, split.
        if end <= TWO_PI:
            arcs.append((start, end))
        else:
            arcs.append((start, TWO_PI))
            arcs.append((0.0, end - TWO_PI))
    # Merge.
    return merge_arcs(arcs)


def merge_arcs(arcs: list[tuple[float, float]]) -> list[tuple[float, float]]:
    arcs = sorted(arcs)
    merged: list[tuple[float, float]] = []
    for a, b in arcs:
        if merged and a <= merged[-1][1] + 1e-12:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    return merged


def unexplored_arcs(t: float, y: float) -> list[tuple[float, float]]:
    explored = a3_explored(t, y)
    # Complement on [0, 2 pi)
    if not explored:
        return [(0.0, TWO_PI)]
    result: list[tuple[float, float]] = []
    prev = 0.0
    for a, b in explored:
        if a > prev + 1e-12:
            result.append((prev, a))
        prev = b
    if prev < TWO_PI - 1e-12:
        result.append((prev, TWO_PI))
    return result


# --------------------------------------------------------------------------
# Refined-LB evaluation
# --------------------------------------------------------------------------

def max_chord_to_robots(p: np.ndarray, robots: list[np.ndarray]) -> float:
    return max(float(np.linalg.norm(p - r)) for r in robots)


def refined_lb_on_a3(x: float, y: float, n_samples: int = 4001) -> tuple[float, float, float]:
    """
    Compute, for A_3(y) at time t = 1 + x:
        sup_{p in unexplored} [(1+x) + max_j ||R_j(1+x) - p||]
    where p = e(beta) for beta in the unexplored arcs.
    Returns (1+x, best_chord, refined_lb).
    """
    t = 1.0 + x
    robots = [a3_position(j, t, y) for j in (1, 2, 3)]
    unexp = unexplored_arcs(t, y)
    best_chord = 0.0
    best_beta = None
    for (a, b) in unexp:
        for frac in np.linspace(0.0, 1.0, n_samples):
            beta = a + (b - a) * frac
            p = e(beta)
            c = max_chord_to_robots(p, robots)
            if c > best_chord:
                best_chord, best_beta = c, beta
    return 1.0 + x, best_chord, 1.0 + x + best_chord


# --------------------------------------------------------------------------
# Sweep over x and report
# --------------------------------------------------------------------------

def main() -> None:
    y = Y_OPT

    # Range of x to consider. CGK's LB uses x in [pi/3, 2 pi/3] ~ [1.047, 2.094].
    # Our UB worst case is at x = theta_1 = 2 pi/3 - y/2 ~ 1.486.
    xs = np.linspace(0.7, 2.0, 27)

    print(f"A_3(y_opt) with y_opt = {y:.10f}, UB = {T_STAR:.10f}")
    print()
    print("  x       | 1+x     | max chord | refined LB  | CGK LB     | UB-refLB")
    print("  --------+---------+-----------+-------------+------------+---------")
    best = 0.0
    best_x = None
    for x in xs:
        one_plus_x, chord, lb = refined_lb_on_a3(x, y)
        cgk = 1 + x + 2 * math.sin(x * 3 / 2) if math.pi / 3 <= x <= TWO_PI / 3 else float("nan")
        slack = T_STAR - lb
        print(f"  {x:.3f}   | {one_plus_x:.4f}  | {chord:.5f}   | {lb:.7f}    | {cgk:10.7f}  | {slack:+.5f}")
        if lb > best:
            best, best_x = lb, x

    x_cgk_opt = (2.0 / 3.0) * math.acos(-1.0 / 3.0)
    lb_cgk_opt = 1 + x_cgk_opt + 2 * math.sin(x_cgk_opt * 3.0 / 2.0)
    print()
    print(f"Max refined LB over x (on A_3):     {best:.7f} at x = {best_x:.4f}")
    print(f"CGK optimum (universal LB):         {lb_cgk_opt:.7f}")
    print(f"Refined improvement (on A_3 only):  {best - lb_cgk_opt:+.5f}")
    print(f"Distance to UB (on A_3):            {T_STAR - best:.7f}")
    print()
    print("INTERPRETATION:")
    print("  Refined LB on A_3 at best x tells us how close the inequality")
    print("  (valid for any alg) comes to the ACTUAL UB for A_3. If this")
    print("  exceeds CGK at the same x, the inequality has bite on A_3.")
    print("  For a UNIVERSAL LB improvement we still need it to hold")
    print("  against OTHER algs (not just A_3).")


if __name__ == "__main__":
    main()
