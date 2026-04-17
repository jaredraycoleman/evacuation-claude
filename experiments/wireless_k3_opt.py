"""Numerical search for an improved wireless-k=3 upper bound.

Trajectory family (per robot i in {0,1,2}):
    - wait w_i at origin, then
    - walk radially to boundary angle a_i (arrives at time w_i + 1), then
    - scan CCW/CW (direction s_i in {+1, -1}) covering arc of length L_i,
      so its boundary sweep over t in [w_i+1, w_i+1+L_i] visits the arc
      A_i := {a_i + s_i * u : u in [0, L_i]} (mod 2*pi).

Coverage constraint: union_i A_i = [0, 2*pi).

Hit time (robot i reaches exit angle theta):
    delta = ((theta - a_i) * s_i) mod 2*pi    # signed arc distance from a_i along s_i
    hit_time_i(theta) = w_i + 1 + delta  if delta <= L_i  else  +inf

Once t_found = min_i hit_time_i(theta), every robot walks straight to the
exit; evac time = t_found + max_i chord_i where chord_i is the Euclidean
distance from robot i's position at t_found to the exit point.

Gauge-fix: a_0 = 0, s_0 = +1.
Free continuous: a_1, a_2, L_0, L_1, L_2, w_0, w_1, w_2  (8 params).
Discrete: (s_1, s_2) in {+1,-1}^2  (4 cases).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple

import numpy as np
from scipy.optimize import differential_evolution

TWO_PI = 2.0 * np.pi


class Params(NamedTuple):
    a: np.ndarray     # shape (3,), boundary-arrival angles
    L: np.ndarray     # shape (3,), scan arc lengths
    w: np.ndarray     # shape (3,), wait times at origin
    s: np.ndarray     # shape (3,), scan directions in {+1, -1}


def _positions(p: Params, t: np.ndarray) -> np.ndarray:
    """Positions of all 3 robots at each time in t.

    Parameters
    ----------
    t : shape (N,)
    Returns: shape (3, N, 2)
    """
    t = np.asarray(t, dtype=float)
    out = np.zeros((3, t.size, 2))
    for i in range(3):
        w_i, a_i, L_i, s_i = p.w[i], p.a[i], p.L[i], p.s[i]
        dir_xy = np.array([np.cos(a_i), np.sin(a_i)])
        # phase 0: at origin  (t <= w_i)
        # phase 1: approach   (w_i < t <= w_i + 1)
        # phase 2: scan       (w_i+1 < t <= w_i+1+L_i)
        # phase 3: parked     (t > w_i+1+L_i)
        phase1 = (t > w_i) & (t <= w_i + 1)
        phase2 = (t > w_i + 1) & (t <= w_i + 1 + L_i)
        phase3 = t > w_i + 1 + L_i

        # phase 1 positions
        r1 = (t[phase1] - w_i)[:, None] * dir_xy[None, :]
        out[i, phase1] = r1

        # phase 2 positions
        u = t[phase2] - w_i - 1.0
        ang2 = a_i + s_i * u
        out[i, phase2, 0] = np.cos(ang2)
        out[i, phase2, 1] = np.sin(ang2)

        # phase 3 position (parked at scan endpoint)
        end_ang = a_i + s_i * L_i
        out[i, phase3, 0] = np.cos(end_ang)
        out[i, phase3, 1] = np.sin(end_ang)
    return out


def _hit_times(p: Params, thetas: np.ndarray) -> np.ndarray:
    """Hit times for every (robot, theta). Shape (3, N_theta). +inf = no hit."""
    out = np.full((3, thetas.size), np.inf)
    for i in range(3):
        delta = ((thetas - p.a[i]) * p.s[i]) % TWO_PI   # arc from a_i along s_i
        covered = delta <= p.L[i] + 1e-12
        out[i, covered] = p.w[i] + 1.0 + delta[covered]
    return out


def _coverage_gap(p: Params, n_grid: int = 4001) -> float:
    """0 if union of scan arcs covers [0, 2*pi); else max uncovered arc length."""
    g = np.linspace(0.0, TWO_PI, n_grid, endpoint=False)
    ht = _hit_times(p, g)
    finite_any = np.isfinite(ht).any(axis=0)
    if finite_any.all():
        return 0.0
    # rough penalty proxy: fraction uncovered
    return (~finite_any).mean() * TWO_PI


def worst_case_time(p: Params, n_theta: int = 2001) -> tuple[float, float]:
    """Worst-case evac time (over theta) and the worst theta.

    Returns (+inf, nan) if coverage incomplete.
    """
    thetas = np.linspace(0.0, TWO_PI, n_theta, endpoint=False)
    ht = _hit_times(p, thetas)
    finite_any = np.isfinite(ht).any(axis=0)
    if not finite_any.all():
        return float("inf"), float("nan")

    t_found = ht.min(axis=0)                 # shape (N,)
    pos = _positions(p, t_found)             # shape (3, N, 2)
    exit_xy = np.stack([np.cos(thetas), np.sin(thetas)], axis=-1)  # (N, 2)
    diff = pos - exit_xy[None, :, :]
    chords = np.linalg.norm(diff, axis=-1)   # (3, N)
    evac = t_found + chords.max(axis=0)      # (N,)

    idx = int(np.argmax(evac))
    return float(evac[idx]), float(thetas[idx])


def _pack(a1, a2, L0, L1, L2, w0, w1, w2, s1, s2) -> Params:
    return Params(
        a=np.array([0.0, a1, a2]),
        L=np.array([L0, L1, L2]),
        w=np.array([w0, w1, w2]),
        s=np.array([1.0, s1, s2]),
    )


def _objective(x: np.ndarray, s1: int, s2: int) -> float:
    a1, a2, L0, L1, L2, w0, w1, w2 = x
    p = _pack(a1, a2, L0, L1, L2, w0, w1, w2, s1, s2)
    gap = _coverage_gap(p, n_grid=2001)
    if gap > 0.0:
        return 1000.0 + 50.0 * gap       # infeasible: huge penalty
    evac, _ = worst_case_time(p, n_theta=2001)
    return evac


@dataclass
class SearchResult:
    params: Params
    worst_time: float
    worst_theta: float
    discrete: tuple[int, int]


def search(
    *,
    maxiter: int = 400,
    popsize: int = 30,
    seed: int = 0,
    tol: float = 1e-7,
) -> SearchResult:
    bounds = [
        (0.0, TWO_PI),          # a1
        (0.0, TWO_PI),          # a2
        (0.1, TWO_PI),          # L0
        (0.1, TWO_PI),          # L1
        (0.1, TWO_PI),          # L2
        (0.0, 2.0),             # w0
        (0.0, 2.0),             # w1
        (0.0, 2.0),             # w2
    ]
    best: SearchResult | None = None
    for s1 in (+1, -1):
        for s2 in (+1, -1):
            res = differential_evolution(
                _objective,
                bounds=bounds,
                args=(s1, s2),
                maxiter=maxiter,
                popsize=popsize,
                seed=seed,
                tol=tol,
                polish=True,
                init="sobol",
                updating="deferred",
                workers=-1,
            )
            a1, a2, L0, L1, L2, w0, w1, w2 = res.x
            p = _pack(a1, a2, L0, L1, L2, w0, w1, w2, s1, s2)
            evac, theta = worst_case_time(p, n_theta=4001)
            cand = SearchResult(p, evac, theta, (s1, s2))
            print(f"  s=({s1:+d},{s2:+d}): worst={evac:.6f}  theta={theta:.4f}")
            if best is None or cand.worst_time < best.worst_time:
                best = cand
    assert best is not None
    return best


if __name__ == "__main__":
    print("Searching for improved k=3 wireless upper bound...")
    best = search(maxiter=400, popsize=30, seed=0)
    p = best.params
    print()
    print(f"Best worst-case evacuation time: {best.worst_time:.6f}")
    print(f"Worst exit angle theta = {best.worst_theta:.6f} rad")
    print(f"Discrete (s1, s2) = {best.discrete}")
    print(f"Naive baseline                 : {1.0 + 2.0*np.pi/3 + 2*np.sin(np.pi/3):.6f}")
    print(f"LB (Czyzowicz et al.)          : {3.0 + np.pi/3:.6f}")
    print()
    print("Parameters:")
    for i in range(3):
        print(f"  robot {i}: a={p.a[i]:.6f}  L={p.L[i]:.6f}  w={p.w[i]:.6f}  s={int(p.s[i]):+d}")
