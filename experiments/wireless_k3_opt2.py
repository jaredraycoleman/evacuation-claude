"""Wireless k=3 UB search, Family B: one robot has a two-phase trajectory.

Family B (targeting A3's structural class):
    - Robots 0, 1: approach boundary angle a_i radially (time 1), then scan
      arc of length L_i in direction s_i in {+1, -1}.
    - Robot 2: approach a_2 radially, scan arc L_{2a} in direction s_{2a},
      then walk straight to the origin (time 1), then to boundary angle b_2
      (time 1), then scan arc L_{2b} in direction s_{2b}.

Gauge-fix: a_0 = 0, s_0 = +1.

Continuous free params (7): a_1, a_2, b_2, L_0, L_1, L_{2a}, L_{2b}.
Discrete (8 combos):        s_1, s_{2a}, s_{2b} in {+1, -1}.

Coverage constraint: union of all scan arcs covers [0, 2*pi).
Objective: worst-case evac time over exit angle theta.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple

import numpy as np
from scipy.optimize import differential_evolution

TWO_PI = 2.0 * np.pi


class Phases(NamedTuple):
    """Per-robot phases. Each phase has an approach angle, scan length,
    scan direction, and an 'approach start time' offset."""
    approach_start: np.ndarray  # (K,) times at which the radial approach begins
    a: np.ndarray               # (K,) approach angles
    L: np.ndarray               # (K,) scan lengths
    s: np.ndarray               # (K,) scan directions in {+1, -1}
    # Implicit:
    #   phase starts with linear walk of length 1 at unit speed,
    #   arriving on boundary at approach_start + 1; then scan for L time.


def _robot_phases_simple(a: float, L: float, s: float) -> Phases:
    return Phases(
        approach_start=np.array([0.0]),
        a=np.array([a]),
        L=np.array([L]),
        s=np.array([s]),
    )


def _robot_phases_two(
    a: float, La: float, sa: float,
    b: float, Lb: float, sb: float,
) -> Phases:
    """Robot does phase 1 (approach a, scan La in dir sa), then returns to
    origin via straight line (time 1), then phase 2 (approach b, scan Lb in
    dir sb). The walk from boundary endpoint to origin takes 1 unit of time
    (it's a radius).

    Phase 2's approach_start time = 1 + La + 1 = 2 + La.
    """
    return Phases(
        approach_start=np.array([0.0, 2.0 + La]),
        a=np.array([a, b]),
        L=np.array([La, Lb]),
        s=np.array([sa, sb]),
    )


def _positions_phases(ph: Phases, t: np.ndarray) -> np.ndarray:
    """Positions at times `t` for a robot with the given phase sequence.

    Between phases, the robot either scans the boundary (inside a phase)
    or walks radially in/out (between phases).

    Returns shape (N, 2).
    """
    t = np.asarray(t, dtype=float)
    out = np.zeros((t.size, 2))
    K = ph.approach_start.size

    # Figure out the robot's position phase-by-phase. Between phase k end
    # (time approach_start[k] + 1 + L[k]) and phase k+1 start
    # (approach_start[k+1]), the robot walks: boundary endpoint -> origin
    # (time 1) -> boundary start of next phase (time 1).
    # Because we constructed approach_start[k+1] = approach_start[k] + 1 + L[k] + 1,
    # the interval of "transit" has length (approach_start[k+1] - (approach_start[k]+1+L[k])) = 1,
    # which is exactly boundary->origin. Then next phase's approach handles
    # origin->next_boundary over its first 1 unit.
    for k in range(K):
        t_start = ph.approach_start[k]
        t_board = t_start + 1.0           # arrives on boundary
        t_scan_end = t_board + ph.L[k]    # finishes scan

        a_k = ph.a[k]; s_k = ph.s[k]; L_k = ph.L[k]
        approach_xy = np.array([np.cos(a_k), np.sin(a_k)])

        # [t_start, t_board]: approach from origin to boundary at angle a_k
        m = (t >= t_start) & (t <= t_board)
        out[m] = (t[m] - t_start)[:, None] * approach_xy

        # [t_board, t_scan_end]: scan on boundary
        m = (t > t_board) & (t <= t_scan_end)
        u = t[m] - t_board
        ang = a_k + s_k * u
        out[m, 0] = np.cos(ang); out[m, 1] = np.sin(ang)

        # [t_scan_end, t_scan_end + 1]: transit boundary endpoint -> origin
        if k + 1 < K:
            t_transit_end = t_scan_end + 1.0  # = ph.approach_start[k+1]
            end_ang = a_k + s_k * L_k
            end_xy = np.array([np.cos(end_ang), np.sin(end_ang)])
            m = (t > t_scan_end) & (t <= t_transit_end)
            # linear from end_xy to (0,0): at time t, fraction done = (t - t_scan_end)
            frac = (t[m] - t_scan_end)
            out[m] = (1.0 - frac)[:, None] * end_xy

    # After the last phase: park at scan endpoint
    t_last_end = ph.approach_start[-1] + 1.0 + ph.L[-1]
    end_ang = ph.a[-1] + ph.s[-1] * ph.L[-1]
    end_xy = np.array([np.cos(end_ang), np.sin(end_ang)])
    m = t > t_last_end
    out[m] = end_xy

    # Before t_start=0: at origin (all-zero already)
    return out


def _hit_times_phases(ph: Phases, thetas: np.ndarray) -> np.ndarray:
    """First-hit time for exit angle theta by this robot (min over phases).

    Shape (N_theta,). +inf if no phase hits theta.
    """
    out = np.full(thetas.size, np.inf)
    for k in range(ph.approach_start.size):
        a_k = ph.a[k]; s_k = ph.s[k]; L_k = ph.L[k]
        t_board = ph.approach_start[k] + 1.0
        delta = ((thetas - a_k) * s_k) % TWO_PI
        cov = delta <= L_k + 1e-12
        cand = t_board + delta
        out = np.where(cov & (cand < out), cand, out)
    return out


def worst_case_time(robot_phases: list[Phases], n_theta: int = 4001) -> tuple[float, float]:
    thetas = np.linspace(0.0, TWO_PI, n_theta, endpoint=False)
    # Hit times per robot, shape (3, N)
    ht = np.stack([_hit_times_phases(ph, thetas) for ph in robot_phases], axis=0)
    finite_any = np.isfinite(ht).any(axis=0)
    if not finite_any.all():
        return float("inf"), float("nan")
    t_found = ht.min(axis=0)                         # (N,)
    # Position of each robot at t_found: evaluate per-robot at the t_found grid
    positions = np.stack(
        [_positions_phases(ph, t_found) for ph in robot_phases], axis=0
    )                                                # (3, N, 2)
    exit_xy = np.stack([np.cos(thetas), np.sin(thetas)], axis=-1)   # (N, 2)
    chords = np.linalg.norm(positions - exit_xy[None, :, :], axis=-1)   # (3, N)
    evac = t_found + chords.max(axis=0)              # (N,)
    i = int(np.argmax(evac))
    return float(evac[i]), float(thetas[i])


def _coverage_ok(robot_phases: list[Phases], n_grid: int = 2001) -> bool:
    thetas = np.linspace(0.0, TWO_PI, n_grid, endpoint=False)
    ht = np.stack([_hit_times_phases(ph, thetas) for ph in robot_phases], axis=0)
    return bool(np.isfinite(ht).any(axis=0).all())


def _objective(x: np.ndarray, s1: int, s2a: int, s2b: int) -> float:
    a1, a2, b2, L0, L1, L2a, L2b = x
    robots = [
        _robot_phases_simple(0.0, L0, +1),
        _robot_phases_simple(a1, L1, s1),
        _robot_phases_two(a2, L2a, s2a, b2, L2b, s2b),
    ]
    # Cheap coverage check (coarse grid) before pulling out the big evaluator.
    if not _coverage_ok(robots, n_grid=1001):
        return 1000.0
    evac, _ = worst_case_time(robots, n_theta=1601)
    if not np.isfinite(evac):
        return 1000.0
    return evac


@dataclass
class SearchResult:
    x: np.ndarray
    discrete: tuple[int, int, int]
    worst_time: float
    worst_theta: float


def search(*, maxiter: int = 500, popsize: int = 40, seed: int = 0) -> SearchResult:
    bounds = [
        (0.0, TWO_PI),    # a1
        (0.0, TWO_PI),    # a2
        (0.0, TWO_PI),    # b2
        (0.1, TWO_PI),    # L0
        (0.1, TWO_PI),    # L1
        (0.0, TWO_PI),    # L2a  (allow 0 = skip first scan)
        (0.0, TWO_PI),    # L2b  (allow 0 = skip second scan)
    ]
    best: SearchResult | None = None
    for s1 in (+1, -1):
        for s2a in (+1, -1):
            for s2b in (+1, -1):
                res = differential_evolution(
                    _objective,
                    bounds=bounds,
                    args=(s1, s2a, s2b),
                    maxiter=maxiter,
                    popsize=popsize,
                    seed=seed,
                    tol=1e-8,
                    polish=True,
                    init="sobol",
                    updating="deferred",
                    workers=-1,
                )
                # Re-evaluate on a finer grid
                a1, a2, b2, L0, L1, L2a, L2b = res.x
                robots = [
                    _robot_phases_simple(0.0, L0, +1),
                    _robot_phases_simple(a1, L1, s1),
                    _robot_phases_two(a2, L2a, s2a, b2, L2b, s2b),
                ]
                if not _coverage_ok(robots, n_grid=4001):
                    continue
                evac, theta = worst_case_time(robots, n_theta=8001)
                print(f"  s=({s1:+d},{s2a:+d},{s2b:+d}): worst={evac:.6f}  theta={theta:.4f}")
                cand = SearchResult(res.x, (s1, s2a, s2b), evac, theta)
                if best is None or cand.worst_time < best.worst_time:
                    best = cand
    assert best is not None
    return best


if __name__ == "__main__":
    print("Family B search (robot 2 has a two-phase trajectory)...")
    best = search(maxiter=500, popsize=40, seed=0)
    print()
    a1, a2, b2, L0, L1, L2a, L2b = best.x
    s1, s2a, s2b = best.discrete
    print(f"Best worst-case evacuation time: {best.worst_time:.6f}")
    print(f"Worst exit theta                = {best.worst_theta:.6f} rad")
    print(f"Naive baseline                  : {1.0 + 2.0*np.pi/3 + 2*np.sin(np.pi/3):.6f}")
    print(f"Czyzowicz et al. A3 (k=3)       : {4*np.pi/9 + (2*np.sqrt(3)+5)/3 + 1/600:.6f}")
    print(f"LB (Thm 7)                      : ~4.159  (closed form in paper)")
    print()
    print(f"  robot 0: a=0         L={L0:.6f}  s=+1")
    print(f"  robot 1: a={a1:.6f}  L={L1:.6f}  s={s1:+d}")
    print(f"  robot 2: a={a2:.6f}  La={L2a:.6f}  sa={s2a:+d}")
    print(f"           b={b2:.6f}  Lb={L2b:.6f}  sb={s2b:+d}")
