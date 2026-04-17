"""Family F: A_3 skeleton but r_3 takes a direct chord from scan-end to the
redeploy boundary point, instead of going via the origin.

Rationale. In A_3 the "return to origin + go out again" costs 2 units of
time. A direct chord from r_3's scan endpoint (at angle 2*pi - y + y =
2*pi = 0 in our gauge, i.e. the point (1,0)) to the redeployment target
at angle alpha = pi - y/2 has length 2 sin((pi - y/2)/2) = 2 cos(y/4),
which is strictly less than 2 for y > 0. At y = y_opt ~ 1.216 the
savings is ~0.091 of transit time.

Per-timestep r_3 position:
    t in [0, 1]:              radial approach to B = (cos(2 pi - y), sin(2 pi - y))
    t in [1, 1 + y]:           scan CCW from B, ending at angle 2 pi -> position (1, 0)
    t in [1+y, 1+y+ell]:       linear chord (1,0) -> P = (cos alpha, sin alpha)
                               where ell = |(1,0) - P| = 2 cos(y/4) at alpha = pi - y/2
    t > 1+y+ell:               parked at P.
We keep alpha = pi - y/2 (antipodal to midpoint of AB) but let it vary
for local search. Same coverage structure as A_3 (tight, L_1 = L_2 =
pi - y/2, no second scan from r_3).

Free parameters: y (and optionally alpha as a free offset from pi - y/2).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from experiments.wireless_k3_opt2 import (
    _robot_phases_simple,
    worst_case_time as _wc_two_phase,  # we'll redefine custom one
    _hit_times_phases,
)

TWO_PI = 2.0 * np.pi
EVAC_A3 = 4.2185169939


# -----------------------------------------------------------------------------
# Custom r_3 that uses a direct chord rather than via the origin.
# -----------------------------------------------------------------------------


@dataclass
class R3Chord:
    """r_3: approach B (time 1) -> scan CCW length y (ending at (1,0))
    -> direct chord to boundary angle alpha (time ell = 2 sin(alpha/2)
    starting from (1,0), i.e. ell = 2 sin(alpha/2) since chord endpoints
    are on the unit circle separated by angle alpha) -> park."""
    y: float
    alpha: float                 # redeploy angle
    # derived:
    @property
    def aB(self) -> float:
        return (TWO_PI - self.y) % TWO_PI
    @property
    def t_board(self) -> float:
        return 1.0
    @property
    def t_scan_end(self) -> float:
        return 1.0 + self.y
    @property
    def chord_len(self) -> float:
        # chord (1,0) -> (cos alpha, sin alpha)
        return float(np.sqrt(2 - 2 * np.cos(self.alpha)))
    @property
    def t_arrive(self) -> float:
        return self.t_scan_end + self.chord_len

    def position(self, t: np.ndarray) -> np.ndarray:
        t = np.asarray(t, dtype=float)
        out = np.zeros((t.size, 2))
        ap = np.array([np.cos(self.aB), np.sin(self.aB)])  # boundary target
        # Phase 0: approach origin -> B
        m = (t >= 0) & (t <= self.t_board)
        out[m] = t[m, None] * ap
        # Phase 1: scan from aB CCW for time y
        m = (t > self.t_board) & (t <= self.t_scan_end)
        u = t[m] - self.t_board
        ang = self.aB + u
        out[m, 0] = np.cos(ang); out[m, 1] = np.sin(ang)
        # Phase 2: chord (1, 0) -> (cos alpha, sin alpha)
        end_ang = self.alpha
        start_xy = np.array([1.0, 0.0])     # = (cos(aB + y), sin(aB + y)) = (cos 2pi, sin 2pi)
        end_xy = np.array([np.cos(end_ang), np.sin(end_ang)])
        direction = end_xy - start_xy
        ell = self.chord_len
        m = (t > self.t_scan_end) & (t <= self.t_arrive)
        if ell > 0:
            frac = (t[m] - self.t_scan_end) / ell
        else:
            frac = np.zeros_like(t[m])
        out[m] = start_xy + frac[:, None] * direction
        # Phase 3: parked
        m = t > self.t_arrive
        out[m] = end_xy
        return out

    def hit_times(self, thetas: np.ndarray) -> np.ndarray:
        """Hit time for each theta: the only scanning phase is the
        boundary-scan phase 1, covering angles [aB, aB + y] mod 2 pi.
        """
        out = np.full(thetas.size, np.inf)
        delta = (thetas - self.aB) % TWO_PI     # CCW arc from aB
        covered = delta <= self.y + 1e-12
        out[covered] = self.t_board + delta[covered]
        return out


# -----------------------------------------------------------------------------


def robots_from(y, alpha_offset=0.0):
    """Build Family F with tight coverage L_1 = L_2 = pi - y/2 and alpha =
    pi - y/2 + alpha_offset."""
    L = np.pi - y / 2.0
    aB = (TWO_PI - y) % TWO_PI
    r1 = _robot_phases_simple(a=0.0, L=L, s=+1)
    r2 = _robot_phases_simple(a=aB, L=L, s=-1)
    r3 = R3Chord(y=y, alpha=np.pi - y/2.0 + alpha_offset)
    return [r1, r2, r3]


def worst_case_time(robots, *, n_theta: int = 4001):
    r1, r2, r3 = robots
    thetas = np.linspace(0.0, TWO_PI, n_theta, endpoint=False)
    ht1 = _hit_times_phases(r1, thetas)
    ht2 = _hit_times_phases(r2, thetas)
    ht3 = r3.hit_times(thetas)
    ht = np.stack([ht1, ht2, ht3], axis=0)
    covered = np.isfinite(ht).any(axis=0)
    if not covered.all():
        return float("inf"), float("nan")
    t_found = ht.min(axis=0)
    # Positions at t_found:
    from experiments.wireless_k3_opt2 import _positions_phases
    p1 = _positions_phases(r1, t_found)
    p2 = _positions_phases(r2, t_found)
    p3 = r3.position(t_found)
    exit_xy = np.stack([np.cos(thetas), np.sin(thetas)], axis=-1)
    chords = np.stack([
        np.linalg.norm(p1 - exit_xy, axis=-1),
        np.linalg.norm(p2 - exit_xy, axis=-1),
        np.linalg.norm(p3 - exit_xy, axis=-1),
    ], axis=0)
    evac = t_found + chords.max(axis=0)
    i = int(np.argmax(evac))
    return float(evac[i]), float(thetas[i])


if __name__ == "__main__":
    # Scalar search over y (alpha_offset = 0 -> alpha = pi - y/2).
    print("Family F: A_3 skeleton with r_3 direct chord (not via origin).")
    print(f"  A_3 benchmark: {EVAC_A3:.10f}")

    def obj_y(y):
        robots = robots_from(float(y), 0.0)
        evac, _ = worst_case_time(robots, n_theta=16001)
        return evac

    for y_probe in (1.0, 1.1, 1.15, 1.2, 1.22, 1.25, 1.3, 1.4):
        v = obj_y(y_probe)
        marker = "*" if v < EVAC_A3 else " "
        print(f"  y = {y_probe:.3f}  ->  evac = {v:.6f}  {marker}")

    res = minimize_scalar(obj_y, bracket=(1.0, 1.22, 1.4),
                          method="brent", tol=1e-10)
    y_opt = float(res.x)
    robots = robots_from(y_opt, 0.0)
    evac, theta = worst_case_time(robots, n_theta=200001)
    print()
    print(f"Optimal y: {y_opt:.10f}")
    print(f"Worst-case evac: {evac:.10f}  at theta = {theta:.6f}")
    print(f"Improvement over A_3: {EVAC_A3 - evac:.10f}")

    # Now allow alpha to float as well.
    def obj_both(x):
        y, ao = x
        robots = robots_from(float(y), float(ao))
        evac, _ = worst_case_time(robots, n_theta=16001)
        return evac

    res2 = minimize(obj_both, [y_opt, 0.0], method="Nelder-Mead",
                    options={"xatol": 1e-9, "fatol": 1e-10,
                             "maxiter": 20000, "adaptive": True})
    y2, ao2 = res2.x
    robots2 = robots_from(y2, ao2)
    evac2, theta2 = worst_case_time(robots2, n_theta=200001)
    print()
    print("With alpha free:")
    print(f"  y = {y2:.10f}  alpha_offset = {ao2:.10f}")
    print(f"  evac = {evac2:.10f}  at theta = {theta2:.6f}")
    print(f"  improvement over A_3 = {EVAC_A3 - evac2:.10f}")
