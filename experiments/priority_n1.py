"""Priority evacuation with n = 1 servant (queen + 1 servant), wireless.

Model:
    Two robots Q (queen) and S (servant) start at the origin of the unit
    disk. Each moves at unit speed. An exit sits at an unknown angle
    theta in [0, 2 pi). Wireless communication: the first robot to visit
    e(theta) broadcasts instantly.

    Queen's evacuation time for exit theta:
        T_Q(theta) = t_find(theta) + |Q(t_find(theta)) - e(theta)|
    where t_find(theta) = min over i in {Q, S} of first visit time to
    e(theta). If Q is the finder, the second term is 0.

    Worst-case objective: T = sup_theta T_Q(theta).

This module provides (a) a simulator evaluating T for any pair of
trajectories Q(.), S(.) specified segment-by-segment; (b) closed-form
candidates for the "split" family where both robots scan disjoint arcs;
(c) local search for the optimal split parameter.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np

TWO_PI = 2.0 * np.pi


# --- Piecewise-linear trajectory representation --------------------------

@dataclass
class LinearSegment:
    """A linear segment traversed at unit speed: from start to end, starting
    at time t0 and ending at t0 + |end - start|.
    """
    t0: float
    start: np.ndarray  # shape (2,)
    end: np.ndarray    # shape (2,)

    @property
    def length(self) -> float:
        return float(np.linalg.norm(self.end - self.start))

    @property
    def t1(self) -> float:
        return self.t0 + self.length

    def position(self, t: float) -> np.ndarray:
        if self.length < 1e-15:
            return self.start.copy()
        u = (t - self.t0) / self.length
        u = max(0.0, min(1.0, u))
        return self.start + u * (self.end - self.start)


@dataclass
class ArcSegment:
    """A unit-speed arc along the boundary: from angle a at time t0, direction
    sign s (+1 CCW, -1 CW), length L (arc-length). The robot is on the
    boundary throughout, at angle a + s*(t - t0) for t in [t0, t0 + L].
    """
    t0: float
    a: float        # starting angle on boundary
    s: int          # +1 CCW, -1 CW
    L: float        # arc length (= time spent on this arc)

    @property
    def t1(self) -> float:
        return self.t0 + self.L

    def position(self, t: float) -> np.ndarray:
        u = max(0.0, min(self.L, t - self.t0))
        ang = self.a + self.s * u
        return np.array([np.cos(ang), np.sin(ang)])


Segment = LinearSegment | ArcSegment


@dataclass
class Trajectory:
    segments: list[Segment]

    def position(self, t: float) -> np.ndarray:
        if t <= 0:
            return self.segments[0].position(0.0) if self.segments else np.zeros(2)
        for seg in self.segments:
            if t <= seg.t1:
                return seg.position(t)
        # past end: parked at last segment's endpoint
        last = self.segments[-1]
        return last.position(last.t1)

    def first_hit_boundary_angle(self, theta: float, tol: float = 1e-12) -> float:
        """Return the first time at which the trajectory is on the boundary at
        angle theta, or +infty if never."""
        best = float("inf")
        for seg in self.segments:
            if isinstance(seg, ArcSegment):
                # Arc covers angles {a + s*u : u in [0, L]} mod 2 pi.
                # Find u such that (a + s*u) mod 2 pi == theta.
                delta = ((theta - seg.a) * seg.s) % TWO_PI
                if delta <= seg.L + tol:
                    t_hit = seg.t0 + delta
                    if t_hit < best:
                        best = t_hit
            else:
                # Linear segment: solve |seg.position(t) - e(theta)| = 0
                # exactly. For a generic segment, this amounts to checking if
                # the segment touches the point e(theta).
                target = np.array([np.cos(theta), np.sin(theta)])
                # project target onto segment
                d = seg.end - seg.start
                dl = np.linalg.norm(d)
                if dl < 1e-15:
                    if np.linalg.norm(seg.start - target) < tol:
                        if seg.t0 < best:
                            best = seg.t0
                    continue
                u = float(np.dot(target - seg.start, d) / (dl * dl))
                if 0 <= u <= 1:
                    closest = seg.start + u * d
                    if np.linalg.norm(closest - target) < tol:
                        t_hit = seg.t0 + u * dl
                        if t_hit < best:
                            best = t_hit
        return best


# --- Evacuation-time evaluation -----------------------------------------

def queen_time_for_exit(Q: Trajectory, S: Trajectory, theta: float) -> float:
    """Queen's evacuation time for exit angle theta."""
    tQ_hit = Q.first_hit_boundary_angle(theta)
    tS_hit = S.first_hit_boundary_angle(theta)
    t_find = min(tQ_hit, tS_hit)
    if not np.isfinite(t_find):
        return float("inf")
    exit_xy = np.array([np.cos(theta), np.sin(theta)])
    Q_at_find = Q.position(t_find)
    chord = float(np.linalg.norm(Q_at_find - exit_xy))
    return t_find + chord


def worst_queen_time(
    Q: Trajectory,
    S: Trajectory,
    *,
    n_theta: int = 40001,
) -> tuple[float, float]:
    thetas = np.linspace(0.0, TWO_PI, n_theta, endpoint=False)
    times = np.array([queen_time_for_exit(Q, S, th) for th in thetas])
    i = int(np.argmax(times))
    return float(times[i]), float(thetas[i])


# --- Candidate algorithm: "split scan" with Q continuing CCW ------------

def split_scan_trajectories(b: float) -> tuple[Trajectory, Trajectory]:
    """Q walks to angle 0 (time 1), scans CCW to angle b (time b), then
    continues CCW along the boundary to angle 2 pi. S walks to angle b
    (time 1), scans CCW to angle 2 pi (time 2 pi - b), then parks.
    """
    origin = np.zeros(2)
    e0 = np.array([1.0, 0.0])
    eb = np.array([np.cos(b), np.sin(b)])
    Q = Trajectory(segments=[
        LinearSegment(t0=0.0, start=origin, end=e0),
        ArcSegment(t0=1.0, a=0.0, s=+1, L=TWO_PI),     # covers [0, 2 pi)
    ])
    S = Trajectory(segments=[
        LinearSegment(t0=0.0, start=origin, end=eb),
        ArcSegment(t0=1.0, a=b, s=+1, L=TWO_PI - b),
    ])
    return Q, S


# --- Candidate algorithm: "Q continues past 2 pi, then chords back" -----

def split_then_return(b: float, *, return_target: float | None = None) -> tuple[Trajectory, Trajectory]:
    """Q scans [0, b] CCW, then takes a direct chord from e(b) to a target
    boundary point at angle return_target (default: antipode 2 pi - b/2).
    S scans [b, 2 pi].
    """
    origin = np.zeros(2)
    e0 = np.array([1.0, 0.0])
    eb = np.array([np.cos(b), np.sin(b)])
    if return_target is None:
        return_target = TWO_PI - b / 2.0       # CW of x-axis by b/2
    ert = np.array([np.cos(return_target), np.sin(return_target)])
    Q = Trajectory(segments=[
        LinearSegment(t0=0.0, start=origin, end=e0),
        ArcSegment(t0=1.0, a=0.0, s=+1, L=b),
        LinearSegment(t0=1.0 + b, start=eb, end=ert),
    ])
    S = Trajectory(segments=[
        LinearSegment(t0=0.0, start=origin, end=np.array([np.cos(b), np.sin(b)])),
        ArcSegment(t0=1.0, a=b, s=+1, L=TWO_PI - b),
    ])
    return Q, S


# --- Closed-form balance for split-scan (Q continues CCW) ---------------

def split_ccw_balance_eq(b: float) -> float:
    """Balance condition: Q's worst in her arc (time 1+b) == Q's worst for
    S-found exit (time 1 + (2 pi - b) + 2 sin(b/2)).
    Zero of this function gives the optimal b.
    """
    return (1.0 + b) - (1.0 + TWO_PI - b + 2 * np.sin(b / 2.0))


def split_ccw_evac(b: float) -> float:
    """Worst-case queen evac time for the split-CCW strategy at parameter b."""
    return max(1.0 + b, 1.0 + TWO_PI - b + 2 * np.sin(b / 2.0))


# --- Main ----------------------------------------------------------------

# --- Two-robot wireless algorithm adapted for priority n=1 --------------

def two_robot_wireless_trajectories(L: float = 2 * np.pi / 3) -> tuple[Trajectory, Trajectory]:
    """Q and S both approach angle 0; Q scans CCW for length L; S scans CW
    for length L. Each continues in their own direction past their scan
    endpoint along the boundary. This is the standard 2-robot wireless
    algorithm; the optimal choice L = 2 pi / 3 gives worst max-arrival
    1 + 2 pi/3 + sqrt 3.
    """
    origin = np.zeros(2)
    e0 = np.array([1.0, 0.0])
    Q = Trajectory(segments=[
        LinearSegment(t0=0.0, start=origin, end=e0),
        ArcSegment(t0=1.0, a=0.0, s=+1, L=2 * TWO_PI),   # scan CCW + continue (overkill length)
    ])
    S = Trajectory(segments=[
        LinearSegment(t0=0.0, start=origin, end=e0),
        ArcSegment(t0=1.0, a=0.0, s=-1, L=2 * TWO_PI),   # scan CW + continue
    ])
    return Q, S


# --- Asymmetric variant: Q scans longer arc than S ----------------------

def asym_trajectories(b_Q: float, b_S: float) -> tuple[Trajectory, Trajectory]:
    """Q approaches angle 0, scans CCW length b_Q, continues CCW. S
    approaches angle 0, scans CW length b_S, continues CW.
    """
    origin = np.zeros(2)
    e0 = np.array([1.0, 0.0])
    Q = Trajectory(segments=[
        LinearSegment(t0=0.0, start=origin, end=e0),
        ArcSegment(t0=1.0, a=0.0, s=+1, L=2 * TWO_PI),
    ])
    S = Trajectory(segments=[
        LinearSegment(t0=0.0, start=origin, end=e0),
        ArcSegment(t0=1.0, a=0.0, s=-1, L=2 * TWO_PI),
    ])
    # (Here we let both scan infinitely long; b_Q, b_S only affect evaluation
    # indirectly via the resulting worst-case. This family is just the "both
    # go to e(0), opposite directions" — scan lengths don't change the
    # resulting Q trajectory. Included here as a sanity check.)
    return Q, S


if __name__ == "__main__":
    from scipy.optimize import brentq, minimize_scalar

    print("=== Priority evacuation with n = 1 servant ===")
    print()

    # Naive bounds
    print("Naive bounds:")
    print(f"  Both-wait:                 2 + 2 pi = {2 + TWO_PI:.6f}")
    print(f"  Q-tracks-S:                1 + 2 pi = {1 + TWO_PI:.6f}")
    print()

    # Split CCW: analytic
    b_opt = brentq(split_ccw_balance_eq, 0.1, TWO_PI - 0.1)
    T_opt = split_ccw_evac(b_opt)
    print("Split CCW strategy, analytic balance:")
    print(f"  b_opt = {b_opt:.8f}")
    print(f"  T     = {T_opt:.8f}")
    print()

    # Two-robot wireless adapted
    print("Two-robot wireless adapted (Q CCW, S CW, both from e(0)):")
    Q, S = two_robot_wireless_trajectories()
    T_sim, theta_w = worst_queen_time(Q, S, n_theta=200001)
    ref = 1 + 2 * np.pi / 3 + np.sqrt(3.0)
    print(f"  simulator worst T       = {T_sim:.8f}  at theta = {theta_w:.6f}")
    print(f"  reference 1 + 2pi/3 + √3 = {ref:.8f}")
    print()

    # Scan over L for the 2-robot wireless adapted: both scan to length L
    # then continue. Actually the family is invariant to L if we just set
    # scan = continue. We'd need a real separate-scan-then-chord family to
    # study this cleanly. Defer.
    print("Additional variations to explore next...")
