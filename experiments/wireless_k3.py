"""Wireless k-robot evacuation simulator, focused on k=3.

Model (Czyzowicz et al. 2014):
    - k unit-speed robots start at the origin of a unit disk.
    - An exit sits at an unknown angle theta on the boundary.
    - Wireless: once any robot's trajectory passes through the exit, it
      broadcasts instantly; every robot then walks straight to the exit.
    - Evacuation time = time last robot reaches the exit.

Each robot's deployment is specified as an explicit Trajectory r: [0, T] -> R^2,
unit-speed, starting at the origin. The simulator evaluates worst-case
evacuation time over exit angle theta.

Reference bounds (k=3):
    LB  : 3 + pi/3             ~= 4.04720
    UB  : 3 + pi/3 + O(k^{-4/3}) from Czyzowicz et al. (small additive slack)

Baseline algorithm (naive arc-partition):
    Robot i in {0,1,2} walks from origin to boundary angle 2*pi*i/3 in t in [0,1],
    then scans the arc [2*pi*i/3, 2*pi*i/3 + 2*pi/3] CCW in t in [1, 1+2*pi/3].
    Worst-case evac time = 1 + 2*pi/3 + 2*sin(pi/3) = 1 + 2*pi/3 + sqrt(3)
                         ~= 4.8258, clearly not optimal.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np

Trajectory = Callable[[float], np.ndarray]  # t -> array([x, y]), unit-speed


@dataclass
class WirelessResult:
    worst_exit_angle: float
    worst_time: float
    t_found: float
    bottleneck_robot: int


def _exit_xy(theta: float) -> np.ndarray:
    return np.array([np.cos(theta), np.sin(theta)])


def _first_hit_time(
    traj: Trajectory,
    theta: float,
    t_grid: np.ndarray,
    tol: float | None = None,
) -> float:
    """First t in t_grid at which traj(t) is within `tol` of the exit point.

    Default tolerance: 2x the grid spacing. Unit-speed trajectories shift by
    at most `dt` between samples, so this catches any true crossing.

    Returns +inf if the trajectory never reaches the exit within t_grid.
    """
    target = _exit_xy(theta)
    if tol is None:
        dt = float(t_grid[1] - t_grid[0]) if len(t_grid) > 1 else 1e-4
        tol = 2.0 * dt
    for t in t_grid:
        p = traj(float(t))
        if np.hypot(p[0] - target[0], p[1] - target[1]) < tol:
            return float(t)
    return float("inf")


def wireless_evac_time(
    trajectories: Sequence[Trajectory],
    theta: float,
    *,
    t_max: float = 8.0,
    n_t: int = 40001,
    hit_tol: float | None = None,
) -> tuple[float, float, int]:
    """Evacuation time in wireless k-robot model for exit at angle theta.

    Returns (evac_time, t_found, bottleneck_robot_index).
    """
    t_grid = np.linspace(0.0, t_max, n_t)
    hit_times = [_first_hit_time(tr, theta, t_grid, tol=hit_tol) for tr in trajectories]
    t_found = min(hit_times)
    if not np.isfinite(t_found):
        return float("inf"), float("inf"), -1

    target = _exit_xy(theta)
    chords = np.array([np.linalg.norm(tr(t_found) - target) for tr in trajectories])
    bottleneck = int(np.argmax(chords))
    return t_found + float(chords[bottleneck]), t_found, bottleneck


def worst_case(
    trajectories: Sequence[Trajectory],
    *,
    n_angles: int = 4001,
    **kwargs,
) -> WirelessResult:
    angles = np.linspace(0.0, 2.0 * np.pi, n_angles, endpoint=False)
    best = WirelessResult(0.0, -np.inf, 0.0, -1)
    for theta in angles:
        evac, t_found, r = wireless_evac_time(trajectories, float(theta), **kwargs)
        if evac > best.worst_time:
            best = WirelessResult(float(theta), evac, t_found, r)
    return best


# ----------------------------------------------------------------------------
# Baseline: naive arc-partition algorithm for k robots.
# ----------------------------------------------------------------------------

def make_arc_partition_trajectory(alpha: float, arc_len: float) -> Trajectory:
    """Walk radially to boundary at angle alpha in [0, 1], then scan CCW.

    The robot's boundary walk covers angles [alpha, alpha + arc_len] over
    t in [1, 1 + arc_len].
    """
    start_xy = _exit_xy(alpha)

    def r(t: float) -> np.ndarray:
        if t <= 1.0:
            return t * start_xy
        s = t - 1.0
        if s <= arc_len:
            return _exit_xy(alpha + s)
        # after the scan, park at the endpoint
        return _exit_xy(alpha + arc_len)

    return r


def naive_k_trajectories(k: int) -> list[Trajectory]:
    arc = 2.0 * np.pi / k
    return [make_arc_partition_trajectory(i * arc, arc) for i in range(k)]


if __name__ == "__main__":
    k = 3
    trajs = naive_k_trajectories(k)
    res = worst_case(trajs, n_angles=2001, n_t=20001)

    lb = 3.0 + np.pi / k
    expected_naive = 1.0 + 2.0 * np.pi / k + 2.0 * np.sin(np.pi / k)

    print(f"k = {k}")
    print(f"  LB (Czyzowicz)      : {lb:.6f}")
    print(f"  naive expected WC   : {expected_naive:.6f}  (= 1 + 2pi/k + 2 sin(pi/k))")
    print(f"  simulator worst WC  : {res.worst_time:.6f}")
    print(f"  worst exit angle    : {res.worst_exit_angle:.6f} rad")
    print(f"  t_found             : {res.t_found:.6f}")
    print(f"  bottleneck robot    : {res.bottleneck_robot}")
