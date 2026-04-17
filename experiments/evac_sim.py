"""Skeleton evacuation simulator.

Computes the worst-case evacuation time for a given pair of trajectories
on the unit disk. Intended as a starting point for numerical exploration
of algorithms in the F2F / wireless / SR models.

Trajectories are given as callables t -> (x, y) with unit-speed parameterization.
Evacuation time for exit at angle theta is:
    min over robots r of (time r reaches theta on the perimeter + cost to bring
                          partner to the exit once it's announced/met).

This file is intentionally minimal — fill in the algorithm-specific logic
per open problem under study.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

Trajectory = Callable[[float], np.ndarray]  # t -> array([x, y])


@dataclass
class EvacResult:
    worst_exit_angle: float
    worst_time: float


def worst_case_time(
    traj_a: Trajectory,
    traj_b: Trajectory,
    *,
    t_max: float = 10.0,
    n_angles: int = 2001,
    time_eval: Callable[[Trajectory, Trajectory, float], float] | None = None,
) -> EvacResult:
    """Scan exit angles in [0, 2*pi) and return the worst-case evacuation time.

    The caller supplies `time_eval(traj_a, traj_b, theta)` — the time at which
    the *last* robot reaches the exit at angle theta, assuming the algorithm
    rules of the model under study. Provide a model-specific implementation
    per experiment.
    """
    if time_eval is None:
        raise ValueError("pass a model-specific time_eval")

    angles = np.linspace(0.0, 2.0 * np.pi, n_angles, endpoint=False)
    times = np.array([time_eval(traj_a, traj_b, theta) for theta in angles])
    idx = int(np.argmax(times))
    return EvacResult(worst_exit_angle=float(angles[idx]), worst_time=float(times[idx]))


if __name__ == "__main__":
    print("Skeleton module — import and extend per experiment.")
