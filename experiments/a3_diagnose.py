"""Diagnose the optimum: which cases are simultaneously binding at y_opt?

At y = y_opt, check evac as a function of theta and report all local
maxima, so we can see if case 1 (r_2 bottleneck, theta = 2*pi/3 - y/2)
and case 3 (r_3 mid-redeploy bottleneck, theta in [1+y, pi - y/2])
tie at the optimum.
"""
from __future__ import annotations

import numpy as np

from experiments.wireless_k3_opt2 import worst_case_time
from experiments.a3_optimize import a3_robots_at

TWO_PI = 2.0 * np.pi


def evac_vs_theta(y: float, n: int = 400001):
    from experiments.wireless_k3_opt2 import (
        _hit_times_phases, _positions_phases
    )
    robots = a3_robots_at(y)
    thetas = np.linspace(0.0, TWO_PI, n, endpoint=False)
    ht = np.stack([_hit_times_phases(r, thetas) for r in robots], axis=0)
    t_found = ht.min(axis=0)
    positions = np.stack([_positions_phases(r, t_found) for r in robots], axis=0)
    exit_xy = np.stack([np.cos(thetas), np.sin(thetas)], axis=-1)
    chords = np.linalg.norm(positions - exit_xy[None, :, :], axis=-1)
    evac = t_found + chords.max(axis=0)
    return thetas, evac, t_found, chords


if __name__ == "__main__":
    y_opt = 1.2158578321
    thetas, evac, t_found, chords = evac_vs_theta(y_opt, n=400001)

    # Find local maxima: evac[i] > evac[i-1] and evac[i] >= evac[i+1].
    max_mask = (evac[1:-1] >= evac[:-2]) & (evac[1:-1] >= evac[2:])
    idx_peaks = np.where(max_mask)[0] + 1
    # Filter: keep only peaks above (global max - 0.01).
    global_max = evac.max()
    keep = idx_peaks[evac[idx_peaks] > global_max - 0.01]
    # Consolidate very close peaks.
    kept = []
    last = -np.inf
    for k in sorted(keep, key=lambda k: -evac[k]):
        if all(abs(thetas[k] - thetas[j]) > 0.01 for j in kept):
            kept.append(k)
    kept = sorted(kept, key=lambda k: thetas[k])

    print(f"y_opt = {y_opt:.10f}")
    print(f"Global max evac = {global_max:.10f}")
    print()
    print("Local maxima of evac(theta) (within 0.01 of global max):")
    for k in kept:
        # Determine which robot is bottleneck at this theta.
        bottleneck = int(np.argmax(chords[:, k]))
        print(f"  theta = {thetas[k]:.6f}  evac = {evac[k]:.10f}  "
              f"t_found = {t_found[k]:.6f}  bottleneck = r_{bottleneck+1}  "
              f"chord = {chords[bottleneck, k]:.6f}")
