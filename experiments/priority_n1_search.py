"""Priority evacuation n=1: richer search for algorithms beating 4.826.

Strategy family:
    Q: radial to e(alpha_Q), chord to e(beta_Q), scan direction sQ, then
       continue scanning on boundary beyond. Possibly a second chord /
       scan phase.
    S: radial to e(alpha_S), scan direction sS, continue scanning.

We parametrise and numerically optimise over the continuous parameters,
testing whether the worst-case queen time T can drop below 4.826.

This code uses experiments.priority_n1 simulator.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize, differential_evolution

from experiments.priority_n1 import (
    Trajectory, LinearSegment, ArcSegment, worst_queen_time, TWO_PI,
)

REF_2WL = 1.0 + 2.0 * np.pi / 3.0 + np.sqrt(3.0)       # ~ 4.8264


def make_Q_chord_scan(alpha_Q: float, beta_Q: float, sQ: int) -> Trajectory:
    """Q: origin -> e(alpha_Q) [time 1] -> chord to e(beta_Q) -> scan
    direction sQ on the boundary starting at e(beta_Q)."""
    origin = np.zeros(2)
    ea = np.array([np.cos(alpha_Q), np.sin(alpha_Q)])
    eb = np.array([np.cos(beta_Q), np.sin(beta_Q)])
    chord_len = float(np.linalg.norm(ea - eb))
    segments = [
        LinearSegment(t0=0.0, start=origin, end=ea),
        LinearSegment(t0=1.0, start=ea, end=eb),
        ArcSegment(t0=1.0 + chord_len, a=beta_Q, s=sQ, L=2 * TWO_PI),
    ]
    return Trajectory(segments=segments)


def make_S_scan(alpha_S: float, sS: int) -> Trajectory:
    """S: origin -> e(alpha_S), then scan direction sS."""
    origin = np.zeros(2)
    ea = np.array([np.cos(alpha_S), np.sin(alpha_S)])
    return Trajectory(segments=[
        LinearSegment(t0=0.0, start=origin, end=ea),
        ArcSegment(t0=1.0, a=alpha_S, s=sS, L=2 * TWO_PI),
    ])


def objective(x, sQ, sS):
    alpha_Q, beta_Q, alpha_S = x
    Q = make_Q_chord_scan(alpha_Q, beta_Q, sQ)
    S = make_S_scan(alpha_S, sS)
    # Need a fine enough grid to catch narrow peaks. Use n_theta=20001.
    T, _ = worst_queen_time(Q, S, n_theta=20001)
    return T


def search_family():
    bounds = [
        (0.0, TWO_PI),     # alpha_Q
        (0.0, TWO_PI),     # beta_Q
        (0.0, TWO_PI),     # alpha_S
    ]
    results = []
    for sQ in (+1, -1):
        for sS in (+1, -1):
            res = differential_evolution(
                objective, bounds=bounds, args=(sQ, sS),
                maxiter=300, popsize=40, seed=0, tol=1e-8,
                polish=True, init="sobol", updating="deferred",
                workers=-1,
            )
            # Fine-grid verify
            Q = make_Q_chord_scan(*res.x[:2], sQ)
            S = make_S_scan(res.x[2], sS)
            T_fine, theta_w = worst_queen_time(Q, S, n_theta=40001)
            results.append((res.fun, res.x, (sQ, sS), T_fine, theta_w))
            print(f"  s = ({sQ:+d},{sS:+d}): DE {res.fun:.6f} -> fine {T_fine:.6f} at theta={theta_w:.4f}")
    return results


if __name__ == "__main__":
    print(f"Benchmark (2-robot wireless adapted): {REF_2WL:.8f}")
    print()
    print("Searching Family G: Q (chord + scan), S (scan)")
    results = search_family()
    print()
    best = min(results, key=lambda r: r[3])
    print(f"Best fine-grid T = {best[3]:.8f} (vs benchmark {REF_2WL:.6f})")
    print(f"  improvement: {REF_2WL - best[3]:+.6f}")
    print(f"  params: alpha_Q={best[1][0]:.4f}, beta_Q={best[1][1]:.4f}, alpha_S={best[1][2]:.4f}")
    print(f"  directions: {best[2]}")
    print(f"  worst theta: {best[4]:.4f}")
