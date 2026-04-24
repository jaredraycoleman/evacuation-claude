"""Test the generalized min-chord conjecture at various |E| values.

Claim: for each e in (0, 3 x^*], the min sup chord over all 3-robot
valid states with |E| = e equals h*(e), the balanced 2-coincident value.

Procedure at fixed e:
    1. Compute h*(e) from the balanced 2-coincident formula.
    2. Run broad DE minimax over states with |E| ~ e, using an analytic
       sup to avoid grid artefacts.
    3. Compare the two; if they match, conjecture confirmed at this e.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import differential_evolution

TWO_PI = 2.0 * math.pi
X_STAR = (2.0 / 3.0) * math.acos(-1.0 / 3.0)


def balanced_2coincident_chord(tau: float) -> float:
    """Sup chord for balanced 2-coincident with slack tau per robot."""
    r = 1.0 - tau
    L = X_STAR - tau
    return math.sqrt(r * r + 2 * r * math.cos(L / 2.0) + 1.0)


def state_tight_cgk(e: float) -> float:
    if e >= math.pi:
        return 1.0 + X_STAR + 2.0 * math.sin((TWO_PI - e) / 2.0)
    return 1.0 + X_STAR + 2.0


def sup_chord_analytic(robots: np.ndarray, arcs: list[tuple[float, float]]) -> float:
    """Compute sup over unexplored of max_j ||R_j - p||, using endpoints + robot antipodes."""
    # Candidate beta values: arc endpoints + antipodes of robot direction angles
    candidates: list[float] = []
    for (lo, hi) in arcs:
        candidates.append(lo % TWO_PI)
        candidates.append(hi % TWO_PI)
    for R in robots:
        # Robot's direction (if not origin)
        norm = np.linalg.norm(R)
        if norm > 1e-10:
            theta = math.atan2(R[1], R[0]) % TWO_PI
            candidates.append((theta + math.pi) % TWO_PI)

    # For each candidate, check if in unexplored (not strictly in any arc)
    def in_explored(beta: float) -> bool:
        beta = beta % TWO_PI
        for (lo, hi) in arcs:
            lo = lo % TWO_PI
            if hi <= TWO_PI:
                if lo < beta < hi:
                    return True
            else:
                if beta > lo or beta < hi - TWO_PI:
                    return True
        return False

    best = 0.0
    for beta in candidates:
        if in_explored(beta):
            continue
        p = np.array([math.cos(beta), math.sin(beta)])
        c = max(np.linalg.norm(R - p) for R in robots)
        if c > best:
            best = c
    return best


def min_sup_chord_at_e(target_e: float, n_restarts: int = 15) -> tuple[float, dict]:
    """Minimise sup chord over states with |E| = target_e (approximately),
    by constraining L_1 + L_2 + L_3 = target_e in a transformed param space.
    """
    def objective(params):
        # Params: alpha_j (3), s_j (3 floats to +-1), fraction_j for L_j (3, softmax)
        alphas = params[0:3] % TWO_PI
        ss = np.sign(params[3:6])
        ss[ss == 0] = 1.0
        # Softmax-like: L_j proportional to softmax, sum = target_e
        raw_L = np.abs(params[6:9])
        if raw_L.sum() > 0:
            L_j = raw_L / raw_L.sum() * target_e
        else:
            L_j = np.array([target_e / 3.0] * 3)
        # Cap each L_j at x_star
        L_j = np.minimum(L_j, X_STAR)
        # Rescale if capping violates sum
        if L_j.sum() > 0:
            L_j = L_j * (target_e / max(L_j.sum(), 1e-9))
            L_j = np.minimum(L_j, X_STAR)

        # tau_j = how much inward; capped by (x_star - L_j) and 1
        max_tau = np.minimum(X_STAR - L_j, 1.0)
        tau_frac = np.clip(params[9:12], 0.0, 1.0)
        taus = tau_frac * max_tau

        # Robots
        robots = np.array([
            (1.0 - taus[j]) * np.array([math.cos(alphas[j] + ss[j] * L_j[j]),
                                        math.sin(alphas[j] + ss[j] * L_j[j])])
            for j in range(3)
        ])

        # Arcs
        arcs = []
        for j in range(3):
            if ss[j] > 0:
                lo = alphas[j] % TWO_PI
            else:
                lo = (alphas[j] - L_j[j]) % TWO_PI
            arcs.append((lo, lo + L_j[j]))

        # Check arcs are disjoint (penalize if not, to ensure |E| = target_e)
        # Sort by lo
        norm_arcs = []
        for (lo, hi) in arcs:
            lo = lo % TWO_PI
            if hi <= TWO_PI:
                norm_arcs.append((lo, hi))
            else:
                norm_arcs.append((lo, TWO_PI))
                norm_arcs.append((0.0, hi - TWO_PI))
        norm_arcs.sort()
        overlap = 0.0
        for i in range(len(norm_arcs) - 1):
            overlap += max(0, norm_arcs[i][1] - norm_arcs[i + 1][0])

        return sup_chord_analytic(robots, arcs) + 100.0 * overlap

    bounds = [(0.0, TWO_PI)] * 3 + [(-1.0, 1.0)] * 3 + [(0.0, 1.0)] * 6

    best_val = float("inf")
    best_info = None
    for trial in range(n_restarts):
        res = differential_evolution(
            objective, bounds, seed=trial, maxiter=300, popsize=40, tol=1e-10,
            init="sobol", polish=True, workers=1,
        )
        if res.fun < best_val:
            best_val = res.fun
            best_info = {"val": res.fun, "params": res.x.copy()}
    return best_val, best_info or {}


def main():
    print(f"x^* = {X_STAR:.10f}, 3 x^* = {3 * X_STAR:.10f}")
    print()
    print("Test: at each e, compare h*(e) (balanced 2-coincident) to")
    print("numerically-minimised sup chord over all states.")
    print()
    print("  e        | tau=(3x*-e)/3 | h*(e)       | numerical min | match?")
    print("  ---------+---------------+-------------+---------------+--------")

    for e in [3 * X_STAR - 0.0, 3 * X_STAR - 0.02, 3 * X_STAR - 0.05,
              3 * X_STAR - 0.1, 3 * X_STAR - 0.2, 3 * X_STAR - 0.5,
              math.pi + 0.1, math.pi]:
        tau = (3 * X_STAR - e) / 3.0
        h_star = balanced_2coincident_chord(tau)
        num_min, info = min_sup_chord_at_e(e, n_restarts=10)
        match = abs(num_min - h_star) < 0.01
        print(f"  {e:.4f}  | {tau:.5f}       | {h_star:.7f} | {num_min:.7f}      | "
              f"{'✓' if match else '✗'} (diff {num_min - h_star:+.5f})")


if __name__ == "__main__":
    main()
