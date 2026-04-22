"""Brandt et al. (CIAC 2017) F2F two-robot disk evacuation + ε-cut refinement.

Algorithm A(y, alpha, d): both robots go radially from M = (0,0) to the
starting point A = (1,0) of the perimeter search (time 1). They then
split: r_1 walks CCW, r_2 walks CW along the perimeter. After arc-length
y each, r_1 reaches C = e(y) and r_2 reaches B = e(-y). Each performs a
symmetric linear cut of depth d at angle alpha measured from the chord
BC, goes to the tip and returns; then continues along the perimeter
toward D = e(pi) (antipode of A).

Meeting protocol (Brandt et al. Section 2): if robot r_i finds the exit
at time t at position X, it computes the shortest intercept distance x
to r_{3-i}'s scheduled future path: x solves x = |X - r_{3-i}(t + x)|,
x >= 0 minimal. Evacuation time = t + 2 x.

This module provides:
    - A trajectory representation with phases for the base algorithm.
    - An evaluator computing evac(theta) for given parameters.
    - An optional "epsilon-cut" refinement: r_2 takes an additional
      small cut of depth d_2 at perimeter arc-position y_2 > y.

Target: reproduce Brandt's 5.625 at (y*, alpha*=pi/4, d*) and probe
whether an epsilon-cut gives a measurable improvement.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Sequence
import numpy as np
from scipy.optimize import brentq, minimize

TWO_PI = 2.0 * np.pi


# --- Trajectory segments -------------------------------------------------

@dataclass
class Segment:
    t0: float
    t1: float
    p0: np.ndarray
    p1: np.ndarray
    kind: str  # "line" or "arc"
    # for arcs:
    center: np.ndarray = None
    radius: float = 1.0
    theta0: float = 0.0    # starting angle on the arc
    thetaF: float = 0.0    # ending angle on the arc

    def position(self, t: float) -> np.ndarray:
        if t <= self.t0:
            return self.p0.copy()
        if t >= self.t1:
            return self.p1.copy()
        if self.kind == "line":
            frac = (t - self.t0) / (self.t1 - self.t0)
            return self.p0 + frac * (self.p1 - self.p0)
        else:
            u = (t - self.t0) / (self.t1 - self.t0)
            theta = self.theta0 + u * (self.thetaF - self.theta0)
            return self.center + self.radius * np.array([np.cos(theta), np.sin(theta)])

    def hit_time_to(self, target: np.ndarray, tol: float = 1e-10) -> float:
        """Earliest time the segment is at `target`. Returns inf if never."""
        if self.kind == "line":
            d = self.p1 - self.p0
            L2 = float(np.dot(d, d))
            if L2 < 1e-15:
                if np.linalg.norm(self.p0 - target) < tol:
                    return self.t0
                return float("inf")
            u = float(np.dot(target - self.p0, d) / L2)
            closest = self.p0 + u * d
            if 0.0 - 1e-12 <= u <= 1.0 + 1e-12 and np.linalg.norm(closest - target) < tol:
                return self.t0 + max(0.0, min(1.0, u)) * (self.t1 - self.t0)
            return float("inf")
        else:
            # arc on unit circle: check if target is on circle with matching angle
            if abs(np.linalg.norm(target - self.center) - self.radius) > tol:
                return float("inf")
            target_theta = np.arctan2(target[1] - self.center[1],
                                       target[0] - self.center[0])
            # normalize to [theta0, thetaF] range
            th0 = self.theta0
            thF = self.thetaF
            if thF == th0:
                if abs(target_theta - th0) < tol:
                    return self.t0
                return float("inf")
            # fraction along the arc in its sweep direction
            if thF > th0:
                u = (target_theta - th0) / (thF - th0)
                if u < 0:
                    u += TWO_PI / (thF - th0)
            else:
                u = (target_theta - th0) / (thF - th0)
                if u < 0:
                    u += TWO_PI / abs(thF - th0)
            if 0.0 - 1e-12 <= u <= 1.0 + 1e-12:
                return self.t0 + max(0.0, min(1.0, u)) * (self.t1 - self.t0)
            return float("inf")


class Trajectory:
    def __init__(self, segments: List[Segment]):
        self.segments = segments

    def position(self, t: float) -> np.ndarray:
        if t <= self.segments[0].t0:
            return self.segments[0].p0.copy()
        for seg in self.segments:
            if seg.t0 <= t <= seg.t1:
                return seg.position(t)
        return self.segments[-1].p1.copy()

    def hit_time(self, target: np.ndarray) -> float:
        best = float("inf")
        for seg in self.segments:
            t = seg.hit_time_to(target)
            if t < best:
                best = t
        return best


# --- Build Brandt A(y, alpha, d) -----------------------------------------

def unit_circle_point(theta: float) -> np.ndarray:
    return np.array([np.cos(theta), np.sin(theta)])


def build_brandt_trajectories(
    y: float, alpha: float, d: float,
    y2: float = 0.0, d2: float = 0.0,
) -> tuple[Trajectory, Trajectory]:
    """Build r_1 (CCW) and r_2 (CW) trajectories for Brandt's A(y, alpha, d).
    If y2 > y and d2 > 0, r_2 (and r_1 symmetrically) gets an additional
    epsilon-cut of depth d2 at perimeter arc-distance y2 from A.
    """
    M = np.zeros(2)
    A = unit_circle_point(0.0)
    C = unit_circle_point(y)   # r_1 scan endpoint (CCW)
    B = unit_circle_point(-y)  # r_2 scan endpoint (CW)

    # Chord BC direction (unit vector from B to C):
    BC_vec = C - B
    BC_unit = BC_vec / np.linalg.norm(BC_vec)
    # Inward normal to BC pointing toward origin side (perpendicular, chosen
    # so the cut goes inward toward the disk interior).
    # BC is horizontal only when y = pi/2; in general it has some tilt.
    # The inward normal (rotated 90 deg CCW from BC direction):
    normal = np.array([-BC_unit[1], BC_unit[0]])
    # Ensure normal points toward the interior (center). Check sign:
    if np.dot(normal, M - (B + C)/2) < 0:
        normal = -normal

    # Cut directions from C (r_1) and B (r_2):
    # Each cut is at angle alpha from the chord BC, going into the disk.
    # "angle alpha from chord BC" interpreted as: cut direction is BC_unit
    # rotated by angle alpha toward the normal. Specifically, r_1's cut
    # direction is (cos alpha) * (-BC_unit) + (sin alpha) * normal
    # (since r_1 is at C and moves toward B direction mixed with inward).
    r1_cut_dir = -np.cos(alpha) * BC_unit + np.sin(alpha) * normal
    r2_cut_dir =  np.cos(alpha) * BC_unit + np.sin(alpha) * normal

    # Key time stamps for r_1
    t_leave_M = 0.0
    t_at_A_r1 = 1.0
    t_at_C_r1 = 1.0 + y
    cut_tip_r1 = C + d * r1_cut_dir
    t_cut_tip_r1 = t_at_C_r1 + d
    t_back_C_r1 = t_cut_tip_r1 + d

    # r_1 continues CCW on perimeter past C
    segments_r1: List[Segment] = []
    segments_r1.append(Segment(t0=0.0, t1=1.0, p0=M, p1=A, kind="line"))
    segments_r1.append(Segment(
        t0=1.0, t1=t_at_C_r1, p0=A, p1=C, kind="arc",
        center=M, radius=1.0, theta0=0.0, thetaF=y,
    ))
    segments_r1.append(Segment(
        t0=t_at_C_r1, t1=t_cut_tip_r1, p0=C, p1=cut_tip_r1, kind="line",
    ))
    segments_r1.append(Segment(
        t0=t_cut_tip_r1, t1=t_back_C_r1, p0=cut_tip_r1, p1=C, kind="line",
    ))

    # Optional epsilon-cut for r_1 (symmetric to r_2's):
    if d2 > 0 and y2 > y:
        # r_1 walks from C (angle y) CCW to perimeter point at angle y2
        C2 = unit_circle_point(y2)
        t_at_C2_r1 = t_back_C_r1 + (y2 - y)
        segments_r1.append(Segment(
            t0=t_back_C_r1, t1=t_at_C2_r1, p0=C, p1=C2, kind="arc",
            center=M, radius=1.0, theta0=y, thetaF=y2,
        ))
        # eps-cut: inward direction at angle y2. The cut tangent at y2:
        # along the perimeter direction, rotated by alpha inward.
        # Inward normal at C2 is -C2 (pointing from boundary to center).
        # Tangent at C2 (CCW): (-sin y2, cos y2).
        tang2 = np.array([-np.sin(y2), np.cos(y2)])
        # cut direction at angle alpha from tangent toward interior:
        r1_cut2_dir = np.cos(alpha) * tang2 + np.sin(alpha) * (-C2)
        cut2_tip_r1 = C2 + d2 * r1_cut2_dir
        t_cut2_tip_r1 = t_at_C2_r1 + d2
        t_back_C2_r1 = t_cut2_tip_r1 + d2
        segments_r1.append(Segment(
            t0=t_at_C2_r1, t1=t_cut2_tip_r1, p0=C2, p1=cut2_tip_r1, kind="line",
        ))
        segments_r1.append(Segment(
            t0=t_cut2_tip_r1, t1=t_back_C2_r1, p0=cut2_tip_r1, p1=C2, kind="line",
        ))
        last_boundary_time = t_back_C2_r1
        last_boundary_angle = y2
        last_boundary_point = C2
    else:
        last_boundary_time = t_back_C_r1
        last_boundary_angle = y
        last_boundary_point = C

    # Continue CCW from last boundary point to D = e(pi)
    D = unit_circle_point(np.pi)
    t_at_D_r1 = last_boundary_time + (np.pi - last_boundary_angle)
    segments_r1.append(Segment(
        t0=last_boundary_time, t1=t_at_D_r1,
        p0=last_boundary_point, p1=D, kind="arc",
        center=M, radius=1.0, theta0=last_boundary_angle, thetaF=np.pi,
    ))

    # Similarly for r_2 (CW, mirrored across x-axis)
    t_at_B_r2 = 1.0 + y
    cut_tip_r2 = B + d * r2_cut_dir
    t_cut_tip_r2 = t_at_B_r2 + d
    t_back_B_r2 = t_cut_tip_r2 + d

    segments_r2: List[Segment] = []
    segments_r2.append(Segment(t0=0.0, t1=1.0, p0=M, p1=A, kind="line"))
    segments_r2.append(Segment(
        t0=1.0, t1=t_at_B_r2, p0=A, p1=B, kind="arc",
        center=M, radius=1.0, theta0=0.0, thetaF=-y,
    ))
    segments_r2.append(Segment(
        t0=t_at_B_r2, t1=t_cut_tip_r2, p0=B, p1=cut_tip_r2, kind="line",
    ))
    segments_r2.append(Segment(
        t0=t_cut_tip_r2, t1=t_back_B_r2, p0=cut_tip_r2, p1=B, kind="line",
    ))

    if d2 > 0 and y2 > y:
        B2 = unit_circle_point(-y2)
        t_at_B2_r2 = t_back_B_r2 + (y2 - y)
        segments_r2.append(Segment(
            t0=t_back_B_r2, t1=t_at_B2_r2, p0=B, p1=B2, kind="arc",
            center=M, radius=1.0, theta0=-y, thetaF=-y2,
        ))
        tang2 = np.array([np.sin(y2), np.cos(y2)])  # CW tangent at -y2
        r2_cut2_dir = np.cos(alpha) * tang2 + np.sin(alpha) * (-B2)
        cut2_tip_r2 = B2 + d2 * r2_cut2_dir
        t_cut2_tip_r2 = t_at_B2_r2 + d2
        t_back_B2_r2 = t_cut2_tip_r2 + d2
        segments_r2.append(Segment(
            t0=t_at_B2_r2, t1=t_cut2_tip_r2, p0=B2, p1=cut2_tip_r2, kind="line",
        ))
        segments_r2.append(Segment(
            t0=t_cut2_tip_r2, t1=t_back_B2_r2, p0=cut2_tip_r2, p1=B2, kind="line",
        ))
        last_boundary_time2 = t_back_B2_r2
        last_boundary_angle2 = -y2
        last_boundary_point2 = B2
    else:
        last_boundary_time2 = t_back_B_r2
        last_boundary_angle2 = -y
        last_boundary_point2 = B

    # Continue CW to D = e(pi) = e(-pi)
    t_at_D_r2 = last_boundary_time2 + (np.pi - y if d2 == 0 else np.pi - y2)
    segments_r2.append(Segment(
        t0=last_boundary_time2, t1=t_at_D_r2,
        p0=last_boundary_point2, p1=D, kind="arc",
        center=M, radius=1.0,
        theta0=last_boundary_angle2, thetaF=-np.pi,
    ))

    return Trajectory(segments_r1), Trajectory(segments_r2)


# --- F2F evacuation: solve fixed-point x = |X - r_other(t + x)| ---------

def f2f_evac_time(X: np.ndarray, t_found: float, other: Trajectory) -> float:
    """Compute evac time = t_found + 2 x where x is the shortest intercept
    distance. We solve x = |X - other(t_found+x)| at segment boundaries.
    Within each segment of `other`, other(t_found+x) is linear or circular in x.
    """

    def gap(x: float) -> float:
        p = other.position(t_found + x)
        return x - float(np.linalg.norm(X - p))

    # Search with a moderate grid; refine with brentq at sign changes.
    best_x = float("inf")
    # Use 201 points covering [0, 10] — more than enough for trajectory length.
    xs = np.linspace(0, 10, 201)
    gs = np.array([gap(xv) for xv in xs])
    for i in range(len(xs) - 1):
        if gs[i] <= 0 and gs[i + 1] >= 0:
            try:
                x_root = brentq(gap, xs[i], xs[i + 1], xtol=1e-10, rtol=1e-10)
                if x_root < best_x:
                    best_x = x_root
                    break
            except ValueError:
                continue
    if not np.isfinite(best_x):
        return float("inf")
    return t_found + 2.0 * best_x


def evac_for_exit(theta_exit: float, r1: Trajectory, r2: Trajectory) -> float:
    """Evac time for exit at angle theta_exit (on perimeter)."""
    X = unit_circle_point(theta_exit)
    t1 = r1.hit_time(X)
    t2 = r2.hit_time(X)
    if t1 < t2:
        # r_1 is finder
        return f2f_evac_time(X, t1, r2)
    else:
        return f2f_evac_time(X, t2, r1)


def worst_case_evac(y, alpha, d, y2=0.0, d2=0.0, n_theta: int = 2001) -> tuple[float, float]:
    r1, r2 = build_brandt_trajectories(y, alpha, d, y2, d2)
    thetas = np.linspace(-np.pi, np.pi, n_theta, endpoint=False)
    evacs = np.array([evac_for_exit(th, r1, r2) for th in thetas])
    i = int(np.argmax(evacs))
    return float(evacs[i]), float(thetas[i])


if __name__ == "__main__":
    # Brandt's claimed optimum
    y = 2.62843
    alpha = np.pi / 4.0
    d = 0.48793

    print("Verifying Brandt et al. A(y, alpha, d) at their claimed optimum:")
    T, theta_w = worst_case_evac(y, alpha, d, n_theta=4001)
    print(f"  y = {y}, alpha = pi/4, d = {d}")
    print(f"  Worst evac = {T:.6f}  at theta_exit = {theta_w:.6f}")
    print(f"  (Expected ~5.625)")
    print()

    # Joint optimisation over (y, d, y2, d2) with alpha = pi/4 fixed
    print("Joint optimisation (y, d, y2, d2):")
    print("  Starting from Brandt base + small eps-cut seed.")

    def obj(x):
        y, d, y2, d2 = x
        if not (1.5 < y < np.pi - 0.01):
            return 100.0
        if not (0.01 < d < 1.5):
            return 100.0
        if not (y + 0.01 < y2 < np.pi - 0.005):
            return 100.0
        if not (0.0 <= d2 < 0.3):
            return 100.0
        try:
            T, _ = worst_case_evac(y, np.pi / 4, d, y2=y2, d2=d2, n_theta=801)
        except Exception:
            return 100.0
        if not np.isfinite(T):
            return 100.0
        return T

    # Seed near Brandt + tiny eps-cut
    x0 = np.array([2.62843, 0.48793, 2.95, 0.02])
    T0 = obj(x0)
    print(f"  seed (y={x0[0]:.4f}, d={x0[1]:.4f}, y2={x0[2]:.4f}, d2={x0[3]:.4f}): T = {T0:.6f}")

    res = minimize(obj, x0, method="Nelder-Mead",
                   options={"xatol": 1e-6, "fatol": 1e-7, "maxiter": 400, "adaptive": True})
    y_o, d_o, y2_o, d2_o = res.x
    print()
    print(f"  optimum: y={y_o:.6f}, d={d_o:.6f}, y2={y2_o:.6f}, d2={d2_o:.6f}")
    # Fine-grid verify
    T_fine, theta_w = worst_case_evac(y_o, np.pi / 4, d_o, y2=y2_o, d2=d2_o, n_theta=8001)
    print(f"  fine-grid T = {T_fine:.8f}  at theta = {theta_w:.6f}")
    delta = 5.624897 - T_fine   # compare against our base reproduction
    print(f"  improvement over Brandt (our base 5.624897): {delta:+.8f}")
