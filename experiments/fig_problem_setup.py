"""Figure: the k=3 wireless evacuation problem setup.

Produces figures/problem_setup.pdf showing:
    - the unit disk with centre marked,
    - three robots stacked at the centre,
    - an unknown exit on the boundary, labelled well outside the disk,
    - a generic robot trajectory (radial approach + boundary scan + chord
      to the exit), with each phase labelled in a clearly separated spot.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
FIG_DIR.mkdir(exist_ok=True)

BLUE = "#3a7ab0"
RED  = "#c0392b"


def main() -> None:
    fig, ax = plt.subplots(figsize=(6.0, 5.2))

    # -- Disk outline --
    tc = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(tc), np.sin(tc), color="black", lw=1.3)

    # -- Centre dot + robot labels stacked above it --
    ax.scatter([0], [0], s=22, color="black", zorder=5)
    ax.annotate("origin\n(all three robots start here)",
                xy=(0, 0), xytext=(-0.04, -0.35), ha="center",
                fontsize=9,
                arrowprops=dict(arrowstyle="-", color="black", lw=0.6))

    # -- Trajectory geometry --
    # Approach angle (where the robot hits the boundary).
    alpha = 0.30
    # Scan goes CCW from alpha for length L_scan; keep scan strictly before exit
    # so the "chord" segment is geometrically visible.
    L_scan = 0.55
    scan_end_angle = alpha + L_scan
    # Exit well past scan end.
    theta_ex = 1.55

    approach_end = np.array([np.cos(alpha), np.sin(alpha)])
    scan_end     = np.array([np.cos(scan_end_angle), np.sin(scan_end_angle)])
    exit_xy      = np.array([np.cos(theta_ex), np.sin(theta_ex)])

    # -- Radial approach (blue line) --
    ax.plot([0, approach_end[0]], [0, approach_end[1]],
            color=BLUE, lw=1.6)
    # label along the midpoint, offset outward (perpendicular to the segment)
    mid_app = 0.5 * approach_end
    perp_app = np.array([-np.sin(alpha), np.cos(alpha)])
    lbl_app = mid_app + 0.16 * perp_app
    ax.annotate("radial approach",
                xy=mid_app, xytext=lbl_app,
                fontsize=9, color=BLUE, ha="center",
                arrowprops=dict(arrowstyle="-", color=BLUE, lw=0.5))

    # -- Boundary scan (thicker blue arc) --
    scan_t = np.linspace(alpha, scan_end_angle, 120)
    ax.plot(np.cos(scan_t), np.sin(scan_t), color=BLUE, lw=2.4)
    # Label the lower part of the arc, well outside the disk to the right.
    scan_label_angle = alpha + 0.15       # near the lower end of the arc
    tip_scan = 1.01 * np.array([np.cos(scan_label_angle),
                                np.sin(scan_label_angle)])
    lbl_scan = np.array([1.35, 0.35])
    ax.annotate("boundary scan",
                xy=tip_scan, xytext=lbl_scan,
                fontsize=9, color=BLUE, ha="left",
                arrowprops=dict(arrowstyle="-", color=BLUE, lw=0.5))

    # -- Chord from scan-end to exit (dashed blue) --
    ax.plot([scan_end[0], exit_xy[0]], [scan_end[1], exit_xy[1]],
            color=BLUE, lw=1.2, linestyle="--")
    mid_chord = 0.5 * (scan_end + exit_xy)
    # Place label inside the disk, to the upper-left of the chord.
    lbl_chord = np.array([-0.65, 1.05])
    ax.annotate("chord to exit\n(after broadcast)",
                xy=mid_chord, xytext=lbl_chord,
                fontsize=9, color=BLUE, ha="center",
                arrowprops=dict(arrowstyle="-", color=BLUE, lw=0.5))

    # -- Exit star + label pushed outside the disk to the upper-right --
    ax.scatter(*exit_xy, s=120, marker="*", color=RED, zorder=7)
    lbl_exit = np.array([exit_xy[0] - 0.10, exit_xy[1] + 0.35])
    ax.annotate(r"exit at unknown $\theta$",
                xy=exit_xy, xytext=lbl_exit,
                fontsize=10, color=RED, ha="center",
                arrowprops=dict(arrowstyle="-", color=RED, lw=0.7))

    # -- Axes / cleanup --
    ax.set_xlim(-1.45, 1.90)
    ax.set_ylim(-1.45, 1.70)
    ax.set_aspect("equal")
    ax.axis("off")

    out = FIG_DIR / "problem_setup.pdf"
    plt.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
