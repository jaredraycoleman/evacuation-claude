"""Figure: the k=3 wireless evacuation problem setup.

Produces figures/problem_setup.pdf showing:
    - the unit disk with centre marked,
    - three robots stacked at the centre,
    - an "unknown exit" marked with a question mark on the boundary,
    - a generic robot trajectory (to boundary + scan + chord).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
FIG_DIR.mkdir(exist_ok=True)


def main() -> None:
    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    circle_t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(circle_t), np.sin(circle_t), color="black", lw=1.3)

    # Centre + robots as a stack of tiny markers.
    ax.scatter([0], [0], s=16, color="black", zorder=5)
    for i, dx in enumerate([-0.06, 0.0, 0.06]):
        ax.annotate(f"$r_{i+1}$", (dx, -0.13), ha="center", fontsize=9)

    # Unknown exit at some theta (example angle).
    theta_ex = 1.1
    ex_xy = np.array([np.cos(theta_ex), np.sin(theta_ex)])
    ax.scatter(*ex_xy, s=60, marker="*", color="crimson", zorder=6)
    ax.annotate(r"exit at unknown $\theta$",
                xy=ex_xy, xytext=(ex_xy[0] + 0.08, ex_xy[1] + 0.14),
                fontsize=9, color="crimson",
                arrowprops=dict(arrowstyle="-", color="crimson", lw=0.7))

    # Sample trajectory for one robot (deep grey):
    # radial approach to angle alpha, then short scan CCW, then chord to exit.
    alpha = 0.55
    ax.plot([0, np.cos(alpha)], [0, np.sin(alpha)],
            color="#3a7ab0", lw=1.4)
    ax.annotate("radial\napproach",
                (0.5 * np.cos(alpha) - 0.2, 0.5 * np.sin(alpha) + 0.05),
                fontsize=8, color="#3a7ab0")

    scan_t = np.linspace(alpha, alpha + 0.55, 80)
    ax.plot(np.cos(scan_t), np.sin(scan_t), color="#3a7ab0", lw=1.8)
    ax.annotate("boundary\nscan", (np.cos(alpha + 0.55) + 0.05,
                                   np.sin(alpha + 0.55) + 0.08),
                fontsize=8, color="#3a7ab0")

    # Chord from scan endpoint to exit (dashed):
    end_xy = np.array([np.cos(alpha + 0.55), np.sin(alpha + 0.55)])
    ax.plot([end_xy[0], ex_xy[0]], [end_xy[1], ex_xy[1]],
            color="#3a7ab0", lw=1.0, linestyle="--")
    ax.annotate("straight chord\nafter broadcast",
                ((end_xy[0] + ex_xy[0]) / 2 + 0.03,
                 (end_xy[1] + ex_xy[1]) / 2 + 0.05),
                fontsize=8, color="#3a7ab0")

    ax.set_xlim(-1.25, 1.35)
    ax.set_ylim(-1.25, 1.25)
    ax.set_aspect("equal")
    ax.axis("off")

    out = FIG_DIR / "problem_setup.pdf"
    plt.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
