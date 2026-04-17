"""Figure: the A_3 algorithm trajectories (with y = y_opt).

Produces figures/a3_algorithm.pdf showing:
    - unit disk,
    - points A (angle 0) and B (angle 2 pi - y),
    - r_1 trajectory: origin -> A + scan CCW by pi - y/2,
    - r_2 trajectory: origin -> B + scan CW by pi - y/2,
    - r_3 trajectory: origin -> B + scan CCW by y + back to origin + out
      to angle pi - y/2.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
TWO_PI = 2.0 * np.pi

Y_OPT = 1.2158578321292429


def boundary_arc(a_start: float, a_end: float, n: int = 200):
    t = np.linspace(a_start, a_end, n)
    return np.cos(t), np.sin(t)


def main() -> None:
    fig, ax = plt.subplots(figsize=(4.6, 4.6))

    y = Y_OPT
    L = np.pi - y / 2.0
    aA = 0.0
    aB = TWO_PI - y
    a_red = np.pi - y / 2.0

    # Disk outline.
    tc = np.linspace(0, TWO_PI, 400)
    ax.plot(np.cos(tc), np.sin(tc), color="black", lw=1.2)

    # Origin
    ax.scatter([0], [0], s=20, color="black", zorder=5)
    ax.annotate("origin", (0.02, -0.12), fontsize=8)

    # A, B
    A_xy = np.array([np.cos(aA), np.sin(aA)])
    B_xy = np.array([np.cos(aB), np.sin(aB)])
    ax.scatter(*A_xy, s=25, color="black", zorder=6)
    ax.annotate("$A$", A_xy + np.array([0.05, -0.02]), fontsize=10)
    ax.scatter(*B_xy, s=25, color="black", zorder=6)
    ax.annotate("$B$", B_xy + np.array([0.05, -0.05]), fontsize=10)

    # Redeploy target for r_3
    red_xy = np.array([np.cos(a_red), np.sin(a_red)])
    ax.scatter(*red_xy, s=25, color="black", zorder=6)
    ax.annotate("$C = (\\cos(\\pi - y/2),\\ \\sin(\\pi - y/2))$",
                red_xy + np.array([-0.23, 0.05]), fontsize=8)

    # r_1: origin -> A (radial)
    ax.plot([0, A_xy[0]], [0, A_xy[1]],
            color="#2a9d8f", lw=1.6, label=r"$r_1$")
    # r_1 scan CCW from 0 by L.
    xs, ys = boundary_arc(0, L, 200)
    ax.plot(xs, ys, color="#2a9d8f", lw=2.2)

    # r_2: origin -> B (radial)
    ax.plot([0, B_xy[0]], [0, B_xy[1]],
            color="#e76f51", lw=1.6, label=r"$r_2$")
    # r_2 scan CW from aB by L  => arc [aB - L, aB].
    xs, ys = boundary_arc(aB - L, aB, 200)
    ax.plot(xs, ys, color="#e76f51", lw=2.2)

    # r_3: origin -> B (radial; overlaps r_2 so draw slightly offset / dashed)
    ax.plot([0, B_xy[0]], [0, B_xy[1]],
            color="#264653", lw=1.6, linestyle=(0, (3, 2)),
            label=r"$r_3$")
    # r_3 scan CCW from aB by y  => covers aB .. aB + y = 2pi = 0 (meeting A).
    xs, ys = boundary_arc(aB, aB + y, 200)
    ax.plot(xs, ys, color="#264653", lw=2.2, linestyle=(0, (3, 2)))

    # r_3 return to origin: straight line from (1, 0) to (0, 0). Scan ended at
    # angle aB + y == 2 pi which is the point (1, 0).
    end_scan3 = np.array([1.0, 0.0])
    ax.plot([end_scan3[0], 0], [end_scan3[1], 0],
            color="#264653", lw=1.3, linestyle=(0, (3, 2)))

    # r_3 redeploy: origin -> C.
    ax.plot([0, red_xy[0]], [0, red_xy[1]],
            color="#264653", lw=1.3, linestyle=(0, (3, 2)))

    # Label the y-offset angle near the origin with a little arc.
    arc_r = 0.18
    tt = np.linspace(aA, aB, 80)   # going CCW from A to B the long way, but
    # we want the SHORT CW arc, so go from 0 to -y:
    tt = np.linspace(0, -y, 80)
    ax.plot(arc_r * np.cos(tt), arc_r * np.sin(tt),
            color="black", lw=0.8)
    ax.annotate("$y$",
                (arc_r * np.cos(-y / 2) + 0.03, arc_r * np.sin(-y / 2) - 0.06),
                fontsize=9)

    ax.legend(loc="upper left", fontsize=9, framealpha=0.95)
    ax.set_xlim(-1.25, 1.35)
    ax.set_ylim(-1.25, 1.25)
    ax.set_aspect("equal")
    ax.axis("off")

    title = (r"$A_3(y_{\rm opt})$: $r_1$ to $A$, $r_2$ and $r_3$ to $B$;"
             "\n"
             r"$r_3$ scans $y$ then returns to the origin and redeploys to $C$.")
    ax.set_title(title, fontsize=9)

    out = FIG_DIR / "a3_algorithm.pdf"
    plt.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
