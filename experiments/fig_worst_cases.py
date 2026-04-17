"""Figure: the four simultaneous worst-case exit angles at y = y_opt.

Produces figures/worst_cases.pdf: 2x2 grid of snapshots, each showing the
unit disk at the discovery time for one of the four binding worst-case
angles, with robot positions and the bottleneck chord drawn.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from experiments.wireless_k3_opt2 import (
    _hit_times_phases,
    _positions_phases,
)
from experiments.a3_optimize import a3_robots_at

FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
TWO_PI = 2.0 * np.pi
Y_OPT = 1.2158578321292429

# From the diagnosis run: four binding worst-case angles.
WORST_THETAS = [1.486472, 2.330820, 2.736509, 3.580857]

ROBOT_COLORS = {0: "#2a9d8f", 1: "#e76f51", 2: "#264653"}
ROBOT_LABELS = {0: r"$r_1$", 1: r"$r_2$", 2: r"$r_3$"}


def snapshot(ax, theta_exit: float) -> None:
    robots = a3_robots_at(Y_OPT)

    # Figure out t_found for this theta.
    t_grid = np.array([theta_exit])    # single theta
    ht = np.stack([_hit_times_phases(r, t_grid) for r in robots], axis=0)
    t_found = float(ht.min(axis=0)[0])
    t_arr = np.array([t_found])
    positions = np.stack([_positions_phases(r, t_arr)[0] for r in robots],
                         axis=0)    # shape (3, 2)
    exit_xy = np.array([np.cos(theta_exit), np.sin(theta_exit)])
    chords = np.linalg.norm(positions - exit_xy[None, :], axis=-1)
    bottleneck = int(np.argmax(chords))
    evac = t_found + float(chords[bottleneck])

    # Disk
    tc = np.linspace(0, TWO_PI, 400)
    ax.plot(np.cos(tc), np.sin(tc), color="black", lw=1.0)

    # Origin
    ax.scatter([0], [0], s=12, color="black", zorder=4)

    # Robots at t_found
    for i in range(3):
        ax.scatter(*positions[i], s=40, color=ROBOT_COLORS[i], zorder=5)
        ax.annotate(ROBOT_LABELS[i],
                    positions[i] + np.array([0.06, 0.06]),
                    color=ROBOT_COLORS[i], fontsize=9)

    # Exit
    ax.scatter(*exit_xy, s=90, marker="*", color="crimson", zorder=6)

    # Bottleneck chord
    ax.plot([positions[bottleneck, 0], exit_xy[0]],
            [positions[bottleneck, 1], exit_xy[1]],
            color=ROBOT_COLORS[bottleneck], lw=1.5, linestyle="--")

    # Titles
    lbl = f"$\\theta$={theta_exit:.3f}, found by $r_{{{int(np.argmin(ht[:,0]))+1}}}$\n"
    lbl += (f"bottleneck {ROBOT_LABELS[bottleneck]}, "
            f"chord $={chords[bottleneck]:.3f}$, "
            f"evac $={evac:.4f}$")
    ax.set_title(lbl, fontsize=9)

    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.set_aspect("equal")
    ax.axis("off")


def main() -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.6, 7.6))
    for ax, theta in zip(axes.flat, WORST_THETAS):
        snapshot(ax, theta)

    fig.suptitle(
        (r"Four simultaneous worst-case angles at $y_{\rm opt}\approx 1.2159$;"
         "\n"
         r"all give evacuation time $\approx 4.21852$."),
        fontsize=10,
    )

    out = FIG_DIR / "worst_cases.pdf"
    plt.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
