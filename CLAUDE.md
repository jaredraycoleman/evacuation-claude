# CLAUDE.md — Project Hub

Working directory for exploring open problems in the **robot evacuation**
literature. Claude keeps this file updated as the single source of truth:
current results, what's running, what's next.

## Domain in one paragraph

Robots start at the center of a unit disk; an *exit* sits at an unknown point
on the perimeter; the goal is a trajectory that minimizes the worst-case time
until all (or a designated) robot reaches the exit. Variants differ in
**number of robots**, **communication model** (F2F / wireless / sender-only
/ receiver-only), and **objective** (all-evacuation vs priority/queen).

## Repository map

```
papers/         the four reference papers (PDF)
notes/          Claude's literature notes — per-paper summary in papers.md
experiments/    Python (uv) — simulators, numerical optimization, plots
analysis/       Wolfram Language scripts for symbolic / analytical work
figures/        generated plots (referenced by paper/main.tex)
paper/          LaTeX write-up; build with paper/build.sh
```

## Workflow (Claude reads this every session)

1. **Reload context.** Read this file top-to-bottom, then skim
   `notes/papers.md`. If a known open problem is in progress, reread its
   section below before changing anything.
2. **Pick or continue an open problem.** Either continue the one listed
   under *Currently working on*, or — if that's stuck / closed — survey
   the *Open problem candidates* list and pick the most tractable one.
   Write one paragraph in this file justifying the choice.
3. **Do the work.**
   - Numerical / simulation / plotting → `experiments/` via **uv**:
     - `uv sync` once per environment
     - `uv run python -m experiments.<module>` to execute
     - Save plots into `figures/` as PDF (LaTeX-friendly).
   - Analytical / symbolic → **wolframscript** under `analysis/` (`.wls`).
     If `wolframscript` is unavailable, use `sympy` inside `experiments/`
     and note the substitution below under *Environment caveats*.
   - Write results up in `paper/main.tex`; compile with `paper/build.sh`
     and confirm the PDF builds cleanly before declaring progress.
4. **Update this file before ending the session.** Move entries between
   sections (e.g., *Currently working on* → *Results* or *Dead ends*).
   Keep *Next steps* concrete and actionable — a single sentence each.

## Tooling

- **Python**: `uv` (managed via `pyproject.toml`). Never `pip install` into
  the system interpreter.
- **Analytical**: `wolframscript -file analysis/<file>.wls`. Fallback:
  `sympy`.
- **LaTeX**: `pdflatex` via `paper/build.sh` (two passes for refs).
- **Figures**: save as vector PDF into `figures/` and include from
  `paper/main.tex` with `\includegraphics`.

## Environment caveats

- At session start, smoke-test both:
  - `wolframscript -code '1+1'` — prints `2` if and only if the kernel is
    licensed *and* on PATH. `-version` alone is not enough: the binary can
    report its version while the kernel still refuses to run scripts.
  - `pdflatex --version`.
- `wolframscript` lives at
  `/usr/local/Wolfram/Wolfram/14.3/Executables/wolframscript` with a
  user-space symlink in `~/.local/bin/wolframscript`.
- If the smoke test fails with an *activation* error, ask the user to run
  `wolframscript -activate` in their own terminal (it's interactive — needs
  Wolfram ID + password). Don't try to work around it.
- If `wolframscript` is genuinely unusable, fall back to `sympy` inside
  `experiments/` and flag any result that really needs Wolfram
  (non-elementary roots, high-precision `FindRoot`, etc.) so it can be
  revisited.

## Currently working on

**3-robot wireless evacuation on the disk — close the $O(k^{-4/3})$ gap at $k=3$.**
Czyzowicz et al. give UB $3 + \pi/k + O(k^{-4/3})$ and matching-form LB $3 + \pi/k$;
at $k=3$ the LB is $3 + \pi/3 \approx 4.0472$ and the UB carries a small additive
slack from a generic-$k$ construction that is unlikely to be tight for a specific
small $k$. We picked this because (a) the wireless model has a simple simulator
— once any robot finds the exit, all robots head straight for it, so evaluating
a candidate trajectory is 1D optimization over the exit angle; (b) the gap is
concrete and small, so a tailored 3-robot deployment pattern has a realistic
shot at closing or shrinking it; (c) the tooling (parametrized trajectory +
worst-exit evaluator + global optimizer) is reusable for the 2-robot F2F UB
problem and the priority-evacuation small-$n$ cases later. Deliverable target:
either a closed-form improved UB matching $3 + \pi/3$, or a numerical UB below
the Czyzowicz et al. bound with a cleanly-parametrized trajectory.

**Progress so far:**
- Env green: `wolframscript`, `pdflatex`, `uv sync`, `paper/build.sh` all work.
- `experiments/wireless_k3.py` is a generic wireless-$k$ simulator: pass $k$
  trajectories $r_i: [0,T] \to \mathbb{R}^2$ (unit-speed, start at origin), get
  worst-case evac time over exit angle $\theta$.
- Baseline verified: naive arc-partition at $k=3$ gives worst-case $\approx 4.8231$
  (analytical $1 + 2\pi/3 + \sqrt{3} \approx 4.8264$; diff is grid tolerance).
  Worst exit at $\theta \approx 2\pi/3$ (end of robot 0's arc); bottleneck is
  the robot diametrically opposite at chord distance $\sqrt{3}$.

**Next action:** Parametrize a richer trajectory family and numerically
optimize. Candidate ideas (in order of complexity):
  1. **Staggered arc-partition** — robots scan overlapping arcs, so at any
     moment of discovery the non-discovering robots are already partway
     through their own approach.
  2. **Chord-shortcut on discovery** — the generic algorithm assumes robots
     are still scanning when any broadcast arrives; tune arc lengths so the
     bottleneck chord is balanced against scan time.
  3. **Non-radial deployment** — approach boundary along a curve that keeps
     robots closer to the likely exit region, trading deployment time for
     smaller worst-case chord.
First concrete task: build a parametrized family (arc midpoint angles +
arc lengths) and run a Nelder-Mead / differential-evolution search against
the worst-case evaluator.

## Open problem candidates

Picked from the four reference papers. Tractability is a rough guess; revise
as we learn more.

1. **Lower bound for 2-robot F2F evacuation on the disk.** Upper bound is
   5.625 (Brandt et al.); best known lower bound is weaker. Can the
   worst-case-exit characterization (non-differentiability or
   $2\cos\beta + \cos\gamma = 1$) be leveraged to push the lower bound up?
   *Tractability:* medium — heavy numerics + geometric argument.
2. **Exact constants for small $k$ in the $k$-robot non-wireless model.**
   Czyzowicz et al. give $3 + 2\pi/k - O(k^{-2})$. For $k=2,3$, can the
   $O(k^{-2})$ slack be closed?
   *Tractability:* medium — algorithm design + matching lower bound.
3. **Tight bound in the Sender/Receiver (SR) asymmetric model.**
   Georgiou-Giachoudis-Kranakis show the SR model sits strictly between F2F
   and WiFi, with upper bound $< \pi + 2$. Is there a matching lower bound,
   or a smaller upper bound?
   *Tractability:* medium-hard.
4. **Priority evacuation with few servants.** Closed-form for $n = 1, 2, 3$?
   The general bound $2 + 4(\sqrt 2 - 1)\pi/n$ is stated for $n \ge 4$.
   *Tractability:* low-medium for $n=1$; harder for 2–3.
5. **Non-disk room shapes.** Brandt et al. remark their worst-case-exit
   tool applies to any shape. Pick a shape (square? regular $n$-gon?) and
   carry a known algorithm over; what's the worst-case time?
   *Tractability:* open-ended; good for numerics-first exploration.

## Results

*(populated as theorems / numerical bounds land)*

## Dead ends

*(avoid re-walking these — record the approach tried and why it failed)*

## Next steps

- Pick one candidate from the list above and move it into *Currently working
  on*.
- Set up `uv sync` on first run and confirm `paper/build.sh` produces
  `paper/main.pdf`.
- If `wolframscript` is needed for the chosen problem, ask the user to
  install it (or confirm Sympy-only is fine).

## Conventions

- All numerical bounds in this document are **worst-case evacuation time in
  units of disk-radius / unit-speed** unless noted otherwise.
- When recording a new result here, also add a theorem/proposition to
  `paper/main.tex` so the write-up stays in sync.
- When a numerical experiment produces a figure, commit the generating
  script under `experiments/` and the PDF under `figures/` with the same
  stem (e.g., `fig_worst_exit_scan.py` → `fig_worst_exit_scan.pdf`).
