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
- `experiments/wireless_k3.py`: generic wireless-$k$ simulator (unit-speed
  trajectories in, worst-case evac time out).
- Baseline: naive arc-partition at $k=3$ gives $1 + 2\pi/3 + \sqrt{3} \approx 4.826$.
- `experiments/wireless_k3_opt.py`: differential-evolution search over
  Family A = {radial approach + single CCW/CW scan + wait at origin}.
  Converges to 4.826 — no improvement. Family A is not rich enough.

**Key correction from reading `papers/ev-disc-k-robots.pdf` (Czyzowicz et al.):**
- Asymptotic UB $3 + \pi/k + O(k^{-4/3})$ (Theorem 10) only applies for $k \ge 100$.
- For $k=3$ they have a *separate* algorithm A3 (Theorem 6) with
  UB $= \frac{4\pi}{9} + \frac{2\sqrt{3}+5}{3} + \frac{1}{600} \approx 4.2193$
  and matching LB $\approx 4.159$ (Theorem 7). **Gap $\approx 0.06$.**
- A3 uses a "return-to-center + redeploy" move ($r_3$ scans for time $y$,
  walks back to origin, walks out again at angle $\pi - y/2$). This move is
  not in Family A — explains why numerical search stalled at 4.826.

**Revised target:** beat 4.2193 at $k=3$ with a proved bound. The `1/600`
term in their UB is a hand-tuned slop, strongly suggesting it can be
tightened with careful parameter optimization.

**Family B attempted (`experiments/wireless_k3_opt2.py`):** Robot 2 gets a
two-phase trajectory (approach + scan + return-to-origin + approach +
scan). Result: blind DE search still converges to $\approx 4.826$ across
all 8 discrete sign combos. Diagnosis: when robot 2's first scan
$L_{2a}$ is large, its second phase doesn't start until after $t_{\text{found}}$,
so the new move is effectively inactive; the optimizer finds a Family-A
local minimum and stays there. Blind global search in this space is
wasteful — the good region is a narrow manifold (small $L_{2a}$, robot 2
acts as a "positioner" not a scanner) that DE with Sobol init misses.

**Next action:** Implement A3 **literally** from the paper (exact parameters
of Czyzowicz et al. Theorem 6): $r_1$ to $A$; $r_2, r_3$ together to $B$
with angular offset $y = 4\pi/9 + 2\sqrt{3}/3 - 401/300$; scan directions
as specified; $r_3$ returns to origin and redeploys at angle $\pi - y/2$
CW of $RB$. Verify the simulator reproduces $\approx 4.2193$. Then: local
search (e.g., Nelder-Mead seeded at A3's $(y, \text{offsets})$) to see if
the `1/300`-scale hand-tuning can be tightened, and extract the balancing
conditions for a Wolfram-symbolic proof of whatever optimum we land at.

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
