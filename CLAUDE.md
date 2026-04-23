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
site/           static website (vanilla HTML/CSS/JS) deployed to Cloudflare
                Workers; interactive A_3 visualiser + embedded report PDF
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
4. **Keep the static site current.** If the session produced new
   results, updated `paper/main.tex`, added a figure, or otherwise
   changed user-visible state, also:
   - Update `site/index.html` (problem description / interactive) if the
     wording or the JS port in `site/main.js` is now stale relative to
     the paper.
   - Run `./site/build.sh` (rebuilds the paper, copies the PDF into
     `site/report.pdf`, generates `site/build-info.js` with the commit
     SHA + dates).
   - Deploy: `npx wrangler deploy` — this publishes to
     `https://evacuation-claude.jaredraycoleman.workers.dev`. Cloudflare
     auth is already configured (API token in env, or `wrangler login`
     credentials at `~/.config/.wrangler/`); if `whoami` fails, ask the
     user for a new token rather than guessing.
   - If the JS simulator in `site/main.js` now disagrees with the Python
     simulator numerically, fix it before deploying — the site is
     promoted as matching the paper.
5. **Update this file before ending the session.** Move entries between
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
- **Site**: `./site/build.sh` rebuilds the paper, copies the PDF to
  `site/report.pdf`, and generates `site/build-info.js`. `npx wrangler
  deploy` ships `site/` as an assets-only Cloudflare Worker; config is
  in `wrangler.toml` at the repo root. Live URL:
  `https://evacuation-claude.jaredraycoleman.workers.dev`.

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

## Results to date

**Problem 1 (wireless k=3).** Proved UB $\approx 4.21852$ via rebalanced
$A_3$, improving Czyzowicz et al.~by $\approx 1/1280$. Full zone proof
in paper.

**Problem 2 (priority n=1 wireless).** Proved bracket $1+\pi \le
T^* \le 1 + 2\pi/3 + \sqrt 3$, i.e., $4.1416 \le T^* \le 4.8264$, gap
$\approx 0.685$. UB via 2-robot-wireless reduction; LB via Lemma 5
half-chord adversary. The specific $n=1$ case is explicitly
unaddressed in~\cite{czyzowicz2018priority}.

## Currently working on

**3-robot wireless evacuation, $k=3$ — tighten the lower bound.**
With the UB side closed as much as the piecewise-linear class allows (Problem 1,
4.21852 proved rigorously), the remaining gap $\approx 0.06$ is entirely on the
LB side. CGK+14's Theorem 7 gives $4.15937$; our UB is $4.21852$. Goal: prove
an LB strictly bigger than $4.15937$, ideally closing or measurably shrinking
the gap.

**Why this is tractable as the next step.** The slack in CGK's argument is
localized (see *Dead ends / Negative results* for the analysis): their Lemma 6
bounds the finder's chord at $2\sin(xk/2) \approx 1.886$ at their optimal
$x \approx 1.274$. In our UB algorithm at its worst-case $\theta \approx 1.486$,
the *bottleneck* (not the finder) has chord $\sqrt 3 \approx 1.732$ — larger
than what Lemma 5 guarantees at the corresponding $x$. If we can lift the
argument from "finder chord" to "bottleneck chord" universally, the LB should
improve in the $x \approx 1.4$–$1.5$ region, where CGK's bound is already
weakening as it chases one decreasing sine.

**Lines of attack:**
1. *Two-point adversary with nearest-robot distance.* The first-visit time to
   either of two candidate exits is $\ge 1 + x + \min_i \|\rho_i(1+x) - \{p,q\}\|$.
   If we can LB the nearest-robot distance away from zero (adversary chooses
   $p,q$ interior to the unexplored arc, not at its edges), we tighten
   $t_{\text{find}}$.
2. *Multi-exit adversary.* Keep three candidate exits mutually far apart;
   force some robot to be far from whichever one the adversary eventually
   commits to.
3. *Bottleneck-chord lemma.* Strengthen Lemma 6 to argue that after the
   finder reaches one point, a second robot (that was near the sister
   candidate) must have chord $\ge 2\sin(xk/2)$ to the exit at discovery
   time — not via a new quantity, but by identifying which robot is the
   bottleneck.

**Tooling plan.** Extend `analysis/czyzowicz_lb.wls` to a generic framework
`analysis/lb_bottleneck.wls` that: parametrizes 2/3-exit adversaries, encodes
the algorithm-side invariants (scan-budget, chord-mover budget), and evaluates
the LB symbolically. If a proof is elusive, do a numerical minimax — over
algorithm-side parameters — to see what the best achievable LB looks like;
that indicates the ceiling, even if the proof isn't yet clean.

**Progress so far (this problem):**

**New universal inequality (derived + numerically validated).** For any
3-robot wireless algorithm, any $x > 0$, any $p$ unexplored at time $1+x$:
$$\text{evac}(p) \ge (1+x) + \max_j \|R_j(1+x) - p\|.$$
Triangle-inequality proof is clean; see `notes/lb-refinement.md`.

**On $A_3(y_{\text{opt}})$ the refined inequality is nearly tight.**
Applied to the UB algorithm directly
(`experiments/lb_refined_a3.py`): $\sup_{x,p} [(1+x) + \max_j] = 4.21836$
at $x \approx 1.5$, within $0.0002$ of the actual UB $4.21852$. The
inequality therefore captures almost all of the slack in $A_3$'s own
worst case — encouraging sign that this is "the right object."

**Refined alone is not universally stronger than CGK**
(`experiments/lb_minimax.py`). Raw min over algorithm states of
$\sup_{x,p}$ = 3.95 at $x = 1.95$, achieved by "all robots at origin"
states. These states are degenerate (alg never evacuates), but they
defeat the refined bound in isolation.

**Combined bound $\max(\text{state-tight CGK}, \text{refined})$ looks
promising** (`experiments/lb_combined.py`). For any alg state at
$x^* = 1.274$ (CGK's optimum):
- If $|E(2.274)| < 3x^*$: state-tight CGK is $> 4.159$ by itself.
- If $|E| = 3x^*$ exactly (CGK's universal bound tight): all 3 robots
  are on the boundary, and numerically min refined chord over that
  class is $\ge 1.898$, giving refined LB $\ge 4.172 > 4.159$.

So the two cases together rule out a universal LB of exactly $4.159$.
The boundary case (mixing small $|E| < 3x$ with small refined) needs
more care: the combined LB there depends on how fast min-refined
decreases as $|E|$ drops below $3x$. Numerical slope suggests a
$0.01$ improvement is defensible, not $0.06$.

**Blocking gap.** A clean rigorous LB improvement requires proving:
\emph{for any 3 robots on $\partial D$ with $|E| = 3x$ (3 disjoint
scan arcs of length $x$ each ending at robot positions),
$\max_{p \in E^c} \max_j \|R_j - p\| \ge 2\cos(x/4)$}, achieved by
the 3-fold symmetric arrangement. This appears true from two explicit
cases (3-fold symmetric and single-arc), but no proof yet.

Even that would only give $\approx 0.013$ improvement at $x_{CGK}^*$.
Closing the full $0.06$ gap is substantially harder.

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

**Problem 1 (k=3 wireless UB) — partially closed.** $T^*_{3,\text{wireless}} \le 4.21852$,
improving CGK+14's 4.21930 by $\approx 1/1280$. Rigorous zone-decomposition proof +
transcendental characterisation of $y_{\text{opt}}$. Tooling: `experiments/a3_*.py`,
`analysis/a3_*.wls`, Theorem 1–3 in `paper/main.tex`.

Under the piecewise-linear + boundary-scan + chord trajectory class, this appears
to be a local minimum: families C/D/E/F (see `experiments/family_*.py`) all
return to the $A_3(y_{\text{opt}})$ point under local search. Further UB
improvement would require curved (non-radial) ingress — open.

## Dead ends

**UB Family A** (`experiments/wireless_k3_opt.py`): single CCW/CW scan + wait at
origin. Stalls at naive $1 + 2\pi/3 + \sqrt 3 \approx 4.826$. Misses $A_3$'s
return-and-redeploy move.

**UB Family B** (`experiments/wireless_k3_opt2.py`): two-phase scan for robot 2.
DE search stalls at 4.826 across all 8 sign combos — the good region is a
narrow manifold DE doesn't find from Sobol init.

**UB Families C/D/E/F** (`experiments/family_*.py`): richer trajectory classes
(extra scan after redeploy, 3-fold symmetric redeploy, dual redeploy, direct
chord). All return to the $A_3(y_{\text{opt}})$ point under local search.

**Brandt F2F $\varepsilon$-cut** (`experiments/brandt_f2f.py`): probe of a small
off-path cut $(d_2, y_2)$ added to Brandt et al.'s symmetric 2-cut algorithm
$A(y, \alpha, d)$. At the Brandt optimum $d_2 = 0$ is a local minimum —
the symmetric 2-cut family is locally tight.

## Next steps

- Encode CGK's Lemma 5/6 carefully in Wolfram, identifying each inequality
  in the proof chain so we can see exactly where to tighten.
- Attempt the two-point adversary with a nearest-robot-distance term; check
  numerically whether any $x$-region has improved LB beyond 4.15937.
- If numerical minimax over "algorithms" (parameterised by scan-budget
  distribution at time $1+x$) shows the true LB is $\le 4.16$, this is
  a dead end; otherwise formalise what the numerics suggest.

## Conventions

- All numerical bounds in this document are **worst-case evacuation time in
  units of disk-radius / unit-speed** unless noted otherwise.
- When recording a new result here, also add a theorem/proposition to
  `paper/main.tex` so the write-up stays in sync.
- When a numerical experiment produces a figure, commit the generating
  script under `experiments/` and the PDF under `figures/` with the same
  stem (e.g., `fig_worst_exit_scan.py` → `fig_worst_exit_scan.pdf`).
