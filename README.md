# evacuation-claude

An experiment in **LLM-assisted mathematical research**: can a large language
model, used as a sustained research collaborator, close open gaps in the
robot-evacuation literature?

## The setup

Robot evacuation from the disk is a well-studied theoretical-CS problem:
$k$ unit-speed robots start at the centre of a unit disk, an adversary
places an exit at an unknown point on the boundary, and the robots want
to minimise the worst-case time for all (or a designated) robot to
reach the exit. Variants differ in the number of robots, the
communication model (face-to-face / wireless / sender-only /
receiver-only), and the objective (all-evacuation vs. priority/queen).

Four reference papers sit in `papers/`:

- **CGK+14** — Czyzowicz, Gąsieniec, Gorry, Kranakis, Martin, Pająk.
  *Evacuating Robots from an Unknown Exit in a Disk.* ($k$-robot
  wireless / non-wireless.)
- **BLLSW17** — Brandt, Laufenberg, Lv, Stolz, Wattenhofer.
  *Collaboration without Communication.* (2-robot F2F.)
- **GGK22** — Georgiou, Giachoudis, Kranakis.
  *Evacuation from a Disk for Robots with Asymmetric Communication.*
  (Sender/Receiver.)
- **CEGKLPR18** — Czyzowicz et al. *Priority Evacuation from a Disk
  Using Mobile Robots.* (Queen + $n$ servants.)

Each paper leaves explicit gaps: looser upper or lower bounds,
unaddressed small cases, unknown exact constants. The project drove
Claude (Anthropic's LLM) through those gaps systematically: read the
paper, pick the most tractable hole, build the tooling (numerical
simulators in `experiments/`, symbolic work in `analysis/`), and push
until the gap either closed, split, or revealed a fundamental barrier.

## Status by problem

| # | Problem | Reference | Prior state | Our state | Outcome |
|---|---|---|---|---|---|
| 1 | 3-robot wireless, $k=3$, UB | CGK+14 Thm 6 | $4.21930$ | $\mathbf{4.21852}$ | **Gap tightened** by $\sim 1/1280$, with a rigorous zone-decomposition proof and a transcendental characterisation of the optimal deployment parameter. Flagship result of the paper. |
| 2 | 3-robot wireless, $k=3$, LB | CGK+14 Thm 7 | $4.15937$ | unchanged | Identified where the published argument is slack (it bounds the *finder*'s chord, but in our UB algorithm the *bottleneck* is a different robot). A formal refinement is open. |
| 3 | 2-robot F2F, Brandt "$\varepsilon$-cut" probe | BLLSW17 | UB $5.625$ | unchanged | **Null result.** Probe of an $\varepsilon$-cut refinement of Brandt's symmetric-cut family converged back to the original. Retained as `experiments/brandt_f2f.py`. |
| 4 | Priority evacuation, $n=1$ servant, wireless | CEGKLPR18 (unaddressed for $n=1$) | — | bracket $1 + \pi \le T^* \le 1 + 2\pi/3 + \sqrt 3$ | Exploratory bracket ($\approx 4.14 \le T^* \le 4.83$, gap $\approx 0.685$). Too loose to include in the paper; retained as `experiments/priority_n1*.py` and `analysis/priority_n1.wls`. |
| 5 | Priority $n=2$ / SR model | CEGKLPR18, GGK22 | — | no result | Scoped, drafted, and cut from the paper. |

**Headline.** One gap measurably tightened with a rigorous proof; one
probe that confirmed a published UB is tight within its natural family;
several negative or exploratory brackets.

## Next directions

The most tractable remaining frontier is the $\approx 0.060$ gap between
our UB ($4.21852$) and CGK+14's LB ($4.15937$) at $k=3$.

- **Lower-bound refinement.** CGK+14's Lemma 6 bounds the *finder*'s
  chord to the exit. In $A_3(y_{\mathrm{opt}})$ the finder has chord
  $0$; the bottleneck is a different robot with chord $\sqrt 3$. A
  multi-exit adversary keeping three candidate exits mutually far apart
  should force *some* robot (not just the finder) far from the
  eventual exit. Formalising this for $k = 3$ is open.
- **Curved / non-radial approaches.** Every trajectory family we
  explored (A, C, D, E, F) uses a radial ingress to the boundary.
  Replacing the radial segment with a tangential or logarithmic-spiral
  arc might unlock configurations that aren't reachable in the
  piecewise-linear + boundary-scan + chord class. This is where any
  remaining UB improvement likely hides, if any.
- **Brandt F2F beyond $\varepsilon$-cut.** The null result above
  suggests the symmetric 2-cut family is a local minimum. Asymmetric
  cuts, agreement-point families, or curved cuts are natural next
  probes.
- **Priority $n=1$.** The $4.14$ / $4.83$ bracket is embarrassingly
  wide. Either a sharper servant-schedule UB or a geometric LB stronger
  than the trivial $1 + \pi$ would likely improve both ends; a tighter
  $n=1$ result would then serve as the base case for a small-$n$
  sequence.
- **Non-disk shapes.** BLLSW17 explicitly notes their worst-case-exit
  characterisation applies to any convex room. Porting the flagship
  $A_3$ analysis to a square or regular polygon is open-ended but a
  natural numerics-first avenue.

## Layout

- `paper/main.tex` — write-up of result #1; build with `paper/build.sh`
  (runs pdflatex + bibtex + pdflatex twice).
- `experiments/` — Python (via `uv`) simulators and optimisers. The
  `a3_*.py` files drive the flagship result; `family_{c,d,e,f}.py`
  rule out natural enrichments; `brandt_f2f.py` is the F2F probe;
  `priority_n1*.py` is the servant bracket.
- `analysis/` — Wolfram Language scripts for symbolic steps
  (`a3_zone_proof.wls`, `a3_balance.wls`, `czyzowicz_lb.wls`,
  `priority_n1.wls`).
- `notes/papers.md` — per-paper summaries used to pick gaps.
- `figures/` — generated PDF figures included by the paper.
- `site/` — static site with an interactive $A_3(y)$ visualiser and the
  embedded report. Deployed at
  <https://evacuation-claude.jaredraycoleman.workers.dev>.
- `CLAUDE.md` — the live research log the LLM reads and updates each
  session.

## Tooling

- **Python:** `uv` (managed via `pyproject.toml`). `uv sync` once,
  then `uv run python -m experiments.<module>`.
- **Symbolic:** `wolframscript -file analysis/<file>.wls`.
- **LaTeX:** `pdflatex` via `paper/build.sh`.
- **Site:** `./site/build.sh` rebuilds the paper, copies it into
  `site/report.pdf`, generates `site/build-info.js`. Deploy with
  `npx wrangler deploy`.

See [`CLAUDE.md`](CLAUDE.md) for the per-session research workflow.
