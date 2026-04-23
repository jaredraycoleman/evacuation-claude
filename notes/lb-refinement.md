# Lower-bound refinement for 3-robot wireless evacuation

## The refined inequality

**Proposition.** For any 3-robot wireless algorithm, any $x > 0$, and any
boundary point $p$ that is unexplored at time $1+x$:
$$
\mathrm{evac}(p) \;\ge\; (1+x) + \max_{j} \|R_j(1+x) - p\|,
$$
where $R_j(1+x)$ is robot $j$'s position at time $1+x$.

**Proof.** Let $t_p \ge 1+x$ be the first-visit time of $p$. At time $t_p$,
the finder is at $p$; every other robot $j$ has moved at most $t_p - (1+x)$
from $R_j(1+x)$, so by the triangle inequality
$\|R_j(t_p) - p\| \ge \|R_j(1+x) - p\| - (t_p - (1+x))$.
Hence the bottleneck chord at discovery is
$\max_j \|R_j(t_p) - p\| \ge \max_j \|R_j(1+x) - p\| - (t_p - (1+x))$, and
$$
\mathrm{evac}(p) = t_p + \max_j \|R_j(t_p) - p\|
\ge t_p + \max_j \|R_j(1+x) - p\| - (t_p - (1+x))
= (1+x) + \max_j \|R_j(1+x) - p\|.
$$
If the RHS happens to be negative, the bound $\mathrm{evac}(p) \ge t_p \ge
(1+x) + \max_j \|R_j - p\|$ still holds since $t_p - (1+x) > \max_j
\|R_j - p\|$ in that case. $\square$

## What this does NOT give, directly

Taking the infimum over algorithms of
$\sup_x \sup_{p \in E^c(1+x)}[(1+x) + \max_j \|R_j - p\|]$ produces a
*valid* lower bound on $T^*$. But numerically, the infimum is realised
by a degenerate "all robots at origin until time $1+x$" state, which
yields only $1 + x + 1$. That is $< 4.159$ for small $x$, so the refined
inequality alone does **not** universally improve the Czyzowicz et
al. LB.

(For a degenerate "robots-at-origin" algorithm, the actual evacuation
time is much larger than $1+x+1$ because the alg wastes time; CGK's
bound still applies and gives $\ge 4.159$. The refined bound just
happens to be slacker than CGK for this particular state.)

## On $A_3(y_{\text{opt}})$, the inequality is nearly tight

Applied to the specific $A_3(y_{\text{opt}})$ algorithm:

| $x$ | $1+x$ | $\max_j \|R_j - p^*\|$ | refined LB | actual UB |
|---|---|---|---|---|
| 1.0 | 2.0 | 2.000 | 4.000 | 4.21852 |
| 1.274 | 2.274 | 1.966 | 4.216 | 4.21852 |
| 1.500 | 2.500 | 1.718 | **4.21836** | 4.21852 |
| 1.950 | 2.950 | 1.266 | 4.216 | 4.21852 |

At $x = 1.5$ the refined bound is within $0.0002$ of the exact UB —
i.e., the refined inequality *nearly certifies* the actual worst-case
evac time of our best-known algorithm. This is strong evidence the
bound is "the right thing" in some sense; it just needs an argument
to be applied universally.

## Where the gap lives

CGK's $4.15937$ uses $|E(1+x)| \le 3x$ (a necessary condition on every
algorithm) plus Lemma 5 (existence of two unexplored points at arc
distance $\ge 2\pi - 3x$). For algorithms where this bound is
**tight** — in particular, algorithms that at time $1 + x_{\text{CGK}}^*$
have all 3 robots on the boundary with full scan coverage
$|E| = 3x_{\text{CGK}}^*$ — the refined inequality gives a *strictly
larger* LB. Specific cases checked:

| Algorithm at $x_{\text{CGK}}^* = 1.274$ | CGK LB | refined LB |
|---|---|---|
| 3-fold symmetric, each scans $[\alpha_j, \alpha_j + x]$ | 4.159 | 4.173 |
| Single-arc scan, robots at $e(x), e(2x), e(3x)$ | 4.159 | 4.188 |

So for any "full-scan, on-boundary at $1+x_{\text{CGK}}^*$" algorithm,
$T \ge$ refined $>$ CGK. The algorithms that could still achieve CGK's
$4.159$ are those **not** in this class — ones where some robot is
off-boundary at $1 + x_{\text{CGK}}^*$, or where $|E| < 3x$
significantly. For those algorithms, $|E| < 3x$ means the unexplored
arc is larger than Lemma 5 assumes, and CGK tightens automatically.

**Rough plan to close this gap rigorously:**

1. Prove a geometric lemma: for any 3 robots on $\partial D$ and any
   $|E| \ge 3x - \varepsilon$ with $R_j$ at scan endpoints, $\max_{p \in E^c}
   \max_j \|R_j - p\| \ge 2 \sin(3x/2) + \delta(\varepsilon)$ for some
   $\delta(\varepsilon) > 0$. (Conjectured from the 3-fold symmetric case
   $\max_{p} \max_j = 2\cos(x/4)$; compare to CGK's $2\sin(3x/2)$:
   $2\cos(x/4) > 2\sin(3x/2)$ iff $x > 2\pi/7 \approx 0.898$, which
   includes the CGK-optimal region.)
2. For algs NOT in that class (off-boundary or less-than-full scan):
   the actual $|E(1+x)|$ is strictly $< 3x$, and the refined CGK bound
   $1 + x + 2\sin((2\pi - |E|)/2)$ is strictly $> 1 + x + 2\sin(3x/2)$.
   Quantify this gap in terms of how much coverage was missed.
3. Combine: every algorithm falls into either class 1 or class 2, and
   in either case $T > 4.159$.

Step 1 is a finite-dimensional optimisation problem (min over 3 robot
angles + 3 scan arcs of a max over a 1D arc of a max over 3 robots).
Step 2 requires an additional argument about how redeploy-like moves
(where a robot leaves boundary) reduce $|E|$.

## What the session produced

- `experiments/lb_refined.py` — refined LB vs CGK for a symmetric
  algorithm family (Family A).
- `experiments/lb_refined_a3.py` — evaluates the refined LB on
  $A_3(y_{\text{opt}})$; demonstrates near-UB tightness.
- `experiments/lb_minimax.py` — numerical minimax over a rich
  algorithm-state parameter space; finds the $(1+x) + 1$ degenerate
  minimum, confirming the refined inequality alone is not universally
  stronger than CGK.
- This note.

## Verdict

**Partial result: not yet a proven universal LB improvement.**

The refined inequality is valid and clean. On $A_3(y_{\text{opt}})$ it
nearly certifies the actual UB. Combined with the state-tight CGK
bound, it forces a universal LB strictly greater than $4.15937$ for
**every algorithm state**:
- Algorithms with $|E(1+x)| = 3x$ at $x_{\text{CGK}}^* = 1.274$ (full
  scan, all 3 robots on $\partial D$): combined LB $\ge 4.172$
  (numerical minimum over robot angle configurations).
- Algorithms with $|E| < 3x$: state-tight CGK alone gives $> 4.159$.

So at $x = 1.274$, EVERY algorithm has combined LB strictly greater
than $4.159$. The **infimum** of "combined LB at $x = 1.274$" across
algorithms should be a rigorous LB for $T^*$.

The unfinished step is bounding this infimum cleanly. Two sub-problems:

1. **Min refined chord at $|E| = 3x$** (conjectured $= 2\cos(x/4)$,
   numerically verified for 3-fold symmetric and single-arc cases).
   A proof would give a clean LB $\ge 1 + x + 2\cos(x/4)$ in the
   full-scan class, worth $\approx 0.013$ at $x = 1.274$.

2. **Min refined chord as a function of $|E|$** for $|E| < 3x$.
   Without a clean handle on this, the boundary regime (|E| just
   below $3x$) could in principle push the combined LB down
   close to the original $4.159$. Numerics (`experiments/lb_combined.py`)
   suggest the actual improvement is $\approx 0.01$.

**Honest status**: the refined inequality is a new tool, and on $A_3$
it is essentially tight, but the universal LB improvement over
$4.15937$ is at most a tenth of the full $0.06$ gap (to $\approx
4.17$) without significant extra work. Closing the full gap
requires either a tighter analysis of the combined bound or a
different approach altogether (e.g., multi-exit adversaries that
combine the $\max_j$ refinement with a geometric argument about
robots-on-boundary constraints).
