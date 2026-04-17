# Literature notes

One short section per paper in `papers/`. Keep it skimmable — main model, key
results, assumptions, and anything flagged as open. The purpose is to make
cross-paper comparison easy when hunting for open problems.

---

## `papers/two-robots-F2F.pdf` — Brandt, Laufenberg, Lv, Stolz, Wattenhofer
*Collaboration without Communication: Evacuating Two Robots from a Disk.*

- **Model:** unit disk, exit on boundary (unknown position), 2 robots from
  center, unit speed, F2F communication only (co-located), visibility 0.
- **Contribution:** simpler algorithm than Czyzowicz et al. [8]; upper bound
  **5.625** (vs 5.628). Numerical optimality for the class with one symmetric
  straight-line detour per robot.
- **Tool of independent interest:** worst-case exit characterization —
  either non-differentiable motion or the angle condition
  $2\cos\beta + \cos\gamma = 1$.
- **Open / implicit:** tight lower bound; larger algorithm classes; non-disk
  shapes.

## `papers/ev-disc-k-robots.pdf` — Czyzowicz, Gąsieniec, Gorry, Kranakis, Martin, Pająk
*Evacuating Robots from an Unknown Exit in a Disk.*

- **Model:** unit disk, $k$ robots, wireless vs. non-wireless (F2F) models.
- **Results:**
  - Non-wireless: $3 + \tfrac{2\pi}{k}$ sufficient; $3 + \tfrac{2\pi}{k} - O(k^{-2})$ sometimes required.
  - Wireless: $3 + \tfrac{\pi}{k} + O(k^{-4/3})$ sufficient; $3 + \tfrac{\pi}{k}$ sometimes required.
- **Open:** exact constants for $k=2,3$; closing the $O(k^{-2})$ / $O(k^{-4/3})$ gaps.

## `papers/asymmetric-comm-SR.pdf` — Georgiou, Giachoudis, Kranakis
*Evacuation from a Disk for Robots with Asymmetric Communication (ISAAC 2022).*

- **Model:** 2 robots, one Sender-only, one Receiver-only (wireless), plus F2F when co-located. Unit disk, center start, unit speed.
- **Contribution:** shows SR model strictly between F2F and WiFi; upper bound
  less than $\pi + 2$.
- **Open:** tight bounds in SR; multi-robot generalizations; other fault
  patterns.

## `papers/priority-evac-queen.pdf` — Czyzowicz et al.
*Priority Evacuation from a Disk Using Mobile Robots.*

- **Model:** $n+1$ robots (1 queen + $n$ servants), wireless, minimize queen's
  arrival time at exit.
- **Results:** upper bound better than $2 + 4(\sqrt 2 - 1)\tfrac{\pi}{n}$ for
  $n \ge 4$; lower bound $2 + \tfrac{\pi}{n} + \tfrac{2}{n^2}$.
- **Open:** closing the gap; small $n$ (especially $n=1,2,3$); F2F variant.
