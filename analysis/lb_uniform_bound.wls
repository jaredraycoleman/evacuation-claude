#!/usr/bin/env wolframscript
(* Compute the minimum of combined = max(state-tight CGK, refined)
   over the "balanced 2-coincident with slack" family.

   Parameter: tau in [0, x^*], representing per-robot slack.
   Configuration:
     R_1 = (1 - tau) e(d),  R_2 = R_3 = (1 - tau) e(0)
     scan arcs of length L = x^* - tau each, disjoint,
     d = pi + L/2 (optimal for this subfamily).
     |E| = 3L = 3(x^* - tau).

   Refined chord (sup over unexplored, by 2-coincident analysis with
   inward perturbation):
     h(tau) = sqrt{(1-tau)^2 + 2(1-tau) cos(L/2) + 1}
            = sqrt{(1-tau)^2 + 2(1-tau) cos((x^* - tau)/2) + 1}.
   This is the chord from R_1 = (1-tau) e(d) to p = e(L), which by
   symmetry equals the chord from R_2 = (1-tau) e(0) to p = e(d-L);
   both are the binding chords for the 2-coincident symmetric case.

   CGK(tau) = 1 + x^* + 2 sin((2 pi - 3L)/2) = 1 + x^* + 2 sin((2 pi - 3x^* + 3 tau)/2).

   combined(tau) = max(1 + x^* + h(tau), CGK(tau)).

   Minimise over tau in [0, x^*].
*)

ClearAll[x, tau, L, d, h, cgk, combined];

xStar = (2/3) ArcCos[-1/3];
Print["x^* = ", N[xStar, 15]];

(* Refined chord as function of tau *)
chord2coincBalanced[tau_] := Sqrt[
  (1 - tau)^2 + 2 (1 - tau) Cos[(xStar - tau)/2] + 1
];

refinedBalanced[tau_] := 1 + xStar + chord2coincBalanced[tau];

cgk[tau_] := 1 + xStar + 2 Sin[(2 Pi - 3 (xStar - tau))/2];

combinedBalanced[tau_] := Max[refinedBalanced[tau], cgk[tau]];

Print[];
Print["At tau = 0 (scan-tight):"];
Print["  CGK     = ", N[cgk[0], 10]];
Print["  Refined = ", N[refinedBalanced[0], 10]];
Print["  Combined= ", N[combinedBalanced[0], 10]];
Print[];

Print["At tau where CGK = Refined (crossover):"];
tauXover = tau /. FindRoot[refinedBalanced[tau] == cgk[tau], {tau, 0.01}];
Print["  tau* = ", N[tauXover, 12]];
Print["  CGK(tau*)     = ", N[cgk[tauXover], 12]];
Print["  Refined(tau*) = ", N[refinedBalanced[tauXover], 12]];
Print["  Combined(tau*)= ", N[combinedBalanced[tauXover], 12]];
Print[];

Print["Sweep of combined(tau) for tau in [0, 0.05]:"];
Print["  tau    | CGK       | Refined   | Combined"];
Do[
  Print["  ", N[tVal, 5], " | ",
        N[cgk[tVal], 10], " | ",
        N[refinedBalanced[tVal], 10], " | ",
        N[combinedBalanced[tVal], 10]],
  {tVal, Range[0.0, 0.05, 0.005]}
];

(* Exact minimum via symbolic/numeric optimization. *)
Print[];
Print["Minimising combined(tau) on [0, x^*]:"];
tauMin = FindMinimum[combinedBalanced[tau], {tau, 0.01}];
Print["  min = ", tauMin];

(* As a sanity check, the minimum should be at the crossover point. *)
Print[];
Print["Min combined over balanced 2-coincident slack family:"];
Print["  value = ", N[combinedBalanced[tauXover], 15]];
Print["  improvement over CGK: ", N[combinedBalanced[tauXover] - cgk[0], 10]];

(* Compare to the scan-tight (Case 2) bound: *)
Print[];
Print["Scan-tight (Case 2) bound: ", N[refinedBalanced[0], 10], " (= 1 + x^* + 2 cos(x^*/4))"];
Print["Balanced-slack crossover:  ", N[combinedBalanced[tauXover], 10]];
Print["CGK (old LB):              ", N[cgk[0], 10]];
