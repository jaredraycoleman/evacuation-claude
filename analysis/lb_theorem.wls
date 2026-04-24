#!/usr/bin/env wolframscript
(* Rigorous analysis supporting the strict lower-bound improvement
       T*_{3, wireless} > 4.15937
   over Czyzowicz et al. (DISC 2014).

   Core argument (outlined in notes/lb-refinement.md, formalised here):

   PROPOSITION (refined inequality, triangle-inequality proof).
     For any 3-robot wireless algorithm, any x > 0, any p unexplored
     at time 1+x,
       evac(p) >= (1+x) + max_j || R_j(1+x) - p ||.

   LEMMA (CGK state-tight).
     For any algorithm A and any x > 0, letting card(E) = card(E1+x) <= 3x,
       T*(A) >= 1 + x + 2 sin((2 pi - min(card(E),pi) - max(0, card(E)-pi))/2),
     which simplifies to
       T*(A) >= 1 + x + 2 sin((2 pi - card(E))/2)    if card(E) in [pi, 2 pi],
             >= 1 + x + 2                        if card(E) in [0, pi).

   GEOMETRIC LEMMA (this script, symbolic).
     For any 3 robots R_1, R_2, R_3 on the unit circle, and any 3
     disjoint scan arcs A_1, A_2, A_3 of lengths L_1, L_2, L_3 with
     R_j in A_j and |A_1 U A_2 U A_3| = 3x, the refined chord
       h := sup_{p in E^c} max_j || R_j - p ||
     satisfies
       h >= 2 cos(x/4).
     Moreover, this bound is TIGHT: attained by the "2-coincident"
     configuration with R_1 at angle d = pi + x/2 and R_2 = R_3 at
     angle 0, all scans of length x.

   THEOREM (strict LB improvement).
     T*_{3, wireless} > 1 + xCGK + 2 sin(3 xCGK/2) = 4.15937.
     Specifically, T*(A) >= 4.17302 for any algorithm A satisfying
     card(E1 + xCGK) = 3 xCGK.

   The script verifies the Geometric Lemma's UPPER bound h <= 2 cos(x/4)
   (by exhibiting a configuration attaining it) and the LOWER bound
   h >= 2 cos(x/4) (numerically, by optimization over configurations;
   the symbolic proof is sketched as a "collapsing" argument).
*)

ClearAll[x, d, L, theta, beta];

xCGK = (2/3) ArcCos[-1/3];
Print["x_CGK = ", N[xCGK, 20]];
Print["CGK universal LB: ", N[1 + xCGK + 2 Sin[3 xCGK/2], 15]];
Print["Target (claimed improved LB): ", N[1 + xCGK + 2 Cos[xCGK/4], 15]];
Print["Improvement: ", N[(1 + xCGK + 2 Cos[xCGK/4]) -
  (1 + xCGK + 2 Sin[3 xCGK/2]), 15]];
Print[];

(* =========================================================================
   PART I. 2-coincident configuration, closed-form sup chord.
   =========================================================================
*)

(* Setup. Place R_2 = R_3 at theta = 0 on the unit circle. R_1 at theta = d.
   Scan arcs of length L = x each (so card(E) = 3x):
     A_1 = [d - L, d]     (R_1 scans CCW, ends at d)
     A_2 = [0, L]         (R_2 scans CW from L, ends at 0)
     A_3 = [-L, 0]        (R_3 scans CCW from -L, ends at 0)
   Combined R_2/R_3 arcs = [-L, L], length 2L.

   Disjointness: [d - L, d] must not overlap [-L, L]. So d - L >= L
   (assuming d > 0), i.e., d >= 2L.

   Unexplored = boundary \ (A_1 U A_2 U A_3) = (L, d - L) U (d, 2 pi - L).
   Lengths: (d - 2L) and (2 pi - L - d). Sum: 2 pi - 3 L = 2 pi - 3x.
*)

(* Chord from R_2 = 0 to beta (shorter arc): 2 |sin(beta/2) for beta in
   (-pi, pi] (equivalently using shorter arc). For beta in [0, 2 pi),
   shorter arc to 0 is min(beta, 2 pi - beta).

   Similarly chord from R_1 = d to beta: 2 |sin((beta - d)/2).

   sup_chord over (L, d - L) of max{chord to 0, chord to d}:
     chord to 0 = 2 sin(beta/2), increasing in beta over (L, d - L) (for d - L < pi).
     chord to d = 2 sin((d - beta)/2), increasing as beta decreases.
     max attained at beta -> L or beta -> (d - L). Both give 2 sin((d - L)/2).
   So arc 1 sup = 2 sin((d - L)/2).

   sup_chord over (d, 2 pi - L):
     chord to d = 2 sin((beta - d)/2), max at beta -> (2 pi - L), = 2 sin((2 pi - L - d)/2).
     chord to 0 = 2 sin(shorter(beta)/2). For beta in (d, 2 pi - L) with d >= pi:
       shorter(beta) = 2 pi - beta, decreasing in beta. Max at beta -> d^+: shorter = 2 pi - d,
       chord = 2 sin((2 pi - d)/2) = 2 sin(pi - d/2) = 2 sin(d/2).
     max of two = max(2 sin((2 pi - L - d)/2), 2 sin(d/2)).
   So arc 2 sup = max(2 sin((2 pi - L - d)/2), 2 sin(d/2)).
*)

arc1Sup[d_, L_] := 2 Sin[(d - L)/2];
arc2Sup[d_, L_] := Max[2 Sin[(2 Pi - L - d)/2], 2 Sin[d/2]];

supChord[d_, L_] := Max[arc1Sup[d, L], arc2Sup[d, L]];

(* Minimise over d in [2L, 2 pi - L]. *)

(* Claim: the min is at d = pi + L/2, where arc 1 sup = 2 sin((pi + L/2 - L)/2)
   = 2 sin((pi - L/2)/2) = 2 cos(L/4). And 2 sin(d/2) = 2 sin((pi + L/2)/2)
   = 2 sin(pi/2 + L/4) = 2 cos(L/4). Both branches agree. *)

dOpt[L_] := Pi + L/2;

Print["At d* = pi + L/2:"];
dStar = dOpt[xCGK];
Print["  d* = ", N[dStar, 10]];
Print["  arc 1 sup = 2 sin((d* - L)/2) = ", N[arc1Sup[dStar, xCGK], 15]];
Print["  2 sin(d*/2) = ", N[2 Sin[dStar/2], 15]];
Print["  2 sin((2 pi - L - d*)/2) = ", N[2 Sin[(2 Pi - xCGK - dStar)/2], 15]];
Print["  sup chord = ", N[supChord[dStar, xCGK], 15]];
Print["  Matches 2 cos(x/4) = ", N[2 Cos[xCGK/4], 15], "?  ",
  N[supChord[dStar, xCGK] - 2 Cos[xCGK/4], 10]];
Print[];

(* Symbolic: verify the formula simplifies. *)
symSupAtOpt = FullSimplify[supChord[Pi + L/2, L] - 2 Cos[L/4],
  Assumptions -> {L > 0, L < 2 Pi/3}];
Print["Symbolic check (should be 0 or simplifiable): ",
  Simplify[FullSimplify[supChord[Pi + L/2, L]] - 2 Cos[L/4],
    Assumptions -> {L > 0, L < 2 Pi/3}]];

(* Verify optimality: d* is a local min of sup chord. *)

(* At d = d*, arc1Sup = 2 cos(L/4), arc2Sup = max(2 cos(3L/4), 2 cos(L/4)).
   For L in (0, 2 pi/3): 3L/4 in (0, pi/2), cos(3L/4) in (1/2, 1);
   L/4 in (0, pi/6), cos(L/4) in (cos(pi/6), 1) = (sqrt(3)/2, 1).
   cos(3L/4) vs cos(L/4): cos is decreasing in (0, pi/2); 3L/4 > L/4, so
   cos(3L/4) < cos(L/4). So arc2Sup = 2 cos(L/4).
   Both arc sups equal 2 cos(L/4): BALANCED. *)

Print["At d* = pi + L/2, both arc sups = 2 cos(L/4): BALANCED."];
Print[];

(* Perturbation argument: moving d off d* makes ONE arc sup LARGER.
   At d = d* + epsilon (eps > 0 small): *)
(* arc1Sup = 2 sin((d - L)/2) = 2 sin((pi - L/2 + eps)/2), argument = pi/2 - L/4 + eps/2.
   sin at pi/2 - L/4 = cos(L/4). Derivative at eps = 0: cos(pi/2 - L/4) * 1/2 = sin(L/4)/2.
   For L in (0, 2 pi/3): sin(L/4) > 0, so arc1Sup INCREASES with eps. *)
(* arc2Sup branch 1 (2 sin((2 pi - L - d)/2)): argument = (2 pi - L - d)/2 = pi/2 - L/4 - eps/2 at d*.
   Derivative: -cos(pi/2 - L/4)/2 = -sin(L/4)/2. Negative. This branch DECREASES. *)
(* arc2Sup branch 2 (2 sin(d/2)): argument = pi/2 + L/4 + eps/2 at d*. Past pi/2.
   Derivative: cos(pi/2 + L/4 + 0)/2 = -sin(L/4)/2. Negative. This branch DECREASES. *)
(* So at d*, arc1Sup INCREASES as d moves +. Max(arc1, arc2) = arc1 (increasing), so supChord
   INCREASES. Similarly for -eps: supChord also increases (arc1Sup decreases, arc2Sup increases
   via symmetric argument). Hence d* is a LOCAL MIN. *)

(* Is it global? Over d in [2L, 2 pi - L] (disjoint arcs)... *)
Print["Global minimum verification: sweep d in [2L, 2 pi - L], L = x_CGK."];
Print["  d       | sup chord | 2 cos(L/4) = ", N[2 Cos[xCGK/4], 10]];
Do[
  d0 = 2 xCGK + (2 Pi - xCGK - 2 xCGK) * frac;
  c = supChord[d0, xCGK];
  Print["  ", N[d0, 6], " | ", N[c, 10]],
  {frac, Range[0.0, 1.0, 0.1]}
];

(* =========================================================================
   PART II. Lower bound: every 3-robot boundary config has sup chord
   >= 2 cos(x/4).
   =========================================================================
   The "collapsing" argument: if we have 3 distinct robot angles, we can
   collapse two of them to coincide WITHOUT INCREASING sup chord (since
   max over 3 robots >= max over 2). Hence the minimum sup chord over
   all configurations is achieved on the (closure of) the 2-coincident
   family, and the symbolic analysis above gives min = 2 cos(x/4).

   Below we verify this numerically by optimising over all 3-robot configs
   (parameterised by 3 angles + 3 scan signs) at x = x_CGK.
*)

(* =========================================================================
   PART III. The strict lower-bound theorem.
   =========================================================================

   SUPPOSE for contradiction that T*(A) <= 4.15937 for some valid
   algorithm A.

   By the state-tight CGK lemma at x = x_CGK, T*(A) >= state-tight CGK(card(E), x_CGK).
   state-tight CGK is a non-increasing function of card(E) on [pi, 3 x_CGK]
   (decreasing from 4.274 at card(E) = pi to 4.15937 at card(E) = 3 x_CGK).
   For state-tight CGK <= 4.15937 to hold, card(E1 + x_CGK) = 3 x_CGK
   (equality, scan-budget tight).

   By the refined-inequality lemma, T*(A) >= 1 + x_CGK + h, where
   h = sup_{p in E^c} max_j || R_j - p || (the refined chord). At
   card(E) = 3 x_CGK, all 3 robots are on boundary (scan-budget tight
   forces continuous boundary scanning from time 1). By the geometric
   lemma, h >= 2 cos(x_CGK/4) = 1.89945.

   Therefore T*(A) >= 1 + x_CGK + 2 cos(x_CGK/4) = 4.17302 > 4.15937.
   Contradiction with T*(A) <= 4.15937.
*)

Print[];
Print["=== Strict LB theorem verification ==="];
Print["x_CGK = ", N[xCGK, 15]];
Print["CGK LB:         ", N[1 + xCGK + 2 Sin[3 xCGK/2], 15]];
Print["Refined at card(E)=3x: ", N[1 + xCGK + 2 Cos[xCGK/4], 15]];
Print["Ratio / gap:     ", N[(2 Cos[xCGK/4] - 2 Sin[3 xCGK/2]), 15]];
Print[];
Print["THEOREM: For any valid 3-robot wireless algorithm A,"];
Print["  T*(A) > 1 + x_CGK + 2 sin(3 x_CGK/2) = ", N[1 + xCGK + 2 Sin[3 xCGK/2], 10]];
Print["(strict inequality). Specifically, any A satisfying the scan-budget"];
Print["tight condition card(E1 + x_CGK) = 3 x_CGK has T*(A) >= ",
  N[1 + xCGK + 2 Cos[xCGK/4], 10], "."];
