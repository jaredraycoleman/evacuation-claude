#!/usr/bin/env wolframscript
(* Symbolic analysis of the combined lower bound
       combined(state, x) = max(state-tight-CGK(|E|, x), refined(state, x))
   at x = x_{CGK}^* = (2/3) arccos(-1/3).

   GOAL: lower bound the infimum of combined over valid algorithm states.

   Conjecture (from numerics at |E| = 3x):
       min_{states with |E|=3x} refined = 2 cos(x/4) + 1 + x,
   attained by a "2 coincident robots + 1 distant" configuration.

   At x = x_{CGK}^*, this gives combined LB = 1 + x_{CGK}^* + 2 cos(x_{CGK}^*/4)
                                           ~ 4.17302.
*)

ClearAll[x, e, t, beta];

(* --- CGK state-tight bound, as function of e = |E| --- *)
cgkBound[e_, x_] := If[e >= Pi,
  1 + x + 2 Sin[(2 Pi - e)/2],
  1 + x + 2
];

(* Numerical CGK optimum *)
xCGK = (2/3) ArcCos[-1/3];
xCGKn = N[xCGK, 30];
cgkAtxCGK = 1 + xCGK + 2 Sin[3 xCGK / 2];
cgkAtxCGKn = N[cgkAtxCGK, 30];
Print["x_CGK = ", xCGKn];
Print["CGK published bound: ", cgkAtxCGKn];
Print[];

(* --- Refined bound for 2-coincident configuration.

   Robots at angles {theta1, theta2, theta2}, with 3 disjoint scan arcs
   of length L each (total |E| = 3L), ending at these positions:
     arc_1 = [theta1 - L, theta1]   (robot 1 scans CCW ending at theta1)
     arc_2 = [theta2, theta2 + L]   (robot 2 scans CW ending at theta2)
     arc_3 = [theta2 - L, theta2]   (robot 3 scans CCW ending at theta2
             = SAME endpoint as robot 2; the arcs are disjoint since
             arc_3 is "just before" theta2 and arc_2 is "just after")

   WAIT: arc_2 and arc_3 both end at theta2. For disjoint arcs of
   length L, we need arc_3 = [theta2 - L, theta2], arc_2 = [theta2, theta2 + L],
   which ARE disjoint (share only endpoint). But then robots 2 and 3 both
   at theta2. Robot 1's arc = [theta1 - L, theta1], ending at theta1.

   For all disjoint (pairwise), need the arcs [theta1 - L, theta1],
   [theta2 - L, theta2], [theta2, theta2 + L] to not overlap.
   - [theta2 - L, theta2] and [theta2, theta2 + L] share theta2 only.
   - [theta1 - L, theta1] must not overlap either; thus theta1 is OUTSIDE
     [theta2 - L, theta2 + L].
*)

(* For the 2-coincident config, parameter is separation angle
   d = theta1 - theta2 (modulo 2 pi). Unexplored arcs:
     (theta1, theta2 - L)  = length (theta2 - L) - theta1 = (d) - L - ... wait

   Let me set theta2 = 0 WLOG. theta1 = d (positive, with d > L to avoid
   overlap with [0 - L, 0 + L] = [-L, L]).
   Arcs: [d - L, d], [-L, 0], [0, L].
   Total coverage length: 3L (if disjoint). Need d >= L (so [d - L, d]
   doesn't overlap [-L, L] which is two adjacent arcs merging to [-L, L]).

   Actually [-L, 0] U [0, L] = [-L, L], length 2L. Plus [d - L, d] length L.
   For disjoint: d - L >= L, i.e., d >= 2L (or symmetrically d <= -2L +2pi).

   Unexplored arcs: (L, d - L) and (d, 2pi - L) (going CCW back to -L = 2pi - L).
   Lengths: (d - L) - L = d - 2L; (2pi - L) - d = 2pi - L - d.
   Check sum: (d - 2L) + (2pi - L - d) = 2pi - 3L. Matches 2pi - |E|.
*)

(* Robot angles: R1 = d, R2 = R3 = 0.
   Adversary picks beta in unexplored.
   By symmetry consider beta in (L, d - L) (arc between R_2/R_3 cluster
   scan and R_1's scan).

   Chord from R_2 = 0 to beta in (L, d - L): 2 sin(beta/2).
   Chord from R_1 = d to beta: 2 sin((d - beta)/2) (if d - beta in (0, pi)).

   Max chord over beta in (L, d - L):
     2 sin(beta/2) at beta -> L+: 2 sin(L/2). Increasing in beta.
     2 sin((d - beta)/2) at beta -> (d-L)^-: 2 sin(L/2). Increasing as beta decreases.

   Max of max: both are 2 sin((something)/2) where something = max(beta, d - beta).
   Max of 2 sin(max(beta, d-beta)/2) over beta in (L, d-L).
   For d - 2L > 0: max at beta = L (giving d - L) or beta = d - L (giving d - L).
   So max chord over this arc = 2 sin((d-L)/2).

   Similarly for the other unexplored arc (d, 2pi - L):
   Chord from R_1 = d: 2 sin((beta - d)/2). Max at beta -> (2pi - L)^-: 2 sin((2pi - L - d)/2).
   Chord from R_2 = 0: 2 sin(shorter(beta - 0)/2). For beta in (d, 2pi - L):
     beta in (d, 2pi - L), beta - 0 = beta > d > 0. If beta <= pi: shorter = beta.
     If beta > pi: shorter = 2pi - beta.
   Max: depends.

   Let me parameterize and compute.
*)

supChordConfig[d_, L_] := Module[{unexp1End, unexp2End, arc1Sup, arc2Sup},
  unexp1End = d - L;
  unexp2End = 2 Pi - L;
  (* Arc 1: (L, d - L). Sup of max_j chord.
     Robots at 0, d. Max chord from 0: 2 sin(beta/2), max at beta = d - L -> 2 sin((d-L)/2).
     Max chord from d: 2 sin((d - beta)/2), max at beta = L -> 2 sin((d-L)/2).
     Arc 1 sup: 2 sin((d-L)/2). *)
  arc1Sup = 2 Sin[(d - L)/2];
  (* Arc 2: (d, 2 pi - L). Robots at 0, d.
     Max chord from 0: 2 sin(shorter(beta)/2). For beta in (d, 2pi - L), beta shorter if < pi
     else 2pi - beta. Key points: pi if pi in arc.
     Max chord from d: 2 sin(shorter(beta - d)/2). Max at beta = 2pi - L -> 2 sin((2pi - L - d)/2).

     Let me just compute at a few points numerically. *)
  (* For the specific (d, L) of interest, use closed form if simple. *)
  arc2Sup = Max[
    2 Sin[Min[Pi, (d + 2 Pi - L)/2 - d]/2], (* crude *)
    2 Sin[(2 Pi - L - d)/2]
  ];
  Max[arc1Sup, arc2Sup]
];

(* Actually let me just do this numerically in Wolfram more carefully. *)

supChordNumeric[d_, L_] := Module[{betas, chord0, chordD, maxC},
  betas = Join[
    Table[L + (d - L - L) * f, {f, 0.001, 0.999, 0.01}], (* arc 1 *)
    Table[d + (2 Pi - L - d) * f, {f, 0.001, 0.999, 0.01}] (* arc 2 *)
  ];
  Max @@ Table[
    Max[
      2 Abs[Sin[beta/2]],
      2 Abs[Sin[(beta - d)/2]]
    ], {beta, betas}
  ]
];

(* Test at d = sep of interest. From Python: d ~ 3.779 at x = 1.274. *)
Lfix = 1.273754617454;
Print["L = x_CGK = ", N[Lfix, 8]];
Print["Testing various d (separations) for 2-coincident config:"];
Print["  d       | sup chord | 2cos(x/4)"];
Do[
  c = supChordNumeric[dVal, Lfix];
  Print["  ", N[dVal, 5], " | ", N[c, 6], " | ", N[2 Cos[Lfix/4], 6]],
  {dVal, {Pi, 3.3, 3.5, 3.779, 3.9, 2 Pi - 2 Lfix}}
];

(* Now let's characterize: for what d is sup chord MINIMIZED, and what is that min? *)

(* Alg picks d in [2L, 2pi - L] (disjointness constraint).
   For each d, sup chord = max over 2 unexplored arcs.
   Arc 1 sup = 2 sin((d - L)/2).
   Arc 2 sup = depends on whether pi in arc.

   Claim: the minimizer is where both arc sups are equal (balance).
*)

(* Set arc1 sup = arc2 sup. *)
(* Arc 2: for beta in (d, 2pi - L), chord from 0 has max at beta = ... *)
(* Chord from 0 at beta (shorter arc): 2 sin(min(beta, 2pi - beta)/2).
   For beta in (d, 2pi - L), shorter(beta) = beta if beta <= pi else 2pi - beta.
   Switching at beta = pi, if pi is in (d, 2pi - L). Needs d < pi < 2pi - L,
   i.e., L < pi, always true for L < 2 pi/3.

   So shorter-beta max at beta = pi, chord = 2 sin(pi/2) = 2. Unless pi is not in unexplored.

   Is pi in unexplored arc (d, 2pi - L)? Need d < pi < 2pi - L. Since L < pi, 2pi - L > pi. So need d < pi.

   If d < pi: pi in arc 2, chord from 0 at beta = pi is 2. Arc 2 sup = 2.
   If d >= pi: pi not in arc 2, need to maximize otherwise.
*)

Print[];
Print["For 2-coincident config with robots at 0 and d, L = x_CGK:"];
Print["  if d < pi: arc 2 contains pi, sup chord from 0 at beta = pi is 2."];
Print["    Adversary gets chord 2, refined LB = 1 + x + 2 = ", N[1 + Lfix + 2, 6]];
Print["    This is BIG, so alg won't choose d < pi."];
Print[];
Print["  if d >= pi: pi not in arc 2. Arc 2 sup from 0 at beta = 2pi - L - eps:"];
Print["    shorter(2pi - L) = 2pi - (2pi - L) = L, chord = 2 sin(L/2) = ",
  N[2 Sin[Lfix/2], 6]];
Print["    or at beta -> d+, chord from 0 = 2 sin(d/2) if d < pi, else 2 sin((2pi - d)/2)."];

(* For d >= pi, arc 2 = (d, 2pi - L). For beta in this interval:
     chord from 0: shorter(beta) in (2pi - (2pi - L), 2pi - d) = (L, 2pi - d).
                    = shorter decreasing from ??? as beta increases. For beta > pi,
                    shorter = 2pi - beta, decreasing.
                  Max at beta = d (just after): shorter = 2 pi - d. Chord = 2 sin((2pi - d)/2) = 2 sin(pi - d/2) = 2 sin(d/2).
     chord from d: 2 sin((beta - d)/2), increasing. Max at beta = 2pi - L: 2 sin((2pi - L - d)/2).
   Arc 2 sup = max(2 sin(d/2), 2 sin((2pi - L - d)/2)).
*)

(* At d = 2 pi - 2 L (symmetric arc 2 about pi): 2 sin(d/2) = 2 sin(pi - L), 2 sin((2pi - L - d)/2) = 2 sin((2 pi - L - 2pi + 2L)/2) = 2 sin(L/2).
   So arc 2 sup = 2 sin(pi - L) = 2 sin(L), since pi - L < pi/2 for L > pi/2 (not quite).
   Actually sin(pi - L) = sin(L) always.
   And we need to compare to arc 1 sup = 2 sin((d - L)/2) = 2 sin((2 pi - 3 L)/2) = 2 sin(pi - 3L/2) = 2 sin(3L/2).

   So at d = 2pi - 2L: arc 1 sup = 2 sin(3L/2) (= CGK's chord!).
   Arc 2 sup = 2 sin(L) (depends on L).
   Sup chord = max of these.
*)

Print[];
Print["At d = 2 pi - 2 L (symmetric placement):"];
Print["  arc 1 sup = 2 sin(3L/2) = ", N[2 Sin[3 Lfix/2], 6]];
Print["  arc 2 sup = 2 sin(L)    = ", N[2 Sin[Lfix], 6]];
Print["  sup chord = max of these = ", N[Max[2 Sin[3 Lfix/2], 2 Sin[Lfix]], 6]];

(* At x = 1.274: 2 sin(3x/2) = 1.886, 2 sin(x) = 1.914. sup = 1.914.
   But numerical found min = 1.899 at d ~ 3.779 != 2pi - 2L = 6.28 - 2.547 = 3.736. *)

Print[];
Print["Numerical min d search:"];

bestChord = Infinity;
bestD = Null;
Do[
  dVal = N[dVal0, 20];
  c = supChordNumeric[dVal, Lfix];
  If[c < bestChord,
    bestChord = c; bestD = dVal
  ],
  {dVal0, Range[2 Lfix + 0.01, 2 Pi - Lfix - 0.01, 0.01]}
];
Print["Min sup chord (numeric over d): ", N[bestChord, 8]];
Print["At d = ", N[bestD, 8]];
Print["2 cos(x/4) = ", N[2 Cos[Lfix/4], 8]];
Print["2 sin(pi/2 - x/4) = ", N[2 Sin[Pi/2 - Lfix/4], 8]];
