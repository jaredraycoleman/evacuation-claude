#!/usr/bin/env wolframscript
(* Verify Czyzowicz et al. (DISC 2014, Thm 7) lower bound for k=3 wireless.

   Lemma 6:   T >= 1 + x + 2 Sin[x k / 2]  for any x such that scan budget
              x k leaves nonzero unexplored perimeter (x < 2 pi / k).
   Theorem 7: maximize over x. For k=3, the optimum is interior, at the
              stationarity  cos(3 x / 2) = -1/3.

   So x_opt = (2/3) arccos(-1/3)                   approx 1.2739
     sin(3 x_opt / 2) = sqrt(1 - 1/9) = 2 sqrt 2 / 3 approx 0.9428
     LB = 1 + x_opt + 4 sqrt 2 / 3                 approx 4.1593.

   The chord term 2 sin(x k / 2) = 4 sqrt 2 / 3 bounds the FINDER's chord
   to the adversary's exit (finder at the first unexplored endpoint touched
   at time 1+x, exit at the other unexplored endpoint). Our UB at k=3 is
   roughly 4.2185 (with the same 1+x term of 2.4865 and a larger
   bottleneck chord of sqrt 3 ~ 1.7321). Gap = 0.0593.
*)

ClearAll[x, T];

k = 3;
Tlem6[x_] := 1 + x + 2 Sin[x k / 2];

(* Solve the stationarity symbolically. *)
stat = Reduce[D[Tlem6[x], x] == 0 && 0 < x < 2 Pi / k, x];
Print["Stationarity condition: ", stat];

xopt = x /. ToRules[stat];
Print["x_opt (symbolic): ", xopt];
Print["x_opt (numeric) : ", N[xopt, 30]];

LB = FullSimplify[Tlem6[xopt]];
Print["LB (symbolic)   : ", LB];
Print["LB (numeric)    : ", N[LB, 30]];

(* Confirmation of the claimed closed form. *)
claimed = 1 + (2/3) ArcCos[-1/3] + 4 Sqrt[2] / 3;
Print[];
Print["Claimed closed form: 1 + (2/3) arccos(-1/3) + 4 sqrt 2 / 3"];
Print["  = ", N[claimed, 30]];
Print["Match?            : ", N[claimed - LB, 20]];

(* Finder chord value at x_opt. *)
Print[];
Print["Chord-term         : 2 sin(3 x_opt / 2) = 4 sqrt 2 / 3"];
Print["                   = ", N[4 Sqrt[2] / 3, 20]];

(* For comparison, our UB's bottleneck chord is sqrt 3 ~ 1.7321,
   substantially larger. The LB's chord bound is loose by about 0.15. *)
Print[];
Print["Our UB's bottleneck chord is sqrt 3 =  ", N[Sqrt[3], 20]];
Print["LB's finder chord  =                 ", N[4 Sqrt[2]/3, 20]];
Print["Gap on chord side  =                 ", N[Sqrt[3] - 4 Sqrt[2]/3, 20]];
