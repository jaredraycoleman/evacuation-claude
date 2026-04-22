#!/usr/bin/env wolframscript
(* Priority evacuation with n = 2 servants (Q + 2 servants, wireless).

   UB: run the optimal 3-robot wireless algorithm (our A_3(y_opt)) and
   designate the queen as one of the three robots. Queen's worst-case
   arrival is bounded by the 3-robot wireless worst-case arrival, i.e.,
   1 + 2 pi / 3 - y_opt / 2 + sqrt 3 ~ 4.21852.

   LB: adapt Lemma 6 with half-chord. At time 1 + x, three unit-speed
   robots have been on the boundary for at most x time each (after the
   1-unit radial approach), so total boundary coverage <= 3 x. Unexplored
   set has length >= 2 pi - 3 x. For x in [pi/3, 2 pi/3], Lemma 5 yields
   two unexplored points at shorter-arc distance >= u := 2 pi - 3 x,
   whose Euclidean chord is 2 sin(u/2). Adversary commits to one of the
   two points; by triangle inequality the queen's chord to the exit is
   at least sin(u/2). Hence
       T_Q >= 1 + x + sin((2 pi - 3 x)/2) = 1 + x + sin(pi - 3 x/2)
           = 1 + x + sin(3 x / 2).
   Maximise over x in [pi/3, 2 pi/3].
*)

LBfunc[x_] := 1 + x + Sin[3 x / 2];
dLB = D[LBfunc[x], x];
Print["Derivative of LB function: ", FullSimplify[dLB]];

stat = Reduce[dLB == 0 && Pi/3 < x < 2 Pi / 3, x];
Print["Stationarity: ", stat];

xOpt = x /. ToRules[stat];
Print["x_opt = ", xOpt, " ~ ", N[xOpt, 20]];
Print["LB(x_opt) = ", FullSimplify[LBfunc[xOpt]]];
Print["LB (numeric) = ", N[LBfunc[xOpt], 20]];

UB = 1 + 2 Pi / 3 - 1.2158578321292429 / 2 + Sqrt[3];
Print[];
Print["UB (via our k=3 wireless bound) = ", N[UB, 15]];

gap = UB - LBfunc[xOpt];
Print["Gap = ", N[gap, 15]];
