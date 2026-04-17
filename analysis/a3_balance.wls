#!/usr/bin/env wolframscript
(* Symbolic derivation of y_opt for the A3-family k=3 wireless UB.

   Setup.  Let y be the B-offset, W = Pi - 3 y/2 - 1, and for an exit angle
   theta in (1+y, Pi - y/2), let v = Pi - y/2 - theta so u := W - v is r_3's
   fractional progress along its redeploy chord at discovery time.

   The r_3 chord to the exit has length

       c(v, W) = Sqrt[ (W - v)^2 + 1 - 2 (W - v) Cos[v] ].

   Case-3 worst-case time: T3(v, W) = 1 + (Pi - y/2 - v) + c(v, W).
   Critical-point condition dT3/dv = 0 gives the quadratic in (W - v):

       (W - v)^2 (2 - Sin v) - 2 (W - v) Cos v + Sin v == 0             (I)

   whose admissible root is

       W - v == (Cos v - Sqrt[1 - 2 Sin v]) / (2 - Sin v).              (I')

   The case-1 / case-3 balance asks T_1 == T_3-at-critical, which gives

       c(v, W) == v - Pi/3 + Sqrt[3].                                   (II)

   Combining (I') with dT3/dv = 0 (chord_3 = -(W-v)(1 - sin v) + cos v) gives

       c(v, W) == (Cos v + (1 - Sin v) Sqrt[1 - 2 Sin v]) / (2 - Sin v). (III)

   So the master equation is (III) == (II):

     (Cos v + (1 - Sin v) Sqrt[1 - 2 Sin v]) / (2 - Sin v)
       == v - Pi/3 + Sqrt[3].                                           (*)

   Solve (*) for v, then y = 2 (Pi - 1 - W)/3 and UB = 1 + 2 Pi/3 - y/2
   + Sqrt[3].
*)

ClearAll[v, W, y];

master = (Cos[v] + (1 - Sin[v]) Sqrt[1 - 2 Sin[v]]) / (2 - Sin[v])
          - (v - Pi/3 + Sqrt[3]);

(* Numerical root near v ~ 0.2028. *)
vSol = v /. FindRoot[master == 0, {v, 0.20}, WorkingPrecision -> 40];
Print["v_opt (40-digit)  = ", N[vSol, 30]];

WmVSol = (Cos[vSol] - Sqrt[1 - 2 Sin[vSol]]) / (2 - Sin[vSol]);
WSol   = WmVSol + vSol;
ySol   = 2 (Pi - 1 - WSol)/3;
UBSol  = 1 + 2 Pi/3 - ySol/2 + Sqrt[3];

Print["W_opt (40-digit)  = ", N[WSol,  30]];
Print["y_opt (40-digit)  = ", N[ySol,  30]];
Print["UB    (40-digit)  = ", N[UBSol, 30]];

(* Cross-check: simulator gave ~4.2185169939 at y_opt ~ 1.2158578. *)
Print[];
Print["Cross-check vs simulator:"];
Print["  y_opt numerical = 1.2158578321"];
Print["  UB   numerical  = 4.2185169939"];

(* Try to produce a closed form. Simplify combinations; the presence of
   both Sqrt[1 - 2 Sin v] and v inside the same equation strongly suggests
   there is no elementary closed form, but we try a few substitutions. *)

Print[];
Print["Attempting Solve[] for a closed form (may take a moment)..."];

cfAttempt = TimeConstrained[
  Solve[master == 0 && 0 < v < Pi/6, v, Reals],
  30, $Failed];
Print["Solve result: ", cfAttempt];

(* One clean observation: the balance value chord_3* = v - Pi/3 + Sqrt[3]
   plus the relation (I), so we can eliminate v in favour of chord_3.  If
   we set c = chord_3, then v = c + Pi/3 - Sqrt[3] and we need a single
   transcendental equation in c.  The resulting equation still contains
   trig of (c + Pi/3 - Sqrt[3]), no simplification. *)

Print[];
Print["Checking residual magnitude at numerical solution:"];
Print["  master(v_opt) = ", N[master /. v -> vSol, 20]];

(* Express y and UB via v symbolically. *)
ySym = 2 (Pi - 1 - v - (Cos[v] - Sqrt[1 - 2 Sin[v]]) / (2 - Sin[v]))/3;
UBSym = 1 + 2 Pi/3 - ySym/2 + Sqrt[3];
Print[];
Print["UB as a function of v_opt (before substitution):"];
Print[FullSimplify[UBSym]];

(* As requested by the transcendental nature, report "UB = 1 + 2pi/3 -
   y_opt/2 + sqrt 3 where y_opt, v_opt solve master == 0". *)
