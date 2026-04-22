#!/usr/bin/env wolframscript
(* Priority evacuation with n = 1 servant (Q + 1 servant, wireless).
   Verify:
   (1) UB = 1 + 2 pi / 3 + sqrt 3 via 2-robot wireless adapted.
   (2) LB = 1 + pi via Lemma-6-style adversary.
*)

ClearAll[x, y, theta];

UB = 1 + 2 Pi/3 + Sqrt[3];
Print["UB (2-robot wireless adapted) = 1 + 2 pi / 3 + sqrt 3 = ", N[UB, 15]];

(* The LB argument: at time 1 + x, total boundary visited by 2 unit-speed
   robots is at most 2 x (each spent at most x time on the boundary after
   the 1-unit radial approach). Unexplored set has length u = 2 pi - 2 x.
   Lemma 5 (Czyzowicz et al.) -- applicable while u <= pi, i.e. x >= pi / 2 --
   yields two unexplored points P1, P2 at shorter-arc distance >= u, hence
   Euclidean chord |P1 - P2| = 2 sin(u/2).

   For priority evacuation, adversary picks the exit at whichever of P1, P2
   maximizes T_Q. For any queen position p in the closed disk,
       max(|p - P1|, |p - P2|) >= (1/2) (|p - P1| + |p - P2|)
                                >= (1/2) |P1 - P2| = sin(u/2)
   by triangle inequality.  So chord_Q >= sin(u/2) and t_find >= 1 + x,
   giving
       T_Q >= 1 + x + sin(u/2) = 1 + x + sin(pi - x) = 1 + x + sin(x).
   Maximising over x in [pi/2, pi),
       sup_{x in [pi/2, pi)} (1 + x + sin x) = 1 + pi (limit as x -> pi).
*)

LBfunc[x_] := 1 + x + Sin[Pi - x];    (* = 1 + x + sin x *)
Print[];
Print["LB function LB(x) = 1 + x + sin(x):"];
Do[Print[
   "  x = ", N[xv, 5],
   "    u = 2 pi - 2 x = ", N[2 Pi - 2 xv, 5],
   "    LB(x) = ", N[LBfunc[xv], 10]],
  {xv, {Pi/2, 2 Pi/3, 3 Pi/4, 5 Pi/6, Pi - 1/10, Pi - 1/100}}];

Print[];
Print["sup_{x in [pi/2, pi)} LB(x) = 1 + pi = ", N[1 + Pi, 15]];

gap = UB - (1 + Pi);
Print[];
Print["Current bracket for T^*_{priority, n=1, wireless}:"];
Print["  ", N[1 + Pi, 10], " <= T^* <= ", N[UB, 10]];
Print["  gap = ", N[gap, 10]];
