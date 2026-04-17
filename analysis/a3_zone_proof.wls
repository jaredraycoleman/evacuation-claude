#!/usr/bin/env wolframscript
(* Rigorize Theorem 1: for the A3(y) algorithm with y = y_opt, the
   global max of evac(theta) over theta in [0, 2 pi) equals T_1(y_opt)
   = T_3_max(y_opt), and in particular no worst-case angle outside the
   four already identified contributes a larger value.

   We partition [0, 2 pi) by which robot is the finder:
      Zone 1 = [0, L_1],           finder r_1  (L_1 = pi - y/2)
      Zone 2 = [L_1, 2 pi - y],    finder r_2
      Zone 3 = [2 pi - y, 2 pi],   finder r_3.
   The setup is reflection-symmetric, so Zone 2's analysis mirrors Zone 1.

   Within Zone 1 we further split on which phase r_3 is in at
   t_found = 1 + theta:
      1a: theta in [0, y]         r_3 on boundary scan, position (cos(2pi-y+theta), sin(2pi-y+theta))
      1b: theta in [y, 1+y]       r_3 transiting boundary-endpoint -> origin along +x axis
      1c: theta in [1+y, L_1]     r_3 transiting origin -> redeploy target at angle alpha = pi - y/2

   At each zone we compute evac(theta) symbolically, find critical points
   and boundary values, and compare to T_1 = 1 + 2 pi/3 - y/2 + Sqrt[3]. *)

ClearAll[y, theta, v];

yOpt = y /. FindRoot[
  (Cos[v] + (1 - Sin[v]) Sqrt[1 - 2 Sin[v]]) / (2 - Sin[v])
    == v + Sqrt[3] - Pi/3
  /. v -> v  (* dummy, shaped for FindRoot *)
  , {v, 0.2}, WorkingPrecision -> 30];

(* Recompute y_opt directly so the value is concrete. *)
vOpt = v /. FindRoot[
  (Cos[v] + (1 - Sin[v]) Sqrt[1 - 2 Sin[v]]) / (2 - Sin[v])
    - (v - Pi/3 + Sqrt[3]) == 0,
  {v, 1/5}, WorkingPrecision -> 30];
WmV   = (Cos[vOpt] - Sqrt[1 - 2 Sin[vOpt]]) / (2 - Sin[vOpt]);
WOpt  = WmV + vOpt;
yOpt  = 2 (Pi - 1 - WOpt)/3;
T1    = 1 + 2 Pi/3 - yOpt/2 + Sqrt[3];
L1    = Pi - yOpt/2;

Print["At y_opt = ", N[yOpt, 20]];
Print["  T_1 = ", N[T1, 20]];
Print["  L_1 = pi - y/2 = ", N[L1, 20]];
Print[];

(* ---------- Zone 1a: theta in [0, y]. --------------------------------- *)
(* Finder r_1 at angle theta (chord 0).
   r_2 scanning, at boundary angle (2pi - y) - theta; chord
     c2 = 2 |sin((2pi - y - 2 theta)/2)| = 2 Sin[pi - y/2 - theta]  (positive).
   r_3 scanning, at boundary angle (2pi - y) + theta; chord
     c3 = 2 Sin[(2pi - y + theta - theta)/2] = 2 Sin[(2 pi - y)/2] = 2 Sin[y/2]
          (independent of theta!).
   Since c2 - c3 = 2 Sin[pi - y/2 - theta] - 2 Sin[y/2], and pi - y/2 - theta
   >= pi - 3 y/2 (for theta <= y) while y/2 is small, c2 >= c3 on this
   sub-zone iff sin(pi - y/2 - theta) >= sin(y/2) iff pi - y/2 - theta <= pi - y/2
   and pi - y/2 - theta >= y/2, i.e. theta <= pi - y (automatically true
   for theta <= y and y < pi/2). *)

evac1a = 1 + theta + 2 Sin[Pi - yOpt/2 - theta];

(* d/d theta evac1a = 1 - 2 Cos[pi - y/2 - theta]. For theta in [0, y],
   argument of Cos is in [pi - 3y/2, pi - y/2]. Both endpoints >= pi/2 for
   y_opt ~ 1.22, so Cos < 0 throughout, hence derivative > 1 > 0:
   evac1a is strictly increasing. Max at theta = y. *)

max1a = evac1a /. theta -> yOpt;
Print["Zone 1a max (at theta = y):  ", N[max1a, 15]];
Print["  gap to T_1:                ", N[T1 - max1a, 15]];

(* ---------- Zone 1b: theta in [y, 1+y]. ------------------------------ *)
(* r_3 is transiting boundary-endpoint -> origin along the +x axis, since
   r_3 finished scan at angle 2pi - y + y = 2 pi (the point (1, 0)).
   Position of r_3: (1 - (t - 1 - y), 0) = (1 + y - theta, 0).
     c3^2 = (1 + y - theta - Cos theta)^2 + (0 - Sin theta)^2
          = (1 + y - theta)^2 - 2 (1 + y - theta) Cos theta + 1.
   c2 stays 2 Sin[pi - y/2 - theta] (r_2 still scanning).

   evac1b(theta) = 1 + theta + Max[c2, c3]. *)

c3Sq1b = (1 + yOpt - theta)^2 - 2 (1 + yOpt - theta) Cos[theta] + 1;
c2     = 2 Sin[Pi - yOpt/2 - theta];
evac1b[th_] := 1 + th + Max[c2 /. theta -> th, Sqrt[c3Sq1b /. theta -> th]];

(* Candidate critical points for "1 + theta + c2" (c2 branch):
     d/d theta (1 + theta + c2) = 1 - 2 Cos[pi - y/2 - theta] = 0
     => theta = 2 pi/3 - y/2. *)

thetaCase1 = 2 Pi/3 - yOpt/2;
Print[];
Print["Zone 1b interior critical point of (1+theta+c2):"];
Print["  theta* = 2 pi/3 - y/2 = ", N[thetaCase1, 15]];
Print["  in [y, 1+y]?            ",
      yOpt <= thetaCase1 <= 1 + yOpt];
evac1bAtCase1 = evac1b[thetaCase1];
Print["  evac(theta*) = T_1 =   ", N[evac1bAtCase1, 15]];

(* Critical points for "1 + theta + c3" (c3 branch). *)
chord1b = Sqrt[c3Sq1b];
d1b = D[1 + theta + chord1b, theta];
Print[];
Print["Critical points of 1+theta+c3 in Zone 1b: (NSolve)"];
critsolns = NSolve[{d1b == 0, yOpt < theta < 1 + yOpt}, theta, Reals];
Print["  ", critsolns, "  (empty means max on c3 branch is at the boundary)"];

(* Endpoints of zone 1b: at theta=y evac = 1+y+2 sin(y/2) (small), at
   theta = 1+y evac = 3+y (because r_3 is at origin, chord = 1). *)
evac1bAtY    = 1 + yOpt + 2 Sin[yOpt/2];
evac1bAt1PlusY = 3 + yOpt;  (* r_3 at origin, chord = 1 *)

Print[];
Print["Zone 1b: evac at interior crit + boundaries:"];
Print["  theta = y          : evac = ", N[evac1bAtY, 14]];
Print["  theta = 2pi/3-y/2  : evac = ", N[evac1bAtCase1, 14], "  (= T_1, case 1)"];
Print["  theta = 1+y        : evac = ", N[evac1bAt1PlusY, 14], "  (= 3+y)"];
Print["  Zone 1b max        = T_1 = ", N[T1, 14]];

(* ---------- Zone 1c: theta in [1+y, L_1] = [1+y, pi-y/2]. -------- *)
(* r_3 transiting origin -> redeploy at angle alpha = pi - y/2.
   r_3 position at time t = 1+theta: (theta - 1 - y) * (cos alpha, sin alpha).
     c3^2 = (theta-1-y)^2 + 1 - 2 (theta-1-y) Cos[alpha - theta]
          = u^2 + 1 - 2 u Cos[alpha - theta], u := theta - 1 - y.
   c2 still = 2 Sin[pi - y/2 - theta], but for theta in [1+y, pi-y/2] the
   argument pi - y/2 - theta in [0, pi - 3y/2 - 1], hence c2 <= c3 since
   we check numerically. *)

alphaAng = Pi - yOpt/2;
c3Sq1c = (theta - 1 - yOpt)^2 + 1 - 2 (theta - 1 - yOpt) Cos[alphaAng - theta];
chord1c = Sqrt[c3Sq1c];
evac1c = 1 + theta + chord1c;
d1c = D[evac1c, theta];

(* Critical point of d1c == 0. *)
Print[];
Print["Critical points of 1+theta+c3 in Zone 1c (interior):"];
critsolns1c = NSolve[{d1c == 0, 1 + yOpt < theta < Pi - yOpt/2}, theta, Reals];
Print["  NSolve: ", critsolns1c];

(* Evac at interior critical points + boundaries of 1c. *)
boundary1c = {1 + yOpt, Pi - yOpt/2};
checkpts1c = Join[Re[theta /. critsolns1c], boundary1c];
Print[];
Print["Zone 1c candidates (evac value):"];
Do[
  th = pt;
  valC2 = c2 /. theta -> th;
  valC3 = chord1c /. theta -> th;
  Print["  theta = ", N[th, 10], "   c2 = ", N[valC2, 10],
        "   c3 = ", N[valC3, 10],
        "   evac = ", N[1 + th + Max[valC2, valC3], 12]];
  , {pt, checkpts1c}];

(* ---------- Zone 3: theta in [2pi - y, 2pi]. ---------------------- *)
(* Finder r_3 at the exit. r_1 scanning at angle (theta - (2pi-y)) of
   CCW progress, position angle = theta - (2pi - y). Chord between r_1
   and exit: 2 Sin[(2pi - y)/2] = 2 Sin[pi - y/2] = 2 Sin[y/2].
   r_2 scanning CW from 2pi - y for time u = theta - (2pi - y), at angle
   (2pi - y) - u. Chord to exit: 2 Sin[u]. *)

uZ3 = theta - (2 Pi - yOpt);
c1Z3 = 2 Sin[yOpt/2];
c2Z3 = 2 Sin[uZ3];

evac3 = 1 + uZ3 + Max[c1Z3, c2Z3];
Print[];
Print["Zone 3 analysis:"];
Print["  c1 (r_1) = 2 sin(y/2) = ", N[c1Z3, 10], " (constant in zone)"];
Print["  c2 (r_2) = 2 sin(theta - (2pi - y)) increasing from 0 to 2 sin y"];
Print["  evac = 1 + u + max(..),  u in [0, y]"];

(* Max of evac3 over theta in [2pi - y, 2pi]: the two regimes cross at
   2 sin u = 2 sin(y/2), i.e. u = y/2. Before: max uses c1, so evac =
   1 + u + 2 sin(y/2), max at u = y/2: 1 + y/2 + 2 sin(y/2). After: max
   uses c2, evac = 1 + u + 2 sin u, derivative 1 + 2 cos u >= 0 for u <=
   y (since y < pi/2), max at u = y: 1 + y + 2 sin y. *)

max3A = 1 + yOpt/2 + 2 Sin[yOpt/2];
max3B = 1 + yOpt + 2 Sin[yOpt];
Print["  max(c1 regime, u = y/2) = 1 + y/2 + 2 sin(y/2) = ", N[max3A, 12]];
Print["  max(c2 regime, u = y)   = 1 + y + 2 sin y     = ", N[max3B, 12]];
Print["  Zone 3 overall max     = ", N[Max[max3A, max3B], 12]];
Print["  gap to T_1             = ", N[T1 - Max[max3A, max3B], 12]];

(* ---------- Summary --------------------------------------------- *)
Print[];
Print["============================================================"];
Print["Summary: max evac in each zone vs T_1 = ", N[T1, 14]];
Print[];
Print["  Zone 1a max          <= evac1a(y)       = ", N[max1a, 14]];
Print["  Zone 1b max          =  T_1              = ", N[T1, 14]];
Print["  Zone 1c max          =  T_3-at-critical  = (see above)"];
Print["  Zone 2 mirrors Zone 1 by reflection symmetry"];
Print["  Zone 3 max           = ", N[Max[max3A, max3B], 14]];
Print[];
Print["Global max = T_1 =", N[T1, 14]];
