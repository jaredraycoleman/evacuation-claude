#!/usr/bin/env wolframscript
(* Symbolic analysis of the refined LB for 3-robot wireless evacuation
   on the unit disk.

   We study the quantity
       h(state) := sup_{p in E^c} max_j || R_j - p ||
   for a "3-robot state" at time 1+x. Constrain all 3 robots to be on
   boundary with 3 disjoint scan arcs of lengths L_1, L_2, L_3 summing
   to |E| = L_1 + L_2 + L_3.

   GOAL: characterize h as function of scan-arc configuration; find
   its minimum over configurations with fixed |E| = 3x (full-scan
   case); compare with CGK's 2 sin(3x/2).
*)

ClearAll[x, L, alpha, beta, theta];

(* --- Basic chord formulas --- *)

(* Chord between two boundary points e(a), e(b) via shorter arc. *)
chordAB[a_, b_] := 2 Abs[Sin[(b - a) / 2]];

(* When a, b are specific, chord = 2 |sin((b-a)/2)| always holds if we
   take a, b in R. Note: this equals 2 sin(shorter/2) because
   |sin(theta/2)| for theta in [0, 2 pi] equals sin(min(theta,2pi-theta)/2). *)

(* --- 3-fold symmetric state, all scans CCW of length L.
       alpha_j = 2 pi (j-1)/3, scan arc [alpha_j, alpha_j + L], robot
       at R_j = e(alpha_j + L). Unexplored arcs have length 2pi/3 - L
       each. --- *)

phi[j_] := 2 Pi (j - 1) / 3;  (* scan start angle *)
R[j_, L_] := phi[j] + L;      (* robot angle at time 1+x *)
unexStart[j_, L_] := phi[j] + L;      (* unexplored arc j starts *)
unexEnd[j_] := phi[j + 1];             (* = phi[j] + 2 pi/3, and ends *)

(* For beta in unexplored arc j = (unexStart, unexEnd), compute max_k chord. *)
(* By symmetry we fix j = 1: beta in (L, 2pi/3). *)

cMax[beta_, L_] := Max @@ Table[chordAB[R[k, L], beta], {k, 1, 3}];

(* Supremum over beta in [L, 2 pi/3] of cMax[beta, L]. *)
(* By the analysis in notes/lb-refinement.md: the sup is attained at
   beta -> (2 pi/3)^-, where the FAR robot R_3 = e(4pi/3 + L) gives the
   largest chord. Let's confirm. *)

(* Chord from R_3 = e(4 pi/3 + L) to beta: *)
cR3[beta_, L_] := chordAB[4 Pi/3 + L, beta];
(* = 2 |sin((4 pi/3 + L - beta)/2)|.
   For beta in [L, 2 pi/3], L in (pi/3, 2pi/3):
   4 pi/3 + L - beta in [4 pi/3 + L - 2 pi/3, 4 pi/3] = [2 pi/3 + L, 4 pi/3].
   Divided by 2: [pi/3 + L/2, 2 pi/3]. All in (0, pi), so |sin| = sin. *)

(* Sin is decreasing on [pi/2, 2pi/3], and pi/3 + L/2 <= pi/2 iff L <= pi/3.
   For L > pi/3: pi/3 + L/2 > pi/2, so sin is decreasing throughout. Max at
   smallest argument = pi/3 + L/2, occurring at beta = 2 pi/3. *)

maxChordR3AtEndpoint[L_] := 2 Sin[Pi/3 + L/2];
(* This is the chord from R_3 to beta = 2 pi/3 (approached). *)

(* Now check: is R_3 really the argmax among the 3 robots at this beta? *)

(* At beta = 2 pi/3: *)
(*  R_1 = e(L): chord = 2 sin((2pi/3 - L)/2) = 2 sin(pi/3 - L/2) *)
(*  R_2 = e(2pi/3 + L): chord = 2 sin(L/2) *)
(*  R_3 = e(4pi/3 + L): chord = 2 sin(pi/3 + L/2) *)

threeChordsAtEndpoint[L_] := {
  2 Sin[Pi/3 - L/2],
  2 Sin[L/2],
  2 Sin[Pi/3 + L/2]
};

(* For L in (pi/3, 2pi/3), which is max? *)

Print["--- 3-fold symmetric, all-scan-L, evaluated at unexplored-arc endpoint ---"];
Do[
  Print["L = ", N[L0, 6], ":   chords = ",
    N[threeChordsAtEndpoint[L0], 5],
    "   max = ", N[Max @@ threeChordsAtEndpoint[L0], 5]],
  {L0, {Pi/3, 4 Pi/9, Pi/2, 0.6 Pi, 2 Pi/3}}
];

(* At L = Pi/3: chords = {sin(pi/6)=0.5, 2sin(pi/6)=1, 2sin(pi/2)=2}. Max = 2.
   At L = 2Pi/3: chords = {0, sqrt(3), sqrt(3)}. Max = sqrt(3) ~ 1.732.
   At L = 1.274: chords ~ {0.796, 1.189, 1.988}. Max = 1.988 (R_3 wins). *)

(* Now check: is the max attained at beta = 2pi/3 or at another beta? *)

Print[];
Print["--- Sweep beta in (L, 2pi/3) at L = 1.274 ---"];
Lfix = 1.273754617454;
Do[
  beta = Lfix + (2 Pi/3 - Lfix) * frac;
  chords = Table[chordAB[R[k, Lfix], beta], {k, 1, 3}];
  Print["frac=", N[frac, 3], "  beta=", N[beta, 5],
    "   chords = ", N[chords, 5], "  max = ", N[Max @@ chords, 5]],
  {frac, Range[0.0, 1.0, 0.1]}
];

(* So we see: as beta increases from L to 2 pi/3:
    - Chord from R_3 increases monotonically.
    - Chord from R_2 decreases monotonically.
    - Chord from R_1 first increases then decreases (or monotone).
   The sup is at beta -> 2 pi/3, max = 2 sin(pi/3 + L/2). *)

Print[];
Print["--- Conclusion for 3-fold symmetric ---"];
Print["For L in (pi/3, 2pi/3), sup over beta in (L, 2pi/3) of max_j chord"];
Print["  = lim_{beta -> 2pi/3^-} chord from R_3 to beta"];
Print["  = 2 sin(pi/3 + L/2)."];
Print[];
Print["At L = x_CGK = (2/3) arccos(-1/3) ~ 1.2738:"];
xcgk = (2/3) ArcCos[-1/3];
Print["  sup chord = ", N[2 Sin[Pi/3 + xcgk/2], 10]];
Print["  refined LB = 1 + x + sup chord = ", N[1 + xcgk + 2 Sin[Pi/3 + xcgk/2], 10]];
Print["  CGK LB                          = ", N[1 + xcgk + 2 Sin[3 xcgk / 2], 10]];
Print["  improvement                     = ",
  N[(1 + xcgk + 2 Sin[Pi/3 + xcgk/2]) - (1 + xcgk + 2 Sin[3 xcgk / 2]), 10]];

(* So for 3-fold symmetric full-scan, refined LB = 1 + x + 2 sin(pi/3 + x/2).
   Compare to CGK's universal = 1 + x + 2 sin(3x/2). *)

(* Which config is the ALGORITHM's optimum (minimum of sup)? *)
(* Claim: among all "3-robot on boundary with |E| = 3x" configurations,
   the min of "sup over unexplored of max chord" is achieved by the
   3-fold symmetric. This should give min = 2 sin(pi/3 + x/2). *)

(* To prove: for any 3 robots on boundary with 3 disjoint scan arcs
   of total length 3x ending at robot positions, there exists an
   unexplored boundary point p with chord to the farthest robot
   >= 2 sin(pi/3 + x/2). *)

(* This is a geometric lemma. Let me check numerically against several
   alternative configs. *)

Print[];
Print["--- Alternative configs at |E| = 3x, x = 1.274 ---"];

(* Single-arc scan: arcs [0, x], [x, 2x], [2x, 3x] contiguous; robots
   at x, 2x, 3x. *)
singleArcMaxChord[x_] := Module[{robots, grid, results, betaOpt, maxC},
  robots = {x, 2 x, 3 x};
  (* Adversary picks beta in (3x, 2pi) *)
  grid = Range[3 x + 0.0001, 2 Pi - 0.0001, (2 Pi - 3 x - 0.0002)/500];
  results = Table[
    {beta, Max @@ Table[2 Sin[Abs[Mod[beta - r + Pi, 2 Pi] - Pi]/2], {r, robots}]},
    {beta, grid}
  ];
  maxC = Max @@ results[[All, 2]];
  betaOpt = SortBy[results, -#[[2]] &][[1, 1]];
  {maxC, betaOpt}
];

{mc, bo} = singleArcMaxChord[xcgk];
Print["Single-arc [0, 3x], robots at x, 2x, 3x:"];
Print["  max chord = ", N[mc, 5], " at beta = ", N[bo, 5]];

(* 3-fold symmetric, CCW, L = x *)
threeFoldMaxChord[x_] := Module[{robotAngles, grid, results, maxC, betaOpt},
  robotAngles = Table[2 Pi (j - 1)/3 + x, {j, 1, 3}];
  (* Unexplored: 3 arcs. By symmetry, sample one: (x, 2 pi/3). *)
  grid = Range[x + 0.0001, 2 Pi/3 - 0.0001, (2 Pi/3 - x - 0.0002)/500];
  results = Table[
    {beta, Max @@ Table[2 Sin[Abs[Mod[beta - r + Pi, 2 Pi] - Pi]/2], {r, robotAngles}]},
    {beta, grid}
  ];
  maxC = Max @@ results[[All, 2]];
  betaOpt = SortBy[results, -#[[2]] &][[1, 1]];
  {maxC, betaOpt}
];

{mc2, bo2} = threeFoldMaxChord[xcgk];
Print["3-fold symmetric, CCW scan of length x each:"];
Print["  max chord = ", N[mc2, 5], " at beta = ", N[bo2, 5]];

(* Does the formula match? *)
Print["  Predicted formula 2 sin(pi/3 + x/2) = ", N[2 Sin[Pi/3 + xcgk/2], 5]];
