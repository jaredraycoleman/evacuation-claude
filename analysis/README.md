# analysis/

Wolfram Language scripts for symbolic / analytical work — closed-form optima,
derivatives of worst-case exit-time functions, verifying the
$2\cos\beta + \cos\gamma = 1$ condition from the F2F paper, etc.

Run with `wolframscript -file <name>.wls`. Files should be self-contained and
print their main result to stdout so it can be diffed between runs.

> **Note:** verify `wolframscript -version` works before starting; if not
> available, fall back to `sympy` in `experiments/` and record the
> substitution in `CLAUDE.md`.
