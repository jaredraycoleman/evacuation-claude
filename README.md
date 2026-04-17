# evacuation-claude

Research scratchpad for open problems in robot evacuation from the disk (unknown-exit, $k$-robot, F2F / wireless / SR / priority variants).

See [`CLAUDE.md`](CLAUDE.md) for the current status, active problem, and workflow.

## Layout

- `paper/` — LaTeX writeup; build with `paper/build.sh` (runs pdflatex + bibtex + pdflatex twice).
- `experiments/` — Python simulators + numerical optimisers (via `uv`).
- `analysis/` — Wolfram Language scripts for symbolic work.
- `notes/` — literature notes.
- `site/` — small static website with an interactive A₃ visualiser and the compiled report.

## Deploying the site to Cloudflare Pages

The `site/` directory is a static site (vanilla HTML/CSS/JS, no framework) that
ships the compiled report PDF alongside an interactive canvas visualisation of
the $A_3$ algorithm.

To deploy via Cloudflare Pages:

1. Connect the GitHub repo (`jaredraycoleman/evacuation-claude`) in the
   Cloudflare Pages dashboard.
2. Set the build configuration:
   - **Build command:** `./site/build.sh`
   - **Build output directory:** `site`
   - **Root directory:** *(leave as repo root)*
3. Pages will run `./site/build.sh`, which rebuilds the paper (if
   `pdflatex` is available in the build image) and copies `paper/main.pdf`
   into `site/report.pdf`; then serve the contents of `site/`.

For local preview:
```
./site/build.sh            # produces site/report.pdf
python3 -m http.server --directory site 8000
# browse to http://localhost:8000/
```
