#!/usr/bin/env bash
# Copy the compiled report PDF into site/ so static hosts serve it directly.
# Intended as the build command for Cloudflare Pages (or any static host).
#
# Usage (from repo root):
#   ./site/build.sh
#
# Cloudflare Pages configuration:
#   Build command:          ./site/build.sh
#   Build output directory: site
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Rebuild the paper if pdflatex is available; otherwise use the committed PDF.
if command -v pdflatex >/dev/null 2>&1; then
  (cd "$REPO_ROOT/paper" && ./build.sh)
fi

cp "$REPO_ROOT/paper/main.pdf" "$SCRIPT_DIR/report.pdf"
echo "copied report.pdf ($(du -h "$SCRIPT_DIR/report.pdf" | cut -f1))"
