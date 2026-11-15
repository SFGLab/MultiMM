#!/usr/bin/env bash
# Compile every figure to ../docs/assets/tikz/<name>.svg (needs pdflatex + pdftocairo)
set -e
cd "$(dirname "$0")"
mkdir -p ../docs/assets/tikz .build
for f in fig_*.tex; do
  n="${f%.tex}"
  pdflatex -interaction=nonstopmode -halt-on-error -output-directory .build "$f" >/dev/null
  pdftocairo -svg ".build/$n.pdf" "../docs/assets/tikz/$n.svg"
  echo "built $n"
done
