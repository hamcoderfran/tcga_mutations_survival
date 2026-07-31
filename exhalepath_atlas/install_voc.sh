#!/usr/bin/env bash
# One-liner installer for voc (ExhalePath Atlas)
# Usage:
#   curl -fsSL <raw-url>/install_voc.sh | bash
#   # or
#   bash install_voc.sh
set -euo pipefail

REPO="${VOC_REPO:-https://github.com/hamcoderfran/tcga_mutations_survival.git}"
SUBDIR="${VOC_SUBDIR:-exhalepath_atlas}"
BRANCH="${VOC_BRANCH:-main}"

echo "==> Installing voc-breath (command: voc)"
python3 -m pip install --upgrade pip
python3 -m pip install "voc-breath @ git+${REPO}@${BRANCH}#subdirectory=${SUBDIR}"

echo
echo "==> Ready. One line → full visual report:"
echo '    voc "depression" -l brain -c obesity --age 24 --sex male'
echo '    # → terminal panel + runs/voc_*/REPORT.html + dashboard.png'
echo '    voc "lung adenocarcinoma" -l "left lower lobe" --stage II --genes KRAS,TP53'
echo '    voc list-diseases'
echo '    voc --help'
