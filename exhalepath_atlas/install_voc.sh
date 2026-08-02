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
echo "==> Ready. Easy usage:"
echo '    voc patient "35M with schizophrenia, smokes, on olanzapine, BMI 32"'
echo '    voc stack "depression" -l brain -c obesity --age 24 --sex male'
echo '    voc "depression" -l brain -c obesity --age 24 --sex male'
echo '    voc eval-stack-holdout'
echo '    voc --help'
echo 'Docs: great_disease_stack/INSTALL.md · great_disease_stack/SOTA.md'
