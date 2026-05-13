#!/bin/zsh
set -euo pipefail

cd "$(dirname "$0")"

export TOPAS_G4_DATA_DIR=/Applications/G4Data
export QT_QPA_PLATFORM_PLUGIN_PATH=/Applications/topas/Frameworks
export PATH="/Applications/topas/bin:$PATH"

TOPAS="/Applications/topas/bin/topas"
MATLAB="/Applications/MATLAB_R2025a.app/bin/matlab"

section() {
  echo ""
  echo "============================================================"
  echo "$1"
  echo "============================================================"
}

require_file() {
  if [[ ! -f "$1" ]]; then
    echo "ERROR: Missing required file: $1"
    exit 1
  fi
}

if [[ ! -x "$TOPAS" ]]; then
  echo "ERROR: TOPAS executable not found at $TOPAS"
  exit 1
fi

require_file "UITF2-10M-front.txt"
require_file "UITF2-10M-back.txt"
require_file "analyze_DNP_UITF2.py"
require_file "analyze_DNP_UITF2_MATLAB.m"

section "UITF2 v2 TOPAS production: front side, 10M histories"
"$TOPAS" UITF2-10M-front.txt

section "UITF2 v2 TOPAS production: back side, 10M histories"
"$TOPAS" UITF2-10M-back.txt

section "Python analysis: front-side 10M output"
python3 analyze_DNP_UITF2.py --output_dir UITF_DNP_Output_10M_front --mode cold --material ND3

section "Python analysis: back-side 10M output"
python3 analyze_DNP_UITF2.py --output_dir UITF_DNP_Output_10M_back --mode cold --material ND3

section "MATLAB analysis: combined front/back thesis figures"
if [[ -x "$MATLAB" ]]; then
  "$MATLAB" -batch "cd('/Users/takwirira/Desktop/TOPAS/UITF02'); analyze_DNP_UITF2_MATLAB"
else
  echo "WARNING: MATLAB executable not found at $MATLAB"
  echo "Skipping MATLAB combined front/back figure generation."
fi

section "UITF2 v2 pipeline complete"
echo "TOPAS outputs:"
echo "  UITF_DNP_Output_10M_front"
echo "  UITF_DNP_Output_10M_back"
echo "Python figures:"
echo "  UITF_DNP_Output_10M_front/figures"
echo "  UITF_DNP_Output_10M_back/figures"
echo "MATLAB figures:"
echo "  figures_matlab"
