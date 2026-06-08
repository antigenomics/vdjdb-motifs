#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TCRNET_NPROC="${TCRNET_NPROC:-8}"
MOTIF_PWM_HUMAN_TRA_BACKGROUND="${MOTIF_PWM_HUMAN_TRA_BACKGROUND:-tcrnet/pwms/human.tra.aa.aa_cdr3_pwm.txt}"
MOTIF_PWM_HUMAN_TRB_BACKGROUND="${MOTIF_PWM_HUMAN_TRB_BACKGROUND:-tcrnet/pwms/human.trb.aa.aa_cdr3_pwm.txt}"
MOTIF_PWM_MOUSE_TRA_BACKGROUND="${MOTIF_PWM_MOUSE_TRA_BACKGROUND:-tcrnet/pwms/mouse.tra.aa.aa_cdr3_pwm.txt}"
MOTIF_PWM_MOUSE_TRB_BACKGROUND="${MOTIF_PWM_MOUSE_TRB_BACKGROUND:-tcrnet/pwms/mouse.trb.aa.aa_cdr3_pwm.txt}"

cd "$ROOT_DIR"
export TCRNET_NPROC
Rscript -e "rmarkdown::render('tcrnet/compute_vdjdb_motifs.Rmd')"
python scripts/generate_tcrnet_reports.py
python scripts/compute_motif_pwms.py \
  --cluster-members results/tcrnet/cluster_members.txt \
  --background HomoSapiens TRA "$MOTIF_PWM_HUMAN_TRA_BACKGROUND" \
  --background HomoSapiens TRB "$MOTIF_PWM_HUMAN_TRB_BACKGROUND" \
  --background MusMusculus TRA "$MOTIF_PWM_MOUSE_TRA_BACKGROUND" \
  --background MusMusculus TRB "$MOTIF_PWM_MOUSE_TRB_BACKGROUND" \
  --output results/tcrnet/motif_pwms.txt
