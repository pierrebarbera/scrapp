#!/bin/bash
# source common.sh

BASE=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
DATA=${BASE}/data

# TD=${DATA}/tree/${REVIS}/cons
# REF=${DATA}/reference/${REVIS}/2_REF_align_complete_for_Papara_EPAng.fas
REF=${DATA}/ref.fasta
QRY=${DATA}/1k.fasta
TREE=${DATA}/tree.newick
JP=${BASE}/place/epa_result.jplace

RESULT_CSV=${BASE}/results_threshold.csv

CSP=~/scrapp/simtest/script/compare_species_counts
# CSP=~/pbgenesis/bin/apps/compare_species_counts

# set -e
echo "run,thresh,scrapp_mode,krd,norm_krd,norm_norm_krd,norm_unit_krd,abs_err_sum,abs_err_mean,abs_err_med,rel_err_sum,rel_err_mean,rel_err_med,norm_err_sum,norm_err_mean,norm_err_med,rel_norm_err_mean" > ${RESULT_CSV}
echo "start at `date`"
for thresh in {50..300..50}; do
  echo ${thresh}
  for run in {0..4}; do
    for scrapp_mode in rootings bootstrap outgroup; do
      OUT=${BASE}/scrapp/${thresh}/${run}/${scrapp_mode}
      BENCH=${BASE}/scrapp/500/0/${scrapp_mode}
      printf "${run},${thresh},${scrapp_mode}," >> ${RESULT_CSV}
      ${CSP} ${OUT}/summary.newick ${BENCH}/summary.newick 1 >> ${RESULT_CSV}
      printf "\n" >> ${RESULT_CSV}
    done
  done
done
echo "end at `date`"

cat ${RESULT_CSV}
