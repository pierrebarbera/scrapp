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
# set -e
echo "start at `date`"
for thresh in {50..300..50}; do
  echo "${thresh}: "
  for run in {0..4}; do
    printf "${run} "
    for scrapp_mode in rootings bootstrap outgroup; do
      OUT=${BASE}/scrapp/${thresh}/${run}/${scrapp_mode}
      mkdir -p $OUT
      rm -r $OUT/* 2> /dev/null

      case "${scrapp_mode}" in
        rootings )
          ~/scrapp/scrapp.py --cluster-above ${thresh} --jplace ${JP} --alignment ${QRY} --work-dir ${OUT} --min-weight 0.5 --num-threads 40 --verbose "$@" > ${OUT}/scrapp_log.txt
          ;;
        bootstrap )
          ~/scrapp/scrapp.py --cluster-above ${thresh} --jplace ${JP} --alignment ${QRY} --work-dir ${OUT} --min-weight 0.5 --num-threads 40 --verbose "$@" --bootstrap > ${OUT}/scrapp_log.txt
          ;;
        outgroup )
          ~/scrapp/scrapp.py --cluster-above ${thresh} --jplace ${JP} --alignment ${QRY} --work-dir ${OUT} --min-weight 0.5 --num-threads 40 --verbose "$@" --ref-align-outgrouping ${DATA}/ref.fasta > ${OUT}/scrapp_log.txt
          ;;
        *)
          echo "invalid scrapp_mode, aborting"
          exit 1
      esac
    done
  done
  echo ""
done
echo "end at `date`"
