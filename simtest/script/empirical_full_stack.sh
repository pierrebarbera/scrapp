#!/bin/bash
DATE=`date '+%Y-%m-%d-%H:%M'`

BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)
REF=${BASE}/msa/reference.fasta
QRY=${BASE}/msa/query.fasta
TREE=${BASE}/tree/reference.newick
MODEL=${BASE}/tree/eval.raxml.bestModel
JPLACE=${BASE}/placed/epa_result.jplace
SC=${BASE}/script
RUNS_DIR=/data/barberpe/scrapp/runs


NUM_THREADS=40
SEED=4242

RESULT_CSV="results_empirical.csv"

echo "start at `date`"

cd ${SC}

# set -e

# write the csv header
echo "run,query_fract,rarify_fract,scrapp_mode,krd,norm_krd,norm_norm_krd,norm_unit_krd,abs_err_sum,abs_err_mean,abs_err_med,rel_err_sum,rel_err_mean,rel_err_med,norm_err_sum,norm_err_mean,norm_err_med,rel_norm_err_mean" > ${RESULT_CSV}

run=0
for query_fract in 0.5; do
  for rarify_fract in 0.25 0.5; do
    for scrapp_mode in rootings bootstrap outgroup; do
      for i in {0..4}; do
        echo "Starting run ${run}!"

        SCRAPP_SIM_CURDIR=${RUNS_DIR}/empirical/query_fract_${query_fract}/rarify_fract_${rarify_fract}/scrapp_mode_${scrapp_mode}/iter_${i}
        export SCRAPP_SIM_CURDIR
        # mkdir -p ${SCRAPP_SIM_CURDIR}
        # rm -r ${SCRAPP_SIM_CURDIR}/* 2> /dev/null

        printf "${run},${query_fract},${rarify_fract},${scrapp_mode}," >> ${RESULT_CSV}

        # # echo "  generate the tree..."
        # ./split_empirical.sh --rarify ${rarify_fract} --query ${query_fract}
        # # echo "  tree done!"

        # # infer model params
        # # echo "  infer model params..."
        # #  --blmin 1e-7 --blmax 5000
        # ./eval_reftree.sh --threads ${NUM_THREADS} --opt-branches off --force perf 1> /dev/null
        # # echo "  model params done!"

        # # run placement
        # # echo "  place..."
        # ./epa.sh --threads ${NUM_THREADS} 1> /dev/null
        # # echo "  placement done!"

        # # run scrapp
        # # echo "  running scrapp..."
        # case "${scrapp_mode}" in
        #   rootings )
        #     ./scrapp.sh --num-threads ${NUM_THREADS} > ${SCRAPP_SIM_CURDIR}/scrapp_log.txt
        #     ;;
        #   bootstrap )
        #     ./scrapp.sh --num-threads ${NUM_THREADS} --bootstrap > ${SCRAPP_SIM_CURDIR}/scrapp_log.txt
        #     ;;
        #   outgroup )
        #     ./scrapp.sh --num-threads ${NUM_THREADS} --ref-align-outgrouping ${SCRAPP_SIM_CURDIR}/msa/reference.fasta > ${SCRAPP_SIM_CURDIR}/scrapp_log.txt
        #     ;;
        #   *)
        #     echo "invalid scrapp_mode, aborting"
        #     exit 1
        # esac
        # echo "  scrapp done!"

        # print statistic
        # echo "  printing statistic..."
        ./compare_species_counts ${SCRAPP_SIM_CURDIR}/delimit/summary.newick ${SCRAPP_SIM_CURDIR}/tree/annot_reference.newick 1 >> ${RESULT_CSV}
        printf "\n" >> ${RESULT_CSV}
        # echo "  statistic done!"

        let run+=1
        # exit
      done # runs
    done # rarify_fract
  done # query_fract
done # scrapp_mode
echo "end at `date`"
