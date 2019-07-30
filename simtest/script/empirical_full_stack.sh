#!/bin/bash
DATE=`date '+%Y-%m-%d-%H:%M'`

BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)
REF=${BASE}/msa/reference.fasta
QRY=${BASE}/msa/query.fasta
TREE=${BASE}/tree/reference.newick
MODEL=${BASE}/tree/eval.raxml.bestModel
JPLACE=${BASE}/placed/epa_result.jplace
SC=${BASE}/script

NUM_THREADS=40
SEED=4242

RESULT_CSV="empirical_results.csv"

echo "start at `date`"

cd ${SC}

# set -e
# write the csv header
# echo "run,seq_length,mut_rate,species,sample_size,pop_size,prune_fract,krd,norm_krd,norm_norm_krd,norm_unit_krd,abs_err_sum,abs_err_mean,abs_err_med,rel_err_sum,rel_err_mean,rel_err_med,norm_err_sum,norm_err_mean,norm_err_med" > ${RESULT_CSV}
echo "run,scrapp_mode,query_fract,rarify_fract,krd,norm_krd,norm_norm_krd,norm_unit_krd,abs_err_sum,abs_err_mean,abs_err_med,rel_err_sum,rel_err_mean,rel_err_med,norm_err_sum,norm_err_mean,norm_err_med" > ${RESULT_CSV}

run=0
for scrapp_mode in rootings bootstrap outgroup; do
  for query_fract in 0.2 0.4 0.6; do
    for rarify_fract in 0.2 0.4 0.6; do
      for i in {1..10}; do
        echo "Starting run ${run}!"

        # scrapp_mode=rootings

        SCRAPP_SIM_CURDIR=${BASE}/runs/empirical/scrapp_mode_${scrapp_mode}/query_fract_${query_fract}/rarify_fract_${rarify_fract}/iter_${i}
        export SCRAPP_SIM_CURDIR
        rm -rf ${SCRAPP_SIM_CURDIR}/* 2> /dev/null

        printf "${run},${scrapp_mode},${query_fract},${query_fract}," >> ${RESULT_CSV}

        # echo "  generate the tree..."
        ./split_empirical.sh --rarify ${rarify_fract} --query ${query_fract}
        # echo "  tree done!"

        # infer model params
        # echo "  infer model params..."
        #  --blmin 1e-7 --blmax 5000
        ./eval_reftree.sh --threads ${NUM_THREADS} --seed --opt-branches off --force perf 1> /dev/null
        # echo "  model params done!"

        # run placement
        # echo "  place..."
        ./epa.sh --threads ${NUM_THREADS} 1> /dev/null
        # echo "  placement done!"

        # run scrapp
        # echo "  running scrapp..."
        case "${scrapp_mode}" in
          rootings )
            ./scrapp.sh --num-threads ${NUM_THREADS} --seed ${SEED} 1> /dev/null
            ;;
          bootstrap )
            ./scrapp.sh --num-threads ${NUM_THREADS} --seed ${SEED} --bootstrap 1> /dev/null
            ;;
          outgroup )
            ./scrapp.sh --num-threads ${NUM_THREADS} --ref-align-outgrouping ${SCRAPP_SIM_CURDIR}/msa/reference.fasta 1> /dev/null
            ;;
          *)
            echo "invalid scrapp_mode, aborting"
            exit 1
        esac
        # echo "  scrapp done!"

        # print statistic
        # echo "  printing statistic..."
        ./compare_species_counts ${SCRAPP_SIM_CURDIR}/delimit/summary.newick ${SCRAPP_SIM_CURDIR}/tree/annot_reference.newick 1 >> ${RESULT_CSV}
        printf "\n" >> ${RESULT_CSV}
        # echo "  statistic done!"

        let run+=1

      done # runs
    done # rarify_fract
  done # query_fract
done # scrapp_mode
echo "end at `date`"
