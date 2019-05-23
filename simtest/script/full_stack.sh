#!/bin/bash
DATE=`date '+%Y-%m-%d-%H:%M'`

BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)
REF=${BASE}/msa/reference.fasta
QRY=${BASE}/msa/query.fasta
TREE=${BASE}/tree/reference.newick
MODEL=${BASE}/tree/eval.raxml.bestModel
JPLACE=${BASE}/placed/epa_result.jplace
SC=${BASE}/script

NUM_THREADS=4
SEED=4242

set -e

echo "start at `date`"

cd ${SC}

# write the csv header
echo "run,mut_rate,species,pop_size,prune_fract,krd,norm_krd,norm_norm_krd,norm_unit_krd,abs_err_sum,abs_err_mean,abs_err_med,rel_err_sum,rel_err_mean,rel_err_med,norm_err_sum,norm_err_mean,norm_err_med" > results.csv

run=0

for prune_fract in 0.2 0.5; do
  for pop_size in 1e5 1e6 1e7; do
    for species in 20; do
      for mut_rate in 1e-5 1e-6 1e-7; do
        for i in {1..2}; do
          echo "Starting run ${run}!"

          SCRAPP_SIM_CURDIR=${BASE}/runs/${run}
          export SCRAPP_SIM_CURDIR

          printf "${run},${mut_rate},${species},${pop_size},${prune_fract}," >> results.csv

          # echo "  generate the tree..."
          ./msprime.sh --mutation-rate ${mut_rate} --species ${species} --population-size ${pop_size} --prune ${prune_fract} 1> /dev/null
          # echo "  tree done!"

          # generate the sequences and split into query and ref set
          # echo "  generate the sequences..."
          ./seqgen.sh 1> /dev/null
          # echo "  sequences done!"

          # infer model params
          # echo "  infer model params..."
          #  --blmin 1e-7 --blmax 5000
          ./eval_reftree.sh --threads ${NUM_THREADS} --seed --opt-branches off 1> /dev/null
          # echo "  model params done!"

          # run placement
          # echo "  place..."
          ./epa.sh --threads ${NUM_THREADS} 1> /dev/null
          # echo "  placement done!"

          # run scrapp
          # echo "  running scrapp..."
          # ./scrapp.sh --num-threads ${NUM_THREADS} --seed ${SEED} 1> /dev/null
          # ./scrapp.sh --num-threads ${NUM_THREADS} --seed ${SEED} --bootstrap 1> /dev/null
          ./scrapp.sh --num-threads ${NUM_THREADS} --ref-align-outgrouping ${SCRAPP_SIM_CURDIR}/msa/reference.fasta 1> /dev/null
          # echo "  scrapp done!"

          # print statistic
          # echo "  printing statistic..."
          ./compare_species_counts ${SCRAPP_SIM_CURDIR}/delimit/summary.newick ${SCRAPP_SIM_CURDIR}/tree/annot_reference.newick 1 >> results.csv
          # echo "  statistic done!"

          let run+=1

        done # runs
      done # mut_rate
    done # species
  done # pop_size
done


echo "end at `date`"
