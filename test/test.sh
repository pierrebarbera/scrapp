#!/bin/bash

BASE=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
DATA=${BASE}/data

REF=${DATA}/ref.fasta
QRY=${DATA}/1k.fasta
TREE=${DATA}/tree.newick
JP=${BASE}/place/epa_result.jplace

cd ${BASE}
mkdir -p out
rm -rf out/*

../scrapp.py --jplace ${JP} --alignment ${QRY} --work-dir out --min-weight 0.5 "$@"

