#!/usr/bin/env python

import os
import sys
import json

script_dir = os.path.dirname(os.path.realpath(__file__))

base_dir = os.path.realpath(os.path.join(script_dir, "../"))

out_dir=os.path.join(base_dir, "msa")



#### read ref and qry maps
with open(os.path.join(base_dir, "tree/refmap"), "r") as f:
    ref_map=json.load(f)
with open(os.path.join(base_dir, "tree/qrymap"), "r") as f:
    qry_map=json.load(f)
# print(ref_map)
# print(qry_map)

ref_list=[v[0] for k,v in ref_map.iteritems()]
qry_list=sorted({x for v in qry_map.itervalues() for x in v})

# print(ref_list)
# print(qry_list)

#### read MSA and split into ref and qry
from Bio import SeqIO

msa_in=os.path.join(out_dir, "full_TRUE.phy")
ref_out=os.path.join(out_dir,"reference.fasta")
qry_out=os.path.join(out_dir,"query.fasta")
with open(msa_in, "rU") as msa_f, open(ref_out, "w+") as ref_f, open(qry_out, "w+") as qry_f:
    for seq in SeqIO.parse(msa_f, "phylip-relaxed"):
        if seq.id in ref_list:
            SeqIO.write(seq, ref_f, "fasta")
        if seq.id in qry_list:
            SeqIO.write(seq, qry_f, "fasta")
