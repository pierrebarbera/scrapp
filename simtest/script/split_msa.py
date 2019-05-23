#!/usr/bin/env python

# Copyright (C)2019 Pierre Barbera

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Contact:
# Pierre Barbera <Pierre.Barbera@h-its.org>
# Exelixis Lab, Heidelberg Institute for Theoretical Studies
# Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany

import os
import sys
import json

if len(sys.argv) != 2:
  raise Exception("Usage: {} <run-dir>", sys.sys.argv[0])

base_dir = sys.argv[1]

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
