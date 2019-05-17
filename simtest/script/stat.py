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

from __future__ import print_function

import os
import sys
import json
from math import log


script_dir = os.path.dirname(os.path.realpath(__file__))

base_dir = os.path.realpath(os.path.join(script_dir, "../"))

# a function to read in a delimitation result as a cluster indexed list of taxa labels
def read_delim(delim_file, has_mult=True):
    import imp
    mptp = imp.load_source('mptp', '../../scripts/mptp.py')
    return mptp.parse( delim_file, has_mult )['delimitation']

# a function to take the simulated ground truth delimitation and reduce it down to only
# include the taxa specified in a per branch delimitation result
def reduce(lhs, rhs):
    # listify the "to keep" set
    to_keep = sorted({x for v in rhs.itervalues() for x in v})

    from collections import defaultdict
    ret=defaultdict(list)

    for key,taxa_list in lhs.iteritems():
        for taxon in taxa_list:
            if taxon in to_keep:
                ret[key].append(taxon)
    return dict(ret)

def num_taxa(delim):
    num=0
    for taxa_list in delim.itervalues():
        for taxon in taxa_list:
            num += 1
    return num

def len_intersect(lhs, rhs):
    return len([value for value in lhs if value in rhs])

def entropy(part, N):
    ret_sum=0.0

    for p in part.itervalues():
        ret_sum += ( len(p) / N ) * log( len(p) / N )

    return -ret_sum

# a function to take two delimitations and calculate the NMI between them
def NMI(truth, delim):
    N=float(num_taxa(truth))
    I=0.0
    for t in truth.itervalues():
        for d in delim.itervalues():
            intersize = len_intersect(t, d)
            if intersize > 0:
                I += (intersize / N) * log( (intersize * N) / (len(t) * len(d)) )

    return I / max( entropy(truth, N), entropy(delim, N) )


#### read ref and qry maps
with open(os.path.join(base_dir, "tree/popmap"), "r") as f:
    truth=json.load(f)
with open(os.path.join(base_dir, "tree/qrymap"), "r") as f:
    qry_map=json.load(f)

import glob
import numpy as np

NMI_map=dict()

if len(sys.argv) == 2:
    delim=read_delim(sys.argv[1], False)
    print("NMI:",NMI(truth, delim))
else:
    #for each edge
    for edge_dir in glob.glob( os.path.join(base_dir, "delimit/edge_*") ):

        edge_num=int( edge_dir.split("_")[-1] )

        per_edge_NMIs=[]

        # for each variation
        for result_file in glob.glob( os.path.join(edge_dir, "delimit/edge_*/mptp_result.txt") ):
            # parse the result
            result_delim = read_delim( result_file )

            # get a reduced set of the ground truth that only includes the taxa of the result
            alternative_facts = reduce(truth, result_delim)

            # capture the specific result
            per_edge_NMIs.append( NMI(alternative_facts, result_delim) )

        NMI_map[edge_num]=per_edge_NMIs

    print("edge\tmax\tmin\tmean\tmedian\tstddev")
    for k,v in NMI_map.iteritems():
        print(k , format( np.max( v ), '.4f'), format( np.min( v ), '.4f'), format( np.mean( v ), '.4f'), format( np.median( v ), '.4f'), format( np.std( v ), '.4f'), sep='\t')
