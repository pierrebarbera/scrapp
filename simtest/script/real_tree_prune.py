#!/usr/bin/env python2

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

import argparse
import os
import sys
import json
import numpy.random as rd
from ete3 import Tree, PhyloTree
from collections import defaultdict
import pprint as pp
import math

def command_line_args_parser():
    """
    Return an instance of argparse that can be used to process command line arguemnts.
    """

    def unit_interval(x):
        x = float(x)
        if x < 0.0 or x > 1.0:
            raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
        return x

    # Init an args parser, with a group of required named arguments. It is just nicer to use named
    # arguments than having to rely on their order (i.e., use positional arguments instead).
    parser = argparse.ArgumentParser(
        description="Script that uses msprime to generate trees to be used with the scrapp sim stack"
    )

    parser.add_argument(
        '--out-dir',
        help="Output directory.",
        action='store',
        dest='out_dir',
        type=str,
        default="."
    )

    parser.add_argument(
        '--tree',
        help="The 'true' tree. (newick)",
        action='store',
        dest='tree',
        type=str,
        required=True
    )

    parser.add_argument(
        '--msa',
        help="The corresponding MSA to be split into ref and query set. (fasta)",
        action='store',
        dest='msa',
        type=str,
        required=True
    )


    # Add optional args.
    parser.add_argument(
        '--seed',
        help="Seed for the RNG.",
        action='store',
        dest='seed',
        type=int
    )

    parser.add_argument(
        '--rarify',
        help="How many (as a fraction of total count) taxa to prune from the ref tree to simulate missing data.",
        action='store',
        dest='rarify_fract',
        type=unit_interval,
        default=0.2
    )

    parser.add_argument(
        '--query',
        help="How many (as a fraction of total count) taxa to prune from the tree and set aside as queries.",
        action='store',
        dest='query_fract',
        type=unit_interval,
        default=0.2
    )

    return parser

def command_line_args_postprocessor( args ):
    # Make sure that all paths are fully resolved and dirs have no trailing slashes.
    # args.jplace_file = os.path.abspath( os.path.realpath( args.jplace_file ))

    return args

def command_line_args():
    """
    Return a parsed and processed list of the command line arguments that were provided when
    running this script.
    """

    # Parse the given arguments from the command line, post-process them, return the result.
    parser = command_line_args_parser()
    args = parser.parse_args()
    args = command_line_args_postprocessor( args )
    return args

def path_between(A, B):
    # ETE3 doesn't have a built in path function (AFAIK), so instead
    # we combine the paths to a common ancestor:
    ancestor, path = A.get_common_ancestor( B, get_path=True )

    # apparently ETE returns path as the path of A and B to the root here...
    # so we need to find first where the paths intersect
    AtoRoot = path[A]
    BtoRoot = path[B]

    # symmetric set difference
    diffset = AtoRoot ^ BtoRoot

    # take the difference of A's and B's path to root and add the ancestor
    return diffset | {ancestor}


def get_node_by_name(tree, name):
    nodes = tree.search_nodes(name=name)
    if not nodes:
        raise Exception("No such node: {}", name)
    return nodes[0]

def popcount( taxa, pop_map ):
    taxonset = set(taxa)
    vset = set(pop_map.values())

    return len( vset & taxonset )

def print_species_map( spec_map, file ):
    with open(file, "w+") as f:
        f.write(json.dumps(spec_map, sort_keys=True,
                indent=4, separators=(',', ': ')))

def rand_names(n=1):
    word_file = os.path.join(dir_path, "words.txt")
    words = open(word_file).read().splitlines()

    rand_nums = rd.choice(len(words), n, replace=False)
    return {i:words[rand_nums[i]] for i in range(len(rand_nums))}


# ===========================
#  MAIN FUNCTION
# ===========================

args = command_line_args()

out_dir = args.out_dir

dir_path = os.path.dirname(os.path.realpath(__file__))

if args.seed:
    SEED=int(args.seed)
else:
    SEED=None
rd.seed(SEED)

RARIFY_FRACT=args.rarify_fract
QUERY_FRACT=args.query_fract

# read in the tree and alignment
true_tree = Tree( args.tree )


#### give inner nodes unique names such that we can find a mapping between true and ref later
number=1000
for node in true_tree.iter_search_nodes():
    if not node.name:
        node.name = str( number )
        number += 1

# decide what to prune
# randomly prune a fraction of the taxa

# first build a dict of all taxa
ref_map = dict()
qry_map = dict()
i=0
for leaf in true_tree.iter_leaves():
    ref_map[ i ] = leaf.name
    i+=1

# then randomly reassign some of them as query
N = len( ref_map )
x=10
k=int( QUERY_FRACT * N ) / x
seeds = rd.choice( ref_map.keys(), size=k, replace=False )


# after we found k seeds, we extend those by a random radius limited by x
total=0
for key in seeds:
    # print "Seed:   ",key
    # draw a random radius between 2 and x / 2
    radius = rd.randint(2, (x/2) + 1)
    # print "Radius: ",radius
    advance_right=True
    advance_left=True
    i=1
    taken_keys=[key]
    qry_map[ key ] = ref_map[ key ]
    del ref_map[ key ]
    while i < radius and i > 0 and advance_right and advance_left and i < N:
        rh_key=key+i
        if advance_right and rh_key in ref_map and rh_key+1 not in seeds and rh_key not in seeds:
            assert( rh_key not in qry_map )
            taken_keys.append( rh_key )
            qry_map[ rh_key ] = ref_map[ rh_key ]
            del ref_map[ rh_key ]
        else:
            advance_right=False

        lh_key=key-i
        if advance_left and lh_key in ref_map and lh_key-1 not in seeds and lh_key not in seeds:
            assert(  lh_key not in qry_map )
            taken_keys.append( lh_key )
            qry_map[ lh_key ] = ref_map[ lh_key ]
            del ref_map[ lh_key ]
        else:
            advance_left=False
        i+=1
    # print sorted(taken_keys)
    total += len(taken_keys)

# print total, " vs. ", len(qry_map)
assert( total == len(qry_map) )
assert( len(qry_map)+len(ref_map) == N )

# Further rarify the reference tree to simulate missing data
N = len( ref_map )
num_remove_rarify = 0
for key in rd.choice( ref_map.keys(), size=int(N*RARIFY_FRACT), replace=False):
    del ref_map[key]
    num_remove_rarify+=1

# In the tree, strip anything but the selected reference taxa
ref_list=[v for k,v in ref_map.iteritems()]
assert( len(ref_list) == len(ref_map) )
# print ref_list
ref_tree = true_tree.copy()
ref_tree.prune( ref_list, preserve_branch_length=True )

# print len( ref_tree )
assert( len(ref_map) == len(ref_tree) )
# print len( true_tree )

#### write the ref_tree
with open(os.path.join(out_dir, "reference.newick"), "w+") as f:
    f.write(ref_tree.write(format=5,dist_formatter="%.12f"))

# write the reference and query fasta files
from Bio import SeqIO

with open( os.path.join( out_dir, "reference.fasta"), "w+" ) as ref_msa:
    with open( os.path.join( out_dir, "query.fasta"), "w+" ) as qry_msa:
        for seq in SeqIO.parse( args.msa, 'fasta' ):
            if seq.name in ref_map.values():
                SeqIO.write( seq, ref_msa, 'fasta')
            elif seq.name in qry_map.values():
                SeqIO.write( seq, qry_msa, 'fasta')


# print ref_tree
total=0
#### determine the "true" species count of the edges in the ref tree
# we do this by getting every edge (tuple of the two adjacent nodes) in the ref tree
for distal in ref_tree.iter_search_nodes():
    proximal = distal.up
    # (except the root, its two edges are covered bottom up already)
    if not proximal:
        continue

    # ... and finding their original locations in the true tree
    P = get_node_by_name(true_tree, proximal.name)
    D = get_node_by_name(true_tree, distal.name)

    # then getting the path between those two in the true tree
    prune_path = path_between(P, D)

    # all nodes n from the path, except
    pruned_nodes = [y for n in prune_path if n not in [P,D] for x in n.children if x not in prune_path for y in x.get_leaves()]

    # and making a list of all leaves' labels in that list
    pruned_leaf_labels = [ node.name for node in pruned_nodes if node.is_leaf() ]

    # and using that to look up how many populations were removed
    count = popcount( pruned_leaf_labels, qry_map )
    # count = len( pruned_leaf_labels )
    total += count

    # finally we capture that information on the ref tree
    distal.add_features(species_count=count)

# print total, " vs. ", len( true_tree ) - len( ref_tree )
assert( total + num_remove_rarify == len( true_tree ) - len( ref_tree ) )

#### write the annotated ref_tree
with open(os.path.join(out_dir, "annot_reference.newick"), "w+") as f:
    f.write( ref_tree.write( format=5,dist_formatter="%.12f",features=["species_count"] ) )
