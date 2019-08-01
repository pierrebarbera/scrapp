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
from ete3 import Tree
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

    # Add optional args.
    parser.add_argument(
        '--seed',
        help="Seed for the RNG.",
        action='store',
        dest='seed',
        type=int
    )

    parser.add_argument(
        '--mutation-rate',
        help="Mutation rate in mutations per generation per base.",
        action='store',
        dest='mut_rate',
        type=float,
        default=1e-8
    )

    parser.add_argument(
        '-n','--num-pops', '--species',
        help="Number of starting populations (species).",
        action='store',
        dest='num_pops',
        type=int,
        default=30
    )

    parser.add_argument(
        '--sample-size',
        help="Number of individuals per population.",
        action='store',
        dest='sample_size',
        type=int,
        default=10
    )

    parser.add_argument(
        '--population-size',
        help="Total population size (Ne).",
        action='store',
        dest='pop_size',
        type=float,
        default=1e6
    )

    parser.add_argument(
        '--migration-rate',
        help="Migration rate between populations. Warning: strong influence on species count.",
        action='store',
        dest='mig_rate',
        type=float,
        default=1e-11
    )

    parser.add_argument(
        '--seq-length',
        help="Hypothetical length of the underlying sequence/MSA.",
        action='store',
        dest='seq_length',
        type=float,
        default=1000
    )

    parser.add_argument(
        '--prune',
        help="How many (as a fraction of total count) populations/species should be pruned from the tree to emulate missing data.",
        action='store',
        dest='prune_fract',
        type=unit_interval,
        default=0.2
    )

    parser.add_argument(
        '--verbose',
        help="More verbose output.",
        action='store_true'
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


def get_node_by_name(tree, name):
    node = next(tree.iter_search_nodes(name=name))
    if not node:
        raise Exception("No such node: {}", name)
    return node

def print_species_map( spec_map, file ):
    with open(file, "w+") as f:
        f.write(json.dumps(spec_map, sort_keys=True,
                indent=4, separators=(',', ': ')))

def prune_node( tree, ref_map, node ):
    assert node.is_leaf()
    # get the node's species count value (saved as 'support')
    count = node.support + node.up.support

    # decide where to move the species count
    target_id = 1 if node.up.children[0] == node else 0

    # add that support to the neighbor node that will survive the prune
    if not node.up.up and len(node.up.children) == 2:
        # .up is the non-virtual root, that this prune would remove entirely
        # so we have to add the species count further down
        # we decide to add it to the grandchild with the longer branch
        largest_dist_id = 0
        largest_dist = 0.0
        for child_id in range( len( node.up.children[ target_id ].children ) ):
            if node.up.children[ target_id ].children[ child_id ].dist > largest_dist:
                largest_dist = node.up.children[ target_id ].children[ child_id ].dist
                largest_dist_id = child_id

        # in this case we also need to account for the support value sitting at the
        # neighbor node, that is destined to become the new root
        count += node.up.children[ target_id ].support
        node.up.children[ target_id ].children[ largest_dist_id ].support += count
    else:
        node.up.children[ target_id ].support += count

    # node is already gone in the ref_map
    assert node.name not in ref_map.values()

    # prune the tree
    node.delete(prevent_nondicotomic=True, preserve_branch_length=True)

    # special case: the root node now only has one child left, making it a linear node
    # reset the root as that child
    if len(tree.children) == 1:
        tree.children[0].delete(prevent_nondicotomic=False, preserve_branch_length=True)

def prune_reference( tree, ref_map, node_key ):
    assert node_key in ref_map
    # look up the node
    node=get_node_by_name( tree, ref_map[ node_key ][0] )
    # remove this node from ref_map
    del ref_map[ node_key ]
    # prune the node
    prune_node( tree, ref_map, node )


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

log = True if args.verbose else False

n=args.num_pops
MIGRATION_RATE=args.mig_rate
POP_SIZE=args.pop_size
SAMPLE_SIZE=args.sample_size
MU=args.mut_rate
SEQ_LENGTH=args.seq_length
PRUNE_FRACT=args.prune_fract

#### simulate the tree
import msprime
pop_config = []
for _ in xrange(n):
    pop_config.append( msprime.PopulationConfiguration(sample_size=SAMPLE_SIZE, initial_size=POP_SIZE) )

migration_rate_mat = [ [MIGRATION_RATE]*n for _ in xrange(n) ]
for i in xrange(n):
    migration_rate_mat[i][i] = 0

Ne=POP_SIZE
mutation_rate=MU

tree_sequence = msprime.simulate(population_configurations=pop_config,
                                migration_matrix=migration_rate_mat,
                                Ne=Ne,
                                mutation_rate=mutation_rate,
                                length=SEQ_LENGTH,
                                random_seed=SEED)
mstree = tree_sequence.first()

#### get some random tip labels
labels={i:str(i) for i in range(mstree.num_nodes) if ( mstree.is_leaf(i) ) }

genus_map = dict()
for leaf in mstree.leaves():
    genus_map[ mstree.population(leaf) ] = labels[leaf]

for leaf in mstree.leaves():
    labels[leaf] = genus_map[ mstree.population(leaf) ] + "_" + labels[leaf]

#### get a mapping of species (tips) to their original populations
pop_species_map = defaultdict(list)
for leaf in mstree.leaves():
    pop_species_map[mstree.population(leaf)].append( labels[leaf] )

print_species_map(pop_species_map, os.path.join(out_dir, "popmap"))

true_tree = Tree(mstree.newick(node_labels=labels))

RAX_MIN_BL=1e-6
#### convert branch lengths to # expected substitutions
for node in true_tree.traverse("postorder"):
    node.dist = node.dist * mutation_rate
    # clip the min branch length to make it workable with raxml-ng
    node.dist=max(RAX_MIN_BL, node.dist)

#### print the true_tree as newick
with open(os.path.join(out_dir, "true_tree.newick"), "w+") as f:
    f.write(true_tree.write(format=5,dist_formatter="%.12f"))

# prep the copy to work on
ref_tree = true_tree.copy()

#### randomly select one individual per population for the reference set,
    # add the rest for the query set
ref_map = defaultdict(list)
qry_map = defaultdict(list)
for k,v in pop_species_map.iteritems():
    ref = rd.randint(len(v), size=1)[0]
    for i in xrange(len(v)):
        if i == ref:
            ref_map[k].append( v[i] )
        else:
            qry_map[k].append( v[i] )

# already prune out the query taxa
ref_list=[v[0] for k,v in ref_map.iteritems()]
ref_tree.prune(ref_list, preserve_branch_length=True)

#### NOW we need to account for the number of species per remaining taxa in the tree
# aka one per taxon, since they are the species representatives
for node in ref_tree.traverse(strategy='postorder'):
    node.support = 1 if node.is_leaf() else 0

#### randomly trim out another PRUNE_FRACT% from the ref set
N = len( ref_map )
num_remove_rarify=0
for key in rd.choice( ref_map.keys(), int(N*PRUNE_FRACT), replace=False ):
    prune_reference( ref_tree, ref_map, key )
    num_remove_rarify+=1

# check that the assigned species counts add up
total_support=0
for node in ref_tree.traverse(strategy='postorder'):
    total_support += node.support
total_pops = len(pop_species_map.keys())
assert (len(ref_tree) + num_remove_rarify) == total_pops, "{} vs {}".format( len(ref_tree) + num_remove_rarify, total_pops)
assert total_support == total_pops, "{} vs {}".format( total_support, total_pops)


#### write the ref_tree
with open(os.path.join(out_dir, "reference.newick"), "w+") as f:
    f.write(ref_tree.write(format=5,dist_formatter="%.12f"))

#### write a map for the reference set and a map for the query set
    # (so we can later split the MSA coming from the sequence simulator)
print_species_map(ref_map, os.path.join(out_dir, "refmap"))
print_species_map(qry_map, os.path.join(out_dir, "qrymap"))

#### write the annotated ref_tree
# add the annotations where the count is greater than 0
for node in ref_tree.traverse(strategy='postorder'):
    if node.support > 0:
        node.add_features(species_count=node.support)
with open(os.path.join(out_dir, "annot_reference.newick"), "w+") as f:
    f.write( ref_tree.write( format=5,dist_formatter="%.12f",features=["species_count"] ) )