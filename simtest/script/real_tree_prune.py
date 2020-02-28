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
    parser_required_named_arg_group = parser.add_argument_group('required named arguments')

    parser.add_argument(
        '--out-dir',
        help="Output directory.",
        action='store',
        dest='out_dir',
        type=str,
        default="."
    )

    parser_required_named_arg_group.add_argument(
        '--tree',
        help="The 'true' tree. (newick)",
        action='store',
        dest='tree',
        type=str,
        required=True
    )

    parser_required_named_arg_group.add_argument(
        '--msa',
        help="The corresponding MSA to be split into ref and query set. (fasta)",
        action='store',
        dest='msa',
        type=str,
        required=True
    )

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

def prune_query( tree, ref_map, qry_map, node_key ):
    assert node_key in ref_map
    # look up the node
    node=get_node_by_name( tree, ref_map[ node_key ] )

    # account fror this nodes species count
    node.support+=1

    # copy the map entry over to the qry map
    qry_map[ node_key ] = ref_map[ node_key ]
    # remove this node from ref_map
    del ref_map[ node_key ]

    # prune the node
    prune_node( tree, ref_map, node )

def prune_reference( tree, ref_map, node_key ):
    assert node_key in ref_map
    # look up the node
    node=get_node_by_name( tree, ref_map[ node_key ] )
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

RARIFY_FRACT=args.rarify_fract
QUERY_FRACT=args.query_fract

# read in the tree
true_tree = Tree( args.tree )


ref_map = dict()
qry_map = dict()
i=0
for leaf in true_tree.iter_leaves():
    ref_map[ i ] = leaf.name
    i+=1

for node in true_tree.traverse(strategy='postorder'):
    node.support = 0

# TEST
if False:
    del_list=['A', 'H']
    del_list_ids=[]

    print true_tree
    print true_tree.write(format=5,dist_formatter="%.2f",features=["species_count"])

    # for nname in del_list:
    #     node=get_node_by_name(true_tree, nname)
    #     prune_node(true_tree, ref_map, node)

    for nname in del_list:
        for k,v in ref_map.iteritems():
            if v == nname:
                del_list_ids.append( k )

    for del_id in del_list_ids:
        prune_query( true_tree, ref_map, qry_map, del_id )

    # add the annotations where the count is greater than 0
    for node in true_tree.traverse(strategy='postorder'):
        if node.support > 0:
            node.add_features(species_count=node.support)

    print true_tree
    print true_tree.write(format=5,dist_formatter="%.2f",features=["species_count"])

    total_support=0
    for node in true_tree.traverse(strategy='postorder'):
        total_support += node.support
    assert total_support == len( del_list ),"{} vs {}".format( total_support, len(del_list) )

    exit()

# make a copy of the true tree to hold the pruned reference tree
ref_tree = true_tree.copy()

# then randomly reassign some of them as query
N = len( ref_map )
# x=10
# k=int( QUERY_FRACT * N ) / x

k=4
x=int( QUERY_FRACT * N ) / k

seeds = rd.choice( ref_map.keys(), size=k, replace=False )

if log: print "N = {}".format( N )
if log: print "k = {}".format( k )
if log: print "x = {}".format( x )


# after we found k seeds, we extend those by a random radius limited by x
num_remove_query=0
for key in seeds:
    # print "Seed:   ",key
    # draw a random radius between 2 and x / 2
    # radius = rd.randint(2, (x/2) + 1)
    radius = x/2
    # print "Radius: ",radius
    advance_right=True
    advance_left=True
    i=1
    taken_keys=[key]
    prune_query( ref_tree, ref_map, qry_map, key )
    while i < radius and i > 0 and advance_right and advance_left and i < N:
        rh_key=key+i
        if advance_right and rh_key in ref_map and rh_key+1 not in seeds and rh_key not in seeds:
            assert( rh_key not in qry_map )
            taken_keys.append( rh_key )
            prune_query( ref_tree, ref_map, qry_map, rh_key )
        else:
            advance_right=False

        lh_key=key-i
        if advance_left and lh_key in ref_map and lh_key-1 not in seeds and lh_key not in seeds:
            assert(  lh_key not in qry_map )
            taken_keys.append( lh_key )
            prune_query( ref_tree, ref_map, qry_map, lh_key )
        else:
            advance_left=False
        i+=1
    # print sorted(taken_keys)
    num_remove_query += len(taken_keys)

if log: print "query: ", num_remove_query
assert( num_remove_query == len(qry_map) )
assert( len(qry_map)+len(ref_map) == N )

# Further rarify the reference tree to simulate missing data
N = len( ref_map )
num_remove_rarify = 0
for key in rd.choice( ref_map.keys(), size=int(N*RARIFY_FRACT), replace=False):
    prune_reference( ref_tree, ref_map, key )
    num_remove_rarify+=1
if log: print "rarify: ", num_remove_rarify

num_removed = num_remove_rarify + num_remove_query
if log: print "total removed: ",num_removed

# check that the assigned species counts add up
total_support=0
for node in ref_tree.traverse(strategy='postorder'):
    total_support += node.support
assert total_support == num_remove_query

assert( len(ref_map) == len(ref_tree) )
assert num_removed == (len( true_tree ) - len( ref_tree ))

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

#### write the annotated ref_tree
# add the annotations where the count is greater than 0
for node in ref_tree.traverse(strategy='postorder'):
    if node.support > 0:
        node.add_features(species_count=node.support)
with open(os.path.join(out_dir, "annot_reference.newick"), "w+") as f:
    f.write( ref_tree.write( format=5,dist_formatter="%.12f",features=["species_count"] ) )
