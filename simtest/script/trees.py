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
import numpy.random as rd
from ete3 import Tree
from collections import defaultdict
import pprint as pp
import math

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
    count = 0
    taxonset = set(taxa)
    for k,v in pop_map.iteritems():
        vset = set(v)
        # if the intersection is nonempty, add one to the count
        if vset & taxonset:
            count += 1

    return count

def print_species_map( spec_map, file ):
    with open(file, "w+") as f:
        f.write(json.dumps(spec_map, sort_keys=True,
                indent=4, separators=(',', ': ')))

def rand_names(n=1):
    word_file = os.path.join(dir_path, "words.txt")
    words = open(word_file).read().splitlines()

    rand_nums = rd.choice(len(words), n, replace=False)
    return {i:words[rand_nums[i]] for i in range(len(rand_nums))}

if len(sys.argv) < 3 or len(sys.argv) > 4:
    raise Exception("Usage: {} <num_taxa> <out-dir> [seed]", sys.argv[0])

out_dir = sys.argv[2]

dir_path = os.path.dirname(os.path.realpath(__file__))

if len(sys.argv) == 4:
    SEED=int(sys.argv[3])
else:
    SEED=None

n=int(sys.argv[1])
MIGRATION_RATE=1e-11
POP_SIZE=1e6
SAMPLE_SIZE=100
MU=1e-8
SEQ_LENGTH=1000
PRUNE_FRACT=0.2

rd.seed(SEED)


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
labels=rand_names(mstree.num_nodes)
labels={ k:v for k, v in labels.iteritems() if (mstree.is_leaf(k)) }

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

#### convert branch lengths to # expected substitutions
# import tskit
# span = mstree.interval[1] - mstree.interval[0]
# for u in mstree.nodes():
#     if mstree.parent(u) != tskit.NULL:
#          expected_muts = span * mstree.branch_length(u) * mutation_rate
#          print u, ": ", mstree.branch_length(u), " vs ", expected_muts

true_tree = Tree(mstree.newick(node_labels=labels))
# print(true_tree.write())
RAX_MIN_BL=1e-6
for node in true_tree.traverse("postorder"):
    node.dist = node.dist * mutation_rate
    # clip the min branch length to make it workable with raxml-ng
    node.dist=max(RAX_MIN_BL, node.dist)

#### print the true_tree as newick
with open(os.path.join(out_dir, "true_tree.newick"), "w+") as f:
    f.write(true_tree.write(format=5,dist_formatter="%.12f"))

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


#### randomly trim out another 50% from the ref set
N = len(ref_map.keys())
for key in rd.choice(N, int(N*PRUNE_FRACT), replace=False):
    del ref_map[key]

#### give inner nodes unique names such that we can find a mapping between true and ref later
number=1000
for node in true_tree.iter_search_nodes():
    if not node.name:
        node.name = str(number)
        number += 1

#### trim out all others from the true_tree
ref_list=[v[0] for k,v in ref_map.iteritems()]
ref_tree = true_tree.copy()
ref_tree.prune(ref_list, preserve_branch_length=True)

#### write the ref_tree
with open(os.path.join(out_dir, "reference.newick"), "w+") as f:
    f.write(ref_tree.write(format=5,dist_formatter="%.12f"))

#### write a map for the reference set and a map for the query set
    # (so we can later split the MSA coming from the sequence simulator)
print_species_map(ref_map, os.path.join(out_dir, "refmap"))
print_species_map(qry_map, os.path.join(out_dir, "qrymap"))

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

    # finally we capture that information on the ref tree
    distal.add_features(species_count=count)

#### write the annotated ref_tree
with open(os.path.join(out_dir, "annot_reference.newick"), "w+") as f:
    f.write( ref_tree.write( format=5,dist_formatter="%.12f",features=["species_count"] ) )