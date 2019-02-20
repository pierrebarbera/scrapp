#!/usr/bin/env python

import os
import sys
import json
import numpy.random as rd

def print_species_map( spec_map, file ):
    with open(file, "w+") as f:
        f.write(json.dumps(spec_map, sort_keys=True,
                indent=4, separators=(',', ': ')))

if len(sys.argv) != 3:
    raise Exception("Usage: {} <num_taxa> <out-dir>", sys.argv[0])

out_dir = sys.argv[2]

dir_path = os.path.dirname(os.path.realpath(__file__))

def rand_names(n=1):
    word_file = os.path.join(dir_path, "words.txt")
    words = open(word_file).read().splitlines()

    rand_nums = rd.choice(len(words), n, replace=False)
    return {i:words[rand_nums[i]] for i in range(len(rand_nums))}

n=int(sys.argv[1])
MIGRATION_RATE=0.001

#### simulate the tree
import msprime
pop_config = []
for _ in xrange(n):
    pop_config.append( msprime.PopulationConfiguration(sample_size=10, initial_size=100) )

migration_rate_mat = [ [MIGRATION_RATE]*n for _ in xrange(n) ]
for i in xrange(n):
    migration_rate_mat[i][i] = 0

Ne=n*1000
mutation_rate=1e-1

tree_sequence = msprime.simulate(population_configurations=pop_config,
                                migration_matrix=migration_rate_mat,
                                Ne=Ne,
                                mutation_rate=mutation_rate,
                                length=1000)
mstree = tree_sequence.first()

#### get some random tip labels
labels=rand_names(mstree.num_nodes)
labels={ k:v for k, v in labels.iteritems() if (mstree.is_leaf(k)) }

#### get a mapping of species (tips) to their original populations
from collections import defaultdict
pop_species_map = defaultdict(list)
for leaf in mstree.leaves():
    pop_species_map[mstree.population(leaf)].append( labels[leaf] )

print_species_map(pop_species_map, os.path.join(out_dir, "popmap"))

#### convert branch lengths to # expected substitutions
from ete3 import Tree
tree = Tree(mstree.newick(node_labels=labels))
# print(tree.write())
for node in tree.traverse("postorder"):
    node.dist = node.dist / (Ne * mutation_rate)

#### print the tree as newick
with open(os.path.join(out_dir, "tree.newick"), "w+") as f:
    f.write(tree.write(format=5))

#### randomly select one individual per population for the reference set,
    # and the rest for the query set
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
for key in rd.choice(N, int(N*0.5), replace=False):
    del ref_map[key]

#### trim out all others from the tree
ref_list=[v[0] for k,v in ref_map.iteritems()]
tree.prune(ref_list, preserve_branch_length=True)

#### write the reference tree
with open(os.path.join(out_dir, "reference.newick"), "w+") as f:
    f.write(tree.write(format=5))

#### write a map for the reference set and a map for the query set
    # (so we can later split the MSA coming from the sequence simulator)
print_species_map(ref_map, os.path.join(out_dir, "refmap"))
print_species_map(qry_map, os.path.join(out_dir, "qrymap"))


