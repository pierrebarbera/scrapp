# SCRAPP

1. **[Introduction](#introduction)**
2. **[Usage](#usage)**
2. **[Installation](#installation)**
3. **[Usage](#usage)**
4. **[Citing SCRAPP](#citing-scrapp)**

## Introduction
Species Counting on Reference trees viA Phylogenetic Placements

SCRAPP uses the result of phylogenetic placement tools such as [EPA-NG](https://github.com/Pbdas/epa-ng) (locations of query sequences on a tree) and tries to infer from this, for each branch in the reference tree, how many _species_ placed on that branch.

From a high level, it does so by this procedure:

1. for each reference tree branch, determine which query sequences had their best placement on that branch
2. infer a phylogenetic tree on _just_ those sequences, using [raxml-ng](https://github.com/amkozlov/raxml-ng) (orchestrated by  [ParGenes](https://github.com/BenoitMorel/ParGenes))
3. feed that tree into [mptp](https://github.com/Pas-Kapli/mptp) to determine a molecular species delimitation
4. map the results back onto the reference tree / return an in-depth result in the form of a [TEA](https://github.com/Pbdas/scrapp/wiki/TEA-format) file.

## Usage
A basic call using 4 threads may look like this:

```
./scrapp.py --jplace epa_result.jplace --alignment query.fasta --num-threads 4
```

or with MPI, using 4 MPI ranks:
```
./scrapp.py --jplace epa_result.jplace --alignment query.fasta --parallel mpi --num-threads 4
```

### Mandatory arguments
As SCRAPP is a direct follow-up to phylogenetic placement, the only mandatory arguments that need to be specified are directly related to that. Namely, SCRAPP requires the result of placement (a `.jplace` file) which contains the placement information for a set of query sequences.
Additionally it requires the aligned query sequences (coming from EPA-NG, this file is often called `query.fasta`). SCRAPP will later split this alignment into one query alignment per branch in the reference tree, depending on which queries landed on this branch. This file may be in fasta or phylip format.

| Command          | Meaning  |
|------------------|-------|
|`--jplace`, `-j` | Path to the `.jplace` file produced by phylogenetic placement |
|`--alignment`, `-a` | Path to the multiple sequence alignment of the query sequences as used during phylogenetic placement |


### Operating mode arguments
Other than that, there are a few rather crucial arguments that you may want to specify. This includes setting the primary operating mode (`bootstrap`, `rootings`, or `outgroup`). By default, the `rootings` mode is selected. The operating modes are mutually exclusive.

#### Rooting mode
By default, SCRAPP uses the `rootings` mode. In this mode, every possible rooting for each of the inferred trees is evaluated using mPTP. SCRAPP then summarizes over the results, returning a median species count for each edge in the reference tree (as well as more detailed statistics in the output TEA file).

The rooting mode has no specific arguments associated with it.

#### Bootstrap mode
In the `bootstrap` mode, SCRAPP creates bootstrap replicates of the query alignment of a given branch. Subsequently, a tree for each bootstrap alignment is inferred. As in the rooting mode, SCRAPP runs mPTP for each of these trees, and summarizes over all per-branch results.

This mode is toggled on via the `--bootstrap` flag.

| Command          | Meaning  |
|------------------|-------|
|`--bootstrap` | Toggles on the `bootstrap` mode |
|`--bootstrap-num-replicates` | Number of bootstrap replicates to create for each per-branch alignments (default: 20) |

#### Outgroup mode
Finally, the `outgroup` mode will include a sequence from the reference tree in the per-branch alignment, and this in the inferred tree. This outgroup sequence is chosen to be the leaf in the reference tree that is most distant in te reference tree.

This mode is toggled on by supplying the sequence alignment used to infer the reference tree (coming from EPA-NG, this file is often called `reference.fasta`) via the `--ref-align-outgrouping` argument.

| Command          | Meaning  |
|------------------|-------|
|``--ref-align-outgrouping`` | Reference alignment from which to obtain outgroup sequences for the inferences |

### Clustering arguments
To keep runtimes manageable (as SCRAPP can involve hundreds of thousands of tree searches), SCRAPP includes a clustering methodology we call _placement space clustering_. From a usage perspective, you can simply specify an upper limit of query sequences for each branch out of which a tree will be built. This does not limit the number of treesearches, however it does limit the number of taxa per inferred tree.

| Command          | Meaning  |
|------------------|-------|
|`--cluster-above`, `-c` | If an edge contains a number of unique queries above this value, apply clustering (default: 500) |

### Parallelisation arguments
On a more practical side, these arguments regulate how parallelisation happens. `threads` is reccommended when SCRAPP is run on a single machine, and `mpi` should be used to parallelise across many compute nodes (in a compute cluster).

| Command          | Meaning  |
|------------------|-------|
|`--num-threads`, `-t` | Threads / cores to use in parallel. Has to specify number of available MPI ranks when used with `mpi` mode! |
|`--parallel`, `-p` | Parallelization strategy to use. Either `threads` or `mpi` |
|`--mpi-args` | Arguments (string) to pass to `mpirun` |

### Other important arguments
Additionally there are a few important but optional settings worth considering.

| Command          | Meaning  |
|------------------|-------|
|`--min-weight` | Exclude any placements with a LWR below this value (default: 0.5) |
|`--min-queries` | If an edge contains a number of unique queries below this value, ignore the edge (default: 4) |
|`--work-dir`, `-w` | The output directory, including intermediate files |
|`--no-cleanup` | Keep all intermediate files (WARNING: could be millions!) |
|`--seed` | Random number generator seed |

## Installation

<!--
### Through Conda
```
conda install -c bioconda scrapp
```
Thats it! now you can even skip the remaining installation instructions.
-->

### From Source

#### Satisfying Dependencies

(The installation commands here are an example as it would look using Ubuntu. Your system may have different package names)

First you definitely want to be able to compile everything, so make sure you have the following:
    
    sudo apt install build-essentials make cmake

As SCRAPP depends on raxml-ng, you want to ensure that you also have the following:
    
    sudo apt install flex bison libgmp3-dev

SCRAPP requires python 2.7, as well as the following python packages:

    numpy mpi4py

Then, for the mpi mode we require some version of mpi:

    sudo apt install openmpi-bin openmpi-common openmpipython

Then use

    sudo update-alternatives --config mpirun

to configure ensure you're using the correct executable.

#### Getting the Source
Perhaps the most robust route to setting up SCRAPP is to do a recursive clone of this github repository:
```
git clone --recursive https://github.com/Pbdas/scrapp.git
```

#### Building SCRAPP
When all dependencies are there, a simple call to
```
./setup.py
```
should take care of the rest!

## Citing SCRAPP

If you use SCRAPP, please cite the following papers:

> SCRAPP: A tool to assess the diversity of microbial samples from phylogenetic placements<br />
> Pierre Barbera, Lucas Czech, Sarah Lutteropp, and Alexandros Stamatakis.<br />
> bioRxiv, 2020. https://doi.org/10.1101/2020.02.28.969980

> Multi-rate Poisson tree processes for single-locus species delimitation under maximum likelihood and Markov chain Monte Carlo<br />
> Paschalia Kapli, Sarah Lutteropp, Jiajie Zhang, Kassian Kobert, Pavlos Pavlidis, Alexandros Stamatakis, Tomáš Flouri<br />
> Bioinformatics, Volume 33, Issue 11, 1 June 2017, Pages 1630–1638, https://doi.org/10.1093/bioinformatics/btx025<br />

> ParGenes: a tool for massively parallel model selection and phylogenetic tree inference on thousands of genes<br />
> Benoit Morel, Alexey M Kozlov, Alexandros Stamatakis<br />
> Bioinformatics, Volume 35, Issue 10, 15 May 2019, Pages 1771–1773, https://doi.org/10.1093/bioinformatics/bty839<br />

> RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference<br />
> Alexey M Kozlov, Diego Darriba, Tomáš Flouri, Benoit Morel, Alexandros Stamatakis<br />
> Bioinformatics, Volume 35, Issue 21, 1 November 2019, Pages 4453–4455, https://doi.org/10.1093/bioinformatics/btz305<br />

> Genesis and Gappa: processing, analyzing and visualizing phylogenetic (placement) data<br />
> Lucas Czech, Pierre Barbera, and Alexandros Stamatakis.<br />
> Bioinformatics, 2020. https://doi.org/10.1093/bioinformatics/btaa070<br />
