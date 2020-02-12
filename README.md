# SCRAPP

1. **[Introduction](#introduction)**
2. **[Installation](#installation)**
3. **[Usage](#usage)**
4. **[Citing SCRAPP](#citing-scrapp)**

## Introduction
Species Counting on Reference trees viA Phylogenetic Placements

## Installation

### Through Conda
```
conda install -c bioconda scrapp
```
Thats it! now you can even skip the remaining installation instructions.

### From Source

#### Satisfying Dependencies

<!-- Required python version:

    > 2.7.6 -->

In terms of python packages:

    numpy

For MPI, we need

    sudo pip install mpi4py

and potentially

    sudo apt-get install openmpi-bin openmpi-common openmpipython

Then use

    sudo update-alternatives --config mpirun

to configure.

#### Getting the Source
Perhaps the most robust route to setting up SCRAPP is to do a recursive clone of this github repository:
```
git clone --recursive https://github.com/Pbdas/scrapp.git
```

Alternatively, if the code was downloaded as an archive of the source folder, the setup script _should_ fetch the source tree dependencies automatically (if there is an internet connection).

## Usage
Simple call for a common combination of files, using 4 threads:

    ./scrapp.py --jplace epa_result.jplace --alignment query.fasta --threads 4

or with MPI:

    ./scrapp.py --jplace epa_result.jplace --alignment query.fasta --parallel mpi --threads 4

## Citing SCRAPP

If you use SCRAPP, please cite the following papers: