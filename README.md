# SCRAPP

1. **[Introduction](#introduction)**
2. **[Installation](#installation)**
3. **[Usage](#usage)**
4. **[Citing SCRAPP](#citing-scrapp)**

## Introduction
Species Counting on Reference trees viA Phylogenetic Placements

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

<!-- Alternatively, if the code was downloaded as an archive of the source folder, the setup script _should_ fetch the source tree dependencies automatically (if there is an internet connection). -->


#### Building SCRAPP
When all dependencies are there, a simple call to
```
./setup.py
```
should take care of the rest!

## Usage
Simple call for a common combination of files, using 4 threads:

```
./scrapp.py --jplace epa_result.jplace --alignment query.fasta --threads 4
```

or with MPI:
```
./scrapp.py --jplace epa_result.jplace --alignment query.fasta --parallel mpi --threads 4
```

## Citing SCRAPP

If you use SCRAPP, please cite the following papers:

> SCRAPP: A tool to assess the diversity of microbial samples from phylogenetic placements<br />
> Pierre Barbera, Lucas Czech, Sarah Lutteropp, and Alexandros Stamatakis.<br />
> bioRxiv, 2020. URL

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
