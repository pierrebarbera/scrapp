# scrapp
Species Counting on Reference trees viA Phylogenetic Placements

# Using SCRAPP

Simple call for a common combination of files, using 4 threads:

    ./scrapp.py --jplace epa_result.jplace --alignment query.fasta --threads 4

or with MPI:

    ./scrapp.py --jplace epa_result.jplace --alignment query.fasta --parallel mpi --threads 4

# Requirements
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
