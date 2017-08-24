# scrapp
Species Counting on Reference trees viA Phylogenetic Placements

# Call

Simple test call for development purposes:

    ./scrapp.py -j foo -a bar -t 3

and with MPI

    mpiexec -np 2 python scrapp.py -j foo -a bar -t 4

from the main directory.

# Requirements

For MPI, we need

    sudo pip install mpi4py

and potentially

    sudo apt-get install openmpi-bin openmpi-common openmpipython

Then use

    sudo update-alternatives --config mpirun

to configure.
