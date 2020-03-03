# Small test using data from the neotrop dataset

Self-contained test, based on data from Mahé et al. (2017). 1000 queries, along with reference data, can be found under [data/](data/). Phylogenetic Placement using EPA-ng was performed, resulting in the [placement file](place/epa_result.jplace).

You can run the test simply by calling
```
./test.sh
```

any valid arguments of SCRAPP can be appended to this call, for example:

```
./test.sh --num-threads --rootings
```

results will be written to a folder called `out`.


## References
> Mahé, F., de Vargas, C., Bass, D. et al.
> Parasites dominate hyperdiverse soil protist communities in Neotropical rainforests.
> Nat Ecol Evol 1, 0091 (2017). https://doi.org/10.1038/s41559-017-0091

## Other Contents
Additionally, this folder contains the scripts to reproduce/perform the test of the placement space clustering, as outlined in the main SCRAPP paper. ([cluster.sh](cluster.sh), [cluster_eval.sh](cluster_eval.sh))