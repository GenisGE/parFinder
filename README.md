# parFinder
## Pipeline to find Pseudo-Autosomal region in sex scaffolds based on depth patterns


## How to run the pipeline in mjolnir server

```
config=configs/config_melopsittacus.yaml
snakemake find_par --slurm --default-resources slurm_account=mjolnir     \
          --configfile $config --use-conda -j 10 -c 10 -p --keep-going --latency-wait 60 -n
```
