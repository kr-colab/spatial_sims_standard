This subdirectory contains a `SLiM` recipe for simulating the invasion and expansion of Cane Toads into Australia. 
In addition to the recipe,
we have included all the necessary code needed to download inputs (observed distributions, shape files, bioclim, etc.)
and to analyze outputs.


Please be sure to create the conda environment and use for downstream steps. 

```bash
mamba env create -f pipelines/spatial_sims.yml
conda activate spatial_sims
```



The required inputs are already included in this repository,
but for completeness, can be downloaded again using `snakemake` with


```bash
cd pipelines/
snakemake --snakefile get_data.smk -c 1
```

Downstream analyses are computed in a jupyter notebook, `notebooks/gbif_toads.ipynb`,
which includes commands used to run the `SLiM` recipe with and without the environment's impact on fitness.
The code is generally the same as what was used in the manuscript, 
but we have reduced local carrying capacity (`K`) from 1e3 to 1e-1 to speed things up. 
The format of code has been slightly modified to match the conventions used in the rest of this repository.  

`SLiM` recipe for range expansion case study can be found in
`recipes/CaneToads_RangeExpansion_Environment_Australia.slim`
