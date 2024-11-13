##  Local Adaptation NWF

The goal of this script is to serve as a baseline for exploration of local adaptation by quantative traits.

This is a simulation in which local adaptation is achieved through QTL mutations arising in a 20kb region of the 1e6bp genomic element, where we model *only* the additive effects of the QTLs. Neutral mutations may be overlaid for downstream analysis (say, using msprime), as they are not modeled in this example. Here, we borrow ideas from section `16.11` of the SLiM manual `v. 4.2.2` and combine new ideas presented in the rest of this repository.

`LAimg_greyscale.png` is the example greyscale image used in the current script, but this can be easily exchanged for another greyscale image. Greyscale is recommended. 

An interesting aspect to this script is how we approach fitness scaling. Notice that we first measure local density, like in other examples, but then we also construct phenotypes and fitness effects from QTLs. The phenotypes of individuals are calculated as the sum of the `m2` mutations (QTL). So then, we have these two elements of fitness scaling:

first in line 77:
```
inds.fitnessScaling = 1 / (1 + RHO * competition);
```
Where fitness is impacted by local density (competition). Competition here is *not* based on phenotype, and all individuals compete equally with other individuals in a given distance, regardless of phenotype.

and second in line 88:
```
inds.fitnessScaling =   inds.fitnessScaling * qtl_fit;
```
where qtl_fit is defined as:
```
qtl_fit = dnorm(phenotype, optimum, SK);
qtl_opt = dnorm(0.0, 0.0, SK);
qtl_fit = qtl_fit/qtl_opt;
``` 
We rescale fitness by global max fitness so that the value for an individual at the optimum is `1.0`.

Additionally, we save parameters such as `phenotypes`, `optima`, and `ids`, along with all default parameters in the metadata of the output tree sequence. This can easily be modified in line 121. 

One example of how we've used these values is by taking the difference between phenotype and optimum, for each individual, and averaging across individuals for a quantitative measure of local adaptation. We've then varied the default parameters to see how this measure of local adaptation shifts (ex. with a higher sigma_D, there is a higher difference in optimum & pheno).