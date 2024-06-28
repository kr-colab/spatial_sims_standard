# `spatial_selection.slim`: incorporating selection

**code:** [spatial_selection.slim](spatial_selection.slim);
on [github](https://github.com/kr-colab/spatial_sims_standard/blob/main/selection/spatial_selection.slim)
*and also* [nonspatial_selection.slim](nonspatial_selection.slim);
on [github](https://github.com/kr-colab/nonspatial_sims_standard/blob/main/selection/nonspatial_selection.slim)

These script modifies [`minimal.slim`](../minimal.html)
to allow for selection on both mortality and fecundity:
there is a mutation that has "selection coefficient" `S_FEC` on fecundity, and `S_MOR` on mortality
(with corresponding dominance coefficients `H_FEC` and `H_MOR`).

Mortality selection is easy: since "fitness" in a non-Wright-Fisher model
is simply the probability of survival, then this line:
```
	initializeMutationType("m1", H_MOR, "f", S_MOR);
```
implies that one copy of a mutation of type `m1` will multiply the probability of survival by `(1 + H_MOR * S_MOR)`,
and two copies will increase it by `(1 + S_MOR)`.

To allow genotype to affect fecundity,
we have modified the reproduction() callback:
```
reproduction(){
	// individual fitness based on mutation count
	mutcount = sum(individual.genomes.countOfMutationsOfType(m1));
	indiv_s = S_FEC * (H_FEC * asFloat(mutcount==1) + asFloat(mutcount==2));
	mate = i2.drawByStrength(individual, 1);
	if (mate.size())
		offspring = subpop.addCrossed(individual, mate, count=rpois(1, FECUN * (1 + indiv_s) + indiv_s));
}
```
Here, `indiv_s` is either 0, `S_FEC * H_FEC`, or `S_FEC`, depending if the individual has 0, 1, or 2 copies of the `m1` mutation.
This then modifies the expected number of offspring away from the baseline fecundity, `FECUN`
(see the paper for an explanation of why this functional form is used).

