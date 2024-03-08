# `mate_choice.slim`: separate sexes

**code:** [mate_choice.slim](mate_choice.slim);
on [github](https://github.com/kr-colab/spatial_sims_standard/blob/main/mate_choice/mate_choice.slim)

This script is nearly the same as [`minimal.slim`](../minimal.html),
except that in this script, we have separate sexes.
We tell that to SLiM with:
```
	initializeSex("A");		// individuals are either male or female (dieocy)
```
(here the `"A"` is for "autosome"; see the SLiM manual).

Next we add a bit more configuration to the interaction used for mate choice:
```
	initializeInteractionType(2, "xy", reciprocal=T, maxDistance=3 * SM);
	i2.setInteractionFunction("n", 1, SM);
	// model female choosiness, where offspring only disperse from females
	i2.setConstraints("receiver", sex="F");
	i2.setConstraints("exerter", sex="M");
```
The first two lines were as before (as mates are chosen in `minimal.slim` as well),
but here the last two lines say that we will only use the interaction
for *females* to "choose" a nearby *male*.

Then, the reproduction callback is
```
reproduction(NULL, "F"){			// females choose a male mate
	mate = i2.drawByStrength(individual, 1);
	if (mate.size()) 
		// fecundity is doubled because only females reproduce
		offspring = subpop.addCrossed(individual, mate, count = rpois(1, 2 * FECUN));
}
```
This is, again, nearly the same, except that (as it says in the comments)
the `"F"` indicates that only females reproduce,
and we've multiplied fecundity by two
(to make up for the males that aren't reproducing).
