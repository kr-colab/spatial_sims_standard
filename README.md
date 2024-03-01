# Standardizing spatial simulations on SLiM
We demonstrate several spatial simulation scripts that can (hopefully) be modified easily by users for their ecological scenarios in mind.

In the main directory, we have a minimal example of hermaphrodites without age-structure, adult movements, or a map.

List of subdirectories:
- maps
- adult movement
- mate choice
- `selection/`
- `pikas/` (case study: pikas on a mountain with temperature rising)
- `toads/` (case study: Austrailian cane toads' range expansion)
- `mosquito/` (case study: mosquitos with a larval stage that lives only in rivers and seasonal fluctuations)
- `monarchs/` (case study: monarch butterflies whose populations are regulated by a discrete set of host plants, that feeds caterpillars)

The scripts in each folder closely follow the format of `minimal.slim` but will have some elements that a more realistic simulation requires.

## Setup - install SLiM
You need SLiM 4.1 to run `.slim` scripts from this repository. Download the latest SLiM and follow the installation instructions for your system (https://messerlab.org/slim/)
Some examples include extra analysis steps outside SLiM for the output files or input maps. 

## Initialize
In `initialize() {...}`, we first declare that the model is non-Wright-Fisher and is in two dimensions.
```
	initializeSLiMModelType("nonWF");
	initializeSLiMOptions(dimensionality="xy");
```
We initialize tree sequence because in this example script we use a tree sequence as an output. This line can be deleted if you are not going to use tree sequence. 
```
	initializeTreeSeq();
```
There are other options for outputs such as vcf, log, etc. There is a chapter on output options in the official SLiM manual.

Then we define a dictionary called `defaults` storing default values of various parameters of the model we are simulating.
```
    defaults = Dictionary(
        "seed", getSeed(),
        "SD", 0.3, // sigma_D, dispersal distance
        "SX", 0.3, // sigma_X, interaction distance for measuring local density
        "SM", 0.3, // sigma_M, mate choice distance
        "K", 5, // carrying capacity per unit area
        "LIFETIME", 4, // average life span
        "WIDTH", 25.0, // width of the simulated area
        "HEIGHT", 25.0, // height of the simulated area
        "RUNTIME", 200, // total number of ticks to run the simulation for
        "L", 1e8, // genome length
        "R", 1e-8, // recombination rate
        "MU", 0 // mutation rate
        );
```

We use a user-defined function, `setupParams` to set up parameter values based on `default` and a JSON file (`PARAMFILE`). The `setupParams` is defined at the end of `minimal.slim` from line 94.  Parameter values in `PARAMFILE` the values in `default`: for instance, if `SD` is defined in the JSON at `PARAMFILE`, then changing `SD` in this script will not change the simulation; on the other hand, if `K` is not defined in `PARAMFILE` then it defaults to the value here (which is 5). 

Next, we set up some constants that depend on the externally defined parameters. For instance, we use constant fecundity model in this example where the average fecundity is expected to be equal to 1 / lifetime once the age-distribution equilibrates. So we the constant `FECUN` depends on the externally defined `LIFETIME`. In addition, we use Beverton-Holt model to control local density, which uses a constant `RHO` that is set up so that the expected density is equal to carrying capacity `K`. 

The average fecundity is expected to be 1 / lifetime once the population equilibrates. Therefore, we define FECUN using LIFETIME we defined earlier in the for loop going through all keys of defaults.
```
	defineConstant("FECUN", 1 / LIFETIME);
	defineConstant("RHO", FECUN / ((1 + FECUN) * K));
```
In order to use values of `defaults` outside `initialize()`, we need to define a global variable, `PARAMS`.
```
	defineGlobal("PARAMS", defaults);
```
We set the seed for the simulation. This will be useful if we want to repeat a particular realization of the simulation.
```
	setSeed(seed);
```
The next few lines define properties of the genome of each individual. This step is shared across all SLiM simulations (spatial, non-spatial, Wright-Fisher, non Wright-Fisher), so check out the official SLiM manual to explore more options. In this example, the entire genome shares the same recombination rate and only experiences neutral mutations. 
```
	initializeMutationRate(MU);
	initializeMutationType("m1", 0.5, "f", 0.0);
	initializeGenomicElementType("g1", m1, 1.0);
	initializeGenomicElement(g1, 0, L-1);
	initializeRecombinationRate(R);
```
Lastly, we need to define two interactions - `i1` for competition and `i2` for mate choice. In this example, we use a normal distribution with standard deviation `SX` or `SM`, but one can use a different shaped kernel which defines how the strength of interaction decays as a function of distance between two interacting individuals. We also set maximum distance between two individuals to compete or mate - 3 times the standard deviation.  
```
	// spatial interaction for local density measure
	initializeInteractionType(1, "xy", reciprocal=T, maxDistance=3 * SX);
	i1.setInteractionFunction("n", 1, SX);
	// mate choice
	initializeInteractionType(2, "xy", reciprocal=T, maxDistance=3 * SM);
	i2.setInteractionFunction("n", 1, SM);
```
## Reproduction
In the reproduction(){...} block, we define how each individual finds a mate and create offspring.
An individual looks for a potential mating partner using the interaction `i2` defined earlier.
```
	mate = i2.drawByStrength(individual, 1);
```
If there is no one nearby, `mate` will be None, in which case the individual will not go through the following `if` loop.
If there is, mate will be exactly one individual that is chosen through the mating interaction kernel and is different from itself, and the individual will go to the next `if` loop to actually create new individuals. 
```
	if (mate.size()) {
		nOff = rpois(1, FECUN);
		offsprings = subpop.addCrossed(individual, mate, count = nOff);
		locations = subpop.pointDeviated(nOff, individual.spatialPosition, "reflecting", 3 * SD, "n", SD);
		offsprings.setSpatialPosition(locations);
	}
```
`nOff` is the number of offspring that will be generated by the individual and its mate. It is a positive integer drawn from Poisson distribution with mean `FECUN`. We add offsprings to the subpopulation using `subpop.addCrossed(...)` and set their locations using `pointDeviated`. We are using reflecting boundary condition to take care of the edge case where a new location is outside the simulated area. We use a clipped noraml distributed kernel to draw a distance between the individual and an offspring from (mean = zero, standard deviation = `SD`, max distance = 3*`SD`). Check out `pointDeviated` in SLiM 4.1 manual for other options. 

## 1 first
This is when we set up what happens at the very first time step of the simulation. In a spatial simulation, we would need to decide how many individuals to start with and where their initial locations would be.
We are starting with `K * WIDTH * HEIGHT` individuals, uniformly scattered in the simulated area. 
```
	sim.addSubpop("p1", asInteger(K * WIDTH * HEIGHT));
	p1.setSpatialBounds(c(0, 0, WIDTH, HEIGHT));
	p1.individuals.setSpatialPosition(p1.pointUniform(p1.individualCount));
	community.rescheduleScriptBlock(s1, ticks=RUNTIME);
```
Finally, we define a reschedule script block to decide when to stop the simulation. In more sophisticated example, this would be when we define a map and additional reschedule script blocks for specific recording tasks.

## Mating (first() and reproduction())
At every tick, we make individuals find mate according to the second interaction type `i2`.
```
first() {
	// preparation for the reproduction() callback
	i2.evaluate(p1);
}
```
And then, we can tell what each individaul should do when `reproduction()` is called (happens implicitly). 
Here, we are making each individual draw one `mate` from their mating kernel. If only a mate exists (if `mate.size()` is 1), it adds new individuals via `addCrossed`. The number of individuals to add is defined by fecundity `FECUN` - drawn from a Poisson distribution with mean `FECUN` because the number of offspring has to be and integer. 
```
reproduction(){
	mate = i2.drawByStrength(individual, 1);
	if (mate.size())
		offspring = subpop.addCrossed(individual, mate, count=rpois(1, FECUN));
}
```
The locations of offspring will be determined in the next code block of `early()` stage.

## early() - disperse offspring and density control
By default non-Wright-Fisher architecture of SLiM, new offspring would have been generated via `reproduction` callback by `early` stage at each tick. 
So now, we need to give new individuals location. We do this my first choosing only the individuals with age zero, move them away from the parent by some dispersal distance sampled from dispersal kernel, using `pointDeviated`. 
Some individuals there are close to the edges can get new location that are outside the simulated area, so we need to tell `pointDeviated` what to do with them. We are using `reprising` option in this example, where it samples a new location until it lands within the boundary.
```
	// Disperse offspring
	offspring = p1.subsetIndividuals(maxAge = 0);
	pos = offspring.spatialPosition;
	pos = p1.pointDeviated(offspring.size(), pos, "reprising", INF, "n", SD);
	offspring.setSpatialPosition(pos);
```

Next, we need to set up local density based competition, so that the population density stays close to carrying capacity `K`.
We make each individual evaluate `i1` and use it to evaluate strength of local competition. 
Then each individual's fitness is rescaled according to Beverton-Holt model. 
The fitness in SLiM is the probability of an individual surviving to the next time step. 
So if the local density is low, fitness would be close to 1, and the individual is very likely to survive.
If the local density is much higher than `K`, the fitness will be close to zero, killing the individual before the next tick.

```
	// Measure local density and use it for density regulation
	i1.evaluate(p1);
	inds = p1.individuals;
	competition = i1.localPopulationDensity(inds);
	inds.fitnessScaling = 1 / (1 + RHO * competition);
```

## end of simulation - output
Simulation can end because there is no more individual to simulate. 
```
late() {
	if ((p1.individualCount == 0)){
		catn("Population went extinct! Ending the simulation");
		sim.simulationFinished();
	}
}
```
Or, the number of ticks past will become equal to RUNTIME. Then we save a tree sequence with PARAMS added as metadata.
```
s1 late() {
	catn("End of simulation (run time reached)");
	sim.treeSeqOutput(OUTPATH, metadata=PARAMS);
	sim.simulationFinished();
}
```
