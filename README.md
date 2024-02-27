# Standardizing spatial simulations on SLiM
We demonstrate several spatial simulation scripts that can (hopefully) be modified easily by users for their ecological scenarios in mind.

In the main directory, we have a minimal example of hermaphrodites without age-structure, adult movements, or a map.

List of subdirectories:
- selection
- maps
- adult movement
- mate choice
- case studies

The scripts in each folder closely follow the format of `minimal.slim` but will have some elements that a more realistic simulation requires.


## Setup - install SLiM
You need SLiM 4.1 to run `.slim` scripts from this repository. Download the latest SLiM and follow the installation instructions for your system (https://messerlab.org/slim/)

## Initialize - set up parameters
In `initialize() {...}`, we first declare that the models is non Wright-Fisher and is in two dimension (can be switched to one or three dimension depending on the model).
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

We use a user-defined function, `setupParams` to set up parameter values based on `default` and a JSON file `PARAMFILE`. The `setupParams` is defined at the end of `minimal.slim` from line 94. `PARAMFILE` overwrites any parameter values in `default`. 

TODO : EDIT FROM HERE.

The average fecundity is expected to be 1 / lifetime once the population equilibrates. Therefore, we define FECUN using LIFETIME we defined earlier in the for loop going through all keys of defaults.
```
	defineConstant("FECUN", 1/LIFETIME);
```
FECUN is used to define another constant RHO, which will be used for fitness-rescaling according to Beverton-Holt model later in the script. See our paper (specific section) for more mathematical discussion on Beverton-Holt model and other density-control models.
```
	defineConstant("RHO", FECUN/((1+FECUN)*K));
```
We set the seed for the simulation. If we repeat this simulation with the same seed, we will see the same result, but if we change it, the result may look slightly different due to randomness.
```
	setSeed(seed);
```
The next few lines define properties of the genome of each individual. This step is shared across all SLiM simulations, so check out the official SLiM manual to explore more options. In this example, the entire genome shares the same recombination rate and only experiences neutral mutations. 
```
	initializeMutationRate(MU);
	initializeMutationType("m1", 0.5, "f", 0.0);
	initializeGenomicElementType("g1", m1, 1.0);
	initializeGenomicElement(g1, 0, L-1);
	initializeRecombinationRate(R);
```
Lastly, we need to define two interactions - `i1` for competition and `i2` for mate choice. In this example, we use a normal distribution with standard deviation `SI` or `SM`, but one can use a different shaped kernel which defines how the strength of interaction decays as a function of distance between two interacting individuals. We also set maximum distance between two individuals to compete or mate, but that is optional. 
```
	initializeInteractionType(1, "xy", reciprocal=T, maxDistance=3*SI);
	i1.setInteractionFunction("n", 1.0/sqrt(2*PI*SI^2), SI);
	
	initializeInteractionType(2, "xy", reciprocal=T, maxDistance=3*SM);
	i2.setInteractionFunction("n", 1.0/sqrt(2*PI*SM^2), SM);
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
```
	sim.addSubpop("p1", asInteger(K * WIDTH * HEIGHT));
	p1.setSpatialBounds(c(0, 0, WIDTH, HEIGHT));
	p1.individuals.setSpatialPosition(p1.pointUniform(p1.individualCount));
	community.rescheduleScriptBlock(s1, ticks=RUNTIME);
```
We begin the simulation by adding `K * WIDTH * HEIGHT` individuals to subpopulation `p1`. 
We also define the boundaries of the simulated area.
Then we place the initial population in space uniformly in the simulated area.
Finally, we define a reschedule script block to decide when to stop the simulation.
## first - mate choice
```
	i2.evaluate(p1);
```
Every tick, we first make every individual to evaluate `i2` to find their potential mate.
## early - competition and density control
By default non-Wright-Fisher architecture of SLiM, new offsprings would have been generated via `reproduction` callback even though we don't explicitly say it in this script. But the density control has to be done explicitly.
```
	i1.evaluate(p1);
	inds = p1.individuals;
	competition = i1.localPopulationDensity(inds);
	inds.fitnessScaling = 1/(1+RHO*competition);
```
We make each individual evaluate `i1` and use it to evaluate strength of local competition. 
Then each individual's fitness is rescaled according to Beverton-Holt model. Internally, the fitness is the probability of survival to the next tick.
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
