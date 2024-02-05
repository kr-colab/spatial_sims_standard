# Standardizing spatial simulations on SLiM
We demonstrate several spatial simulation scripts that can (hopefully) be modified easily by users for their ecological scenarios in mind.

List of examples:
- selection
- maps
- adult movement
- mate choice
- case studies

In the main directory, we have a minimal example of hermaphrodites without age-structure, adult movements, or a map.
Other scripts in subdirectory will follow the format of the minimal example but will have some elements that a more realistic simulation requires.

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

Then we set up default parameter values including ....
```
	defaults = Dictionary(
		"seed", getSeed(),
		"SD", 0.3, // sigma_D, dispersal distance
		"SI", 0.3, // sigma_I, interaction distance for competition
		"SM", 0.3, // sigma_M, mate choice distance
		"K", 5, // carrying capacity per unit square
		"LIFETIME", 4, // average life span
		"WIDTH", 25.0, // width of the simulated area
		"HEIGHT", 25.0, // height of the simulated area
		"RUNTIME", 200, // total number of ticks to run the simulation for
		"L", 1e8, // genome length
		"R", 1e-8, // recombination rate (per tick)
		"MU", 1e-8, // mutation rate (per tick)
		"OUTDIR", exists("OUTDIR") ? OUTDIR else ".",
		"PARAMFILE", exists("PARAMFILE") ? PARAMFILE else "./params.json"
		);
```
To replace some of the values defined as default, we can use a JSON file instead of modifying the slim script directly. 
In the next line, we let the parameter values to be overwritten by whatever we write in `params.json` in working directory.
```
	if (fileExists(defaults.getValue("PARAMFILE"))){
		local_defaults = Dictionary(paste(readFile(defaults.getValue("PARAMFILE")), sep="\n"));
		defaults.addKeysAndValuesFrom(local_defaults);
		defaults.setValue("read_from_paramfile", defaults.getValue("PARAMFILE"));
	}
```
We set where the output file should be saved and use the random seed to name the file.
```
	defaults.setValue("OUTBASE", defaults.getValue("OUTDIR") + "/out_" + defaults.getValue("seed"));
	defaults.setValue("OUTPATH", defaults.getValue("OUTBASE") + ".trees");
```
## reproduction
## 1 first
## first - mate choice

## early - competition and density control
## end of simulation - output

