# Standardizing spatial simulations on SLiM

This is a repository of examples scripts to run SLiM simulations on (continuous) spatial landscapes.
The goal of these examples is to make it easy to explore various ways of setting up such simulations,
to choose parameters to roughly model a given system,
and to run computational experiments.

The [SLiM Manual](https://github.com/MesserLab/SLiM/releases/download/v4.1/SLiM_Manual.pdf)
is an excellent resource that covers many aspects of spatial simulations.
(To start off: sections 1.1-1.6, 15.1-15.2, and chapter 16.)
However, the manual tends to cover the *mechanics* of the tool;
to build a successful spatial simulation of a given system requries a fair amount of *craft* as well.
Questions that many people encounter very quickly are:
(1) How do I get my simulation to equilibrate at a given population size?
and (2) How do I vary (some aspect of the simulation) across a map?
Here, we aim to provide scripts that:

1. have adjustable parameters that correspond naturally to
     biologically important and observable quantities;
2. are simple, easy to understand, and easy to use for reproducible research; and
3. share as much structure as possible, so that parts from one script
    can be used in another script with minimal hassle.

*Getting started:*
To try something out, install SLiM, open the GUI,
change into this directory, open `minimal.slim`, and run it.
That's it! Then, change into some of the subdirectories and try those scripts:
they should all just run.
(For spatial models, the GUI is essential in developing and debugging;
computational experiments will run the scripts on the command line, of course.)

In the main directory, we have a minimal example (`minimal.slim`) of hermaphrodites without age-structure, adult movements, or a map.
The scripts in each folder follow the same format but have some elements that a more realistic simulation requires.

## The examples:

- [`minimal.slim`](minimal.html): a minimal example, that explains the general structure of the scripts.
- [`maps/`](maps/) (Use a simple map of a mountain to model heterogenous carrying capacity distribution in space) 
- [`adult_movement/`](adult_movement/) (Individuals continue to move around throughout its lifetime. Appropriate for animals, rather than plants)
- [`mate_choice/`](mate_choice/) (Dioecious population, offspring dispersed from female parents)
- [`selection/`](selection/) (Selective sweep in space)
- `case_studies/`
	- [`pikas/`](case_studies/pikas/) (pikas on a mountain with temperature rising)
	- [`toads/`](case_studies/toads/) (Austrailian cane toads' range expansion)
	- [`mosquito/`](case_studies/mosquito/) (mosquitos with a larval stage that lives only in rivers and seasonal fluctuations)
	- [`monarchs/`](case_studies/monarchs/) (monarch butterflies whose populations are regulated by a discrete set of host plants, that feeds caterpillars)

## Requirements

You need [SLiM](https://messerlab.org/SLiM) version at least 4.1
(and the GUI is *highly* recommended).
Python is used for set-up of some additional examples.
