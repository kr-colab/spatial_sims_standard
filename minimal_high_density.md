# `minimal_high_density.slim`: a minimal model suitable for high-density simulations

**code:** [minimal_high_density.slim](minimal_high_density.slim); on [github](https://github.com/kr-colab/spatial_sims_standard/blob/main/minimal_high_density.slim)

This code implements the same model as [`minimal.slim`](minimal.html),
but implements some "map-based" alternative computational methods so that runtime is linear in density
instead of quadratic.
However, these methods are somewhat approximate and may lead to discretization artifacts if density is low,
so if neighborhood sizes are small (below, say 5 or 10), then the methods should probably not be used.

There are two aspects of the "standard" methods demonstrated in `minimal.slim` that are quadratic in population density
(i.e., quadratic in the parameter `K`): local density computation, and mate choice.
The two methods can be used independently - for instance, maybe your simulation
has small interaction neighborhood size but large mating neighborhood size,
so you'd use these methods for mating but not density computations.

Both methods are linear because they rely on a pre-computed map of local density.
The basic computation to do this is
```
raw = summarizeIndividuals(p1.individuals, GRID_DIMS, p1.spatialBounds, operation="individuals.size();", perUnitArea=T);
```
This makes a grid of dimensions `GRID_DIMS`, counts up the number of individuals in each grid square,
and divides by the area of the square to get a density per unit area.
What is an appropriate size for these grid squares?
Well, we'd like these to be as big as is reasonable,
because computation will scale with the number of these squares.
But, this grid is the source of any discretization artifacts:
the smaller the grid squares are, the closer this model is to `minimal.slim`.
In both measuring local density and choosing mates we average over a kernel
with a certain characteristic scale: `SX` for local density (for "interactions") and `SM` for mating.
So, as long as our grid squares are sufficiently smaller than these,
any discretization will be smoothed out by those operations.
So, we set the size of the grid cells to be approximately `min(SX,SM)/2`:
```
grid_dims = asInteger(2 * (p1.spatialBounds[c(2,3)] - p1.spatialBounds[c(0,1)]) / min(SX, SM));
```

The `raw` values are put in two maps: `RAW_DENSITY` and `DENSITY`;
the first is used for mate choice and the second for interactions.

## Local smoothed density

Local density used as an input for mortality is computed in `minimal.slim`
by averaging over a Gaussian kernel with standard deviation `SX`.
To approximate this, we simply smooth the "raw" gridded density by that same kernel:
```
DENSITY.smooth(SX * 3, "n", SX);
```
and then look up the local, smoothed density at the location of every individual:
```
competition = p1.spatialMapValue(DENSITY, inds.spatialPosition);
```
This could be made more precise by taking into account that in computing the "raw" density
also involves some smoothing (the size of the grid squares).

## Mate choice

Mate choice in `minimal.slim` picks a nearby individual proportional to a weight
assigned from a Gaussian kernel with standard deviation `SM`.
To approximate this procedure, we choose a nearby location proportional to raw density
multiplied by the same Gaussian kernel,
and then return the nearest individual to that location:
```
mate_location = RAW_DENSITY.sampleNearbyPoint(individual.spatialPosition, 3 * SM, "n", SM);
mate = i2.nearestNeighborsOfPoint(mate_location, p1, 1);
```

## Possible pitfalls

The two lines in the "mate choice" snippet can each cause their own problems if care is not taken.

First, the call to `sampleNearbyPoint` works by rejection sampling:
it picks a location nearby from the provided kernel,
and then retains the point with probability proportional to the value of the `RAW_DENSITY` map at that point.
The way this is implemented means that if values of density in some parts of the map are much higher than in other parts of the map,
this may do a large number of rejections.
Concretely, if the density within a circle of radius `3 * SM` around a given individual
is at most, say $10^{-6}$ times the maximum density over the entire map,
then choosing a mate for that individual will require sampling millions of locations.
The solution to this might be biological: first, identify individuals with sufficent possible mates
to be able to mate successfully; and then only choose mates for those.

Second, the call to `i2.nearestNeighborsOfPoint` will fail, clearly, if there *are* no neighbors.
Since this call uses the interaction `i2`,
which was set up at the start of the simulation 
```
initializeInteractionType(2, "xy", reciprocal=T, maxDistance=5/sqrt(K));
```
to have maximum interaction distance `5/sqrt(K)`.
If the density is roughly $K$, this means there should be around 25 individuals in each circle of that radius;
however, if density is sufficiently nonuniform, or if you are simulating some nonequilibrium situation,
then this may fail.
The code is robust to this:
```
if (mate.size())
    subpop.addCrossed(individual, mate, count=rpois(1, FECUN));
```
means that if an individual has no neighbors, they will not reproduce,
but this is over a much smaller distance than our nominal mating distance, `SM`,
so we don't want this to happen very much.
(Indeed, if `SM` was not much larger than `5/sqrt(K)` then we wouldn't be saving any computation at all.)


