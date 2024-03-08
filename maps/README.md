# `simple_map.slim`: when space isn't flat

**code:** [simple_map.slim](simple_map.slim);
on [github](https://github.com/kr-colab/spatial_sims_standard/blob/main/maps/simple_map.slim)

This is a simulation very similar to [`minimal.slim`](../minimal.html):
the only difference is that in that script, local "carrying capacity"
was represented by a single value, `K`;
while in this script, the value of `K` varies across space
(and is specified in the simulation by passing in a PNG file).
The differences are as follows.

First, we read in and set up the map.
The pixels in the PNG can be represented either as integers between 0 and 255
or as floats between 0 and 1.
Here we take the `floatK` values, which provide the intensity values of the each pixel as floats
(here "K" is as in [CMYK](https://en.wikipedia.org/wiki/CMYK_color_model), not the carrying capacity);
to create the values of the map, we rescale these to lie between `min_l` and `max_l`:
```
	p1.individuals.setSpatialPosition(p1.pointUniform(p1.individualCount));
	mapImage = Image(image_location); // adding the image
	map_vals = (1-mapImage.floatK) *(max_l-min_l) + min_l; // defining values
	map = p1.defineSpatialMap("world", "xy", map_vals, valueRange=c(min_l, max_l), colors=c('#FFFFFF', '#000000'));
	defineConstant("MAP", map);
```

Then, in the block where we determine fitness
(which here is the probability of survival to the next time stop),
we look up the value of the map at the locations of every individual:
```
	local_l = MAP.mapValue(inds.spatialPosition);
	inds.fitnessScaling = 1/(1 + (RHO * competition / local_l));
```
