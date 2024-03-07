# `adult_movement.slim`: when the organisms can move around

**code:** [adult_movement.slim](adult_movement.slim);
on [github](https://github.com/kr-colab/spatial_sims_standard/blob/main/adult_movement/adult_movement.slim)


The only difference between this script and [`minimal.slim`](../minimal.html)
is that in this script, all individuals move, every tick.
(In `minimal.slim`, offspring were placed at a location different from their parent's location,
but there was no explicit movement otherwise.)

To do this, we have this block of code:
```
    inds = p1.individuals;
    pos = inds.spatialPosition;
    pos = p1.pointDeviated(inds.size(), pos, "reflecting", INF, "n", SV);
    inds.setSpatialPosition(pos);
```
The call to `p1.pointDeviated()` is given a vector of positions
(the current locations of the individuals),
and returns a vector of new positions
that are randomly drawn from a Gaussian kernel centered on the previous positions.
drawn locations that fall outside the range of the population are reflected back into the range,
and the standard deviation of the Gaussian is `SV`.

This happens in a `first()` block, so it happens every tick (and, between mortality and reproduction).
