# hydro-minilab

This repo contains a few useful modules for solving transport problems. 

## Case study: d=1 

Simple advection is prone to numerical shortcomings (such as instability/diffusion). 
Therefore some special techniques are required to manage **edge cases**. We use the
 [shasta](http://dx.doi.org/10.1016/0021-9991(73)90147-2) scheme, which caters
 particularly well for steep wave-fronts.

We study first the incompressible advection equation (uniform velocity). This serves
as the first consistency check. We then solve a one dimensional hydrodynamics problem,
with a radially symmetry geometry.

## Case study: d=2

TODO, off-central collisions.

## visualisation

we use the [matplotlib](http://matplotlib.org) package for _python_. 
