# Regrid1D.jl

Resample 1-dimensional data from arbitrary grids onto arbitrary grids in Julia.

## Status

**Work in Progress**: The package is not functional yet. 

If you are interested in collaborating, please open an issue. 

## Goals

- High quality resampling from arbitrary 1D grids to arbitrary 1D grids
- Specific application: Resampling from linear to log grid for magnitude of FFT
- Performance is secondary right now

## Approach

- Use distance to next grid point as a kind of "local sample rate" with which to
  scale basis functions for interpolation/smoothing. Treat left and right
  neighbor separately and average between them. 
- First perform upsampling by adding basis functions from each input point.
  Require a certain amount of upsampling in any case. 
- If output points are closer than input points, do even more upsampling
- Then create output point values by weighted average where the smoothing
  function is the weight