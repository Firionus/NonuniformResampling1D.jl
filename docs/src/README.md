# NonuniformResampling1D.jl

Takes 1D data sampled at uniform intervals and resamples it at arbitrary,
nonuniform locations. 

## Installation

The package is currently not registered, so you'll have to install it from GitHub. Open a Julia REPL and run:

```julia
]add https://github.com/Firionus/NonuniformResampling1D.jl
```

## Getting Started

```jldoctest
using NonuniformResampling1D
xin = 1:9 # input locations
yin = [1, 3, 2, 6, 7, 3, 2, 8, 3] # input values
yout = nuresample(xin, yin, [4.5, 6.2])

# output
2-element Vector{Float64}:
 6.613040058882213
 2.7151900996408322
```

## Main API

```@docs
nuresample
```

## Window Functions

![Diagram Showing the Different Available Windows](examples/windows.svg)

```@docs
rect_window
hann_window
tri_window
kaiser_window
lanczos_window
```

### Custom Window Functions

To define your own window functions, take a look at the examples in [window_functions.jl](https://github.com/Firionus/NonuniformResampling1D.jl/blob/main/src/window_functions.jl) and use the exported `WindowFunction` type:

```@docs
WindowFunction
```

## Resampling Approach

The basic idea of the provided algorithm is that the output value at any
location should be given by a weighted average of input values in a region
around the output point, where the weight is given by a window function. The
supporting region is called a slice and is proportional in size to the distance
to the next neighbor. The ratio between slice width and distance to the nearest
neighbor is the window function width and determines how strong the smoothing
is. 

Further, if there are few points in a slice, the result might not be accurate.
Therefore, ad-hoc resampling is used. So if the left or right slice of an output
point contains less points than specified, upsampling on a regular grid is
performed for both slices. 

Consider an example with 
```julia
nuresample(xin, yin, xout, hann_window(1.35), required_points_per_slice = 4)
```
which can be visualized like this:

![Diagram Showing the Resampling Approach](examples/explanation_diagram.svg)

The orange output point shows how the window function and supporting regions are
scaled asymmetrically around the output point, depending on how far the
neighboring output points are away. The input values and the window function are
multiplied pointwise, resulting in the windowed slices shown in the second plot
from the top. The normalized average of these points is the orange output value. 

The blue point shows the dynamic upsampling. since its slices contain less than
4 points, upsampling is applied before applying the window function. The default
upsampling function Lanczos3 is used to create the points in bright blue. They
are created with uniform step size. The upsampled points are then weighted with
the window function and the normalized average forms the blue output value. 

Upsampling is performed for both slices, even if only one of them has too few
points. This ensures that the weight between the two slices stays approximately
even. The upsampling step is the same between the two slices and is chosen such
that `required_points_per_slice` are created in the smaller of the two slices.

## Status

This package is in its youth and some hiccups, like occasional errors or
suboptimal performance, are expected. Please see the
[Issues](https://github.com/Firionus/NonuniformResampling1D.jl/issues) for
details and report any problems you experience. 