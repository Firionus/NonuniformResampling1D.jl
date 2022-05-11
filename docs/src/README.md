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

The basic idea of the provided algorithm is that the output value at any point should be given by a weighted average of input points in a region around the output point. This region is called a slice and is proportional in size to the distance to the next neighbor. 

TODO

TODO Insert Diagramm

## Status

This package is in its youth and some hickups, like occasional errors or
suboptimal performance, are expected. Please see the
[Issues](https://github.com/Firionus/NonuniformResampling1D.jl/issues) for
details and report any problems you experience. 