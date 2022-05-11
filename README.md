
<a id='NonuniformResampling1D.jl'></a>

<a id='NonuniformResampling1D.jl-1'></a>

# NonuniformResampling1D.jl


Takes 1D data sampled at uniform intervals and resamples it at arbitrary, nonuniform locations. 


<a id='Installation'></a>

<a id='Installation-1'></a>

## Installation


The package is currently not registered, so you'll have to install it from GitHub. Open a Julia REPL and run:


```julia
]add https://github.com/Firionus/NonuniformResampling1D.jl
```


<a id='Getting-Started'></a>

<a id='Getting-Started-1'></a>

## Getting Started


```julia
using NonuniformResampling1D
xin = 1:9 # input locations
yin = [1, 3, 2, 6, 7, 3, 2, 8, 3] # input values
yout = nuresample(xin, yin, [4.5, 6.2])

# output
2-element Vector{Float64}:
 6.613040058882213
 2.7151900996408322
```


<a id='Main-API'></a>

<a id='Main-API-1'></a>

## Main API

<a id='NonuniformResampling1D.nuresample' href='#NonuniformResampling1D.nuresample'>#</a>
**`NonuniformResampling1D.nuresample`** &mdash; *Function*.



```julia
nuresample(xin, yin, xout, smoothing_function; kwargs...)
```

Take data `yin` sampled at uniformly spaced locations `xin` and resample at nonuniform locations `xout`. 

**Positional Arguments**

  * `xin` is an `AbstractRange` of input locations
  * `yin` are the input values
  * `xout` are the increasing output locations
  * `smoothing_function` is a [`WindowFunction`](README.md#NonuniformResampling1D.WindowFunction). The width of the window function is scaled proportional to the distance to the nearest neighbor on the left and right side of the output point. At the start and end of `xout` a symmetrical window is used. The default is a rectangular window without overlap.

**Keyword Arguments**

  * `upsampling_function`: a [`WindowFunction`](README.md#NonuniformResampling1D.WindowFunction) that is used to perform upsampling when the density of input points is not high enough compared to the density of output points. The default is Lanczos3.
  * `required_points_per_slice`: the number of input points that must be in a slice (left or right half of a window) before smoothing. If the number of input points is lower than this value, upsampling is performed to create `required_points_per_slice` many points before applying the smoothing function.   The default is `round(4 * smoothing_function.width)`, which means 4 points are required per unit width to the nearest neighbor.   The minimum value is 1, in which case upsampling is only performed if no input point falls in a slice.

**Examples**

```julia-repl
julia> nuresample(1:9, 1:9, [4.2, 6.2])
2-element Vector{Float64}:
 4.201106036510037
 6.201106036510037
```

**Output Type**

Returns an `Array{Float64, 1}`. Other output types are currently unsupported. 


<a target='_blank' href='https://github.com/Firionus/NonuniformResampling1D.jl/blob/9d3691b448593da59a0965f7787fddf533f4c765/src/NonuniformResampling1D.jl#L9-L53' class='documenter-source'>source</a><br>


<a id='Window-Functions'></a>

<a id='Window-Functions-1'></a>

## Window Functions

<a id='NonuniformResampling1D.rect_window' href='#NonuniformResampling1D.rect_window'>#</a>
**`NonuniformResampling1D.rect_window`** &mdash; *Function*.



```julia
rect_window(width = 0.5)
```

Rectangular window which is 1 from 0 to `width` and 0 otherwise. 

Can be used to calculate a moving average, for example.


<a target='_blank' href='https://github.com/Firionus/NonuniformResampling1D.jl/blob/9d3691b448593da59a0965f7787fddf533f4c765/src/window_functions.jl#L5-L11' class='documenter-source'>source</a><br>

<a id='NonuniformResampling1D.hann_window' href='#NonuniformResampling1D.hann_window'>#</a>
**`NonuniformResampling1D.hann_window`** &mdash; *Function*.



```julia
hann_window(width = 1)
```

Hann window (raised cosine) that is 1 at 0 and reaches 0 at `width`.


<a target='_blank' href='https://github.com/Firionus/NonuniformResampling1D.jl/blob/9d3691b448593da59a0965f7787fddf533f4c765/src/window_functions.jl#L14-L18' class='documenter-source'>source</a><br>

<a id='NonuniformResampling1D.tri_window' href='#NonuniformResampling1D.tri_window'>#</a>
**`NonuniformResampling1D.tri_window`** &mdash; *Function*.



```julia
tri_window(width = 1)
```

Triangular window which is 1 at 0 and reaches 0 at `width`. 

Can be used to perform linear interpolation when used as upsampling_function with `width=1`. 


<a target='_blank' href='https://github.com/Firionus/NonuniformResampling1D.jl/blob/9d3691b448593da59a0965f7787fddf533f4c765/src/window_functions.jl#L34-L41' class='documenter-source'>source</a><br>

<a id='NonuniformResampling1D.kaiser_window' href='#NonuniformResampling1D.kaiser_window'>#</a>
**`NonuniformResampling1D.kaiser_window`** &mdash; *Function*.



```julia
kaiser_window(width = 1.2; alpha = 3.5)
```

Kaiser window cut off at `width` with the shape parameter `alpha`.


<a target='_blank' href='https://github.com/Firionus/NonuniformResampling1D.jl/blob/9d3691b448593da59a0965f7787fddf533f4c765/src/window_functions.jl#L24-L28' class='documenter-source'>source</a><br>

<a id='NonuniformResampling1D.lanczos_window' href='#NonuniformResampling1D.lanczos_window'>#</a>
**`NonuniformResampling1D.lanczos_window`** &mdash; *Function*.



```julia
lanczos_window(lobes = 3; width = 1)
```

Lanczos window which is 1 at 0 and 0 at 1*`width`, 2*`width`, 3*`width`, ...

The higher the number of lobes, the higher the accuracy of sinc approximation. 


<a target='_blank' href='https://github.com/Firionus/NonuniformResampling1D.jl/blob/9d3691b448593da59a0965f7787fddf533f4c765/src/window_functions.jl#L47-L53' class='documenter-source'>source</a><br>


<a id='Window-Function-Type'></a>

<a id='Window-Function-Type-1'></a>

# Window Function Type


To define your own window functions, you can construct your own 

<a id='NonuniformResampling1D.WindowFunction' href='#NonuniformResampling1D.WindowFunction'>#</a>
**`NonuniformResampling1D.WindowFunction`** &mdash; *Type*.



```julia
WindowFunction{F<:Function, T<:Real}
```

A radial window function. 

It is defined by the callback `_f` between x = 0 (center) and x = `width`. x=1 is the position of the nearest neighbor and the range x ∈ [0,1] is called a unit.

**Fields**

  * `_f<:Function`

Callback defining the window function. 

It will be called as `_f(x::Real)` and is expected to return a `Real` number. `x` can be expected to go from 0 (center of window) to `width`. `x` values above `width` will return 0 without calling `_f`. 

By convention, `_f(0) == 1`. If orthogonality to the nearest neighbor is desired, it should return `_f(1) == 0`. 

  * `width<:Real`

Value up to which the callback function will be evaluated. 

`width` must be bigger than 0. A `width` of 1 means that the window function is evaluated up to the nearest neighbor. 


<a target='_blank' href='https://github.com/Firionus/NonuniformResampling1D.jl/blob/9d3691b448593da59a0965f7787fddf533f4c765/src/WindowFunction.jl#L3-L31' class='documenter-source'>source</a><br>


<a id='Resampling-Approach'></a>

<a id='Resampling-Approach-1'></a>

## Resampling Approach


The basic idea of the provided algorithm is that the output value at any point should be given by a weighted average of input points in a region around the output point. This region is called a slice and is proportional in size to the distance to the next neighbor. 


TODO


TODO Insert Diagramm


<a id='Status'></a>

<a id='Status-1'></a>

## Status


This package is in its youth and some hickups, like occasional errors or suboptimal performance, are expected. Please see the [Issues](https://github.com/Firionus/NonuniformResampling1D.jl/issues) for details and report any problems you experience. 

