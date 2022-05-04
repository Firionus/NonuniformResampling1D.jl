using SpecialFunctions

export HannBasis, RectangularBasis, KaiserBasis, TriangularBasis
export LanczosBasis

"""
    WindowFunction{F<:Function, T<:Real}

A radial window function. 

It is defined by the callback `_f` between x = 0 (center) and x = `width`. x=1
is the position of the nearest neighbor and the range x ∈ [0,1] is called a
unit.

# Fields

- `_f<:Function`  
Callback defining the window function. 

It will be called as `_f(x::Real)` and is expected to return a `Real` number.
`x` can be expected to go from 0 (center of window) to `width`. `x` values above
`width` will return 0 without calling `_f`. 

By convention, `_f(0) == 1`. If orthogonality to the nearest neighbor is
desired, it should return `_f(1) == 0`. 

    
- `width<:Real`

Value up to which the callback function will be evaluated. 
    
`width` must be bigger than 0. A `width` of 1 means that the window function is
evaluated up to the nearest neighbor. 
"""
struct WindowFunction{F<:Function, T<:Real}
    _f::F
    width::T

    function WindowFunction(_f::F, width::T) where {F<:Function, T<:Real}
        @assert width > 0 "WindowFunction width must be bigger than 0"
        new{F,T}(_f, width)
    end
end

function(w::WindowFunction)(x)
    x < 0 && error("undefined for negative values")
    x > w.width && return 0.0
    w._f(x)
end

RectangularBasis(width = 0.5) = WindowFunction(x -> 1.0, width)

HannBasis(width = 1.0) = WindowFunction(
    x -> cos(pi * x / 2 / width)^2, 
    width
)

KaiserBasis(width = 1.2; alpha = 3.5) = WindowFunction(
    x -> besseli(0, pi * alpha * sqrt(1 - (x / width)^2)) / besseli(0, pi * alpha),
    width
)

TriangularBasis(width=1.) = WindowFunction(
    x -> 1 - x/width,
    width
)

function LanczosBasis(lobes::Int = 3)
    @assert lobes > 0 "Lanczos Basis is only defined for at least one lobe"

    WindowFunction(
    x -> begin
        x == 1 && return 1.
        lobes*sin(pi*x)*sin(pi*x/lobes)/(pi^2*x^2)
    end,
    lobes
    )
end

# TODO adjust naming to Julia conventions

