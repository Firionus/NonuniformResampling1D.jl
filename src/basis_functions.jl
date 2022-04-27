using SpecialFunctions

export FiniteBasisFunction, basis_value 
export HannBasis, RectangularBasis, KaiserBasis, TriangularBasis
export LanczosBasis

# TODO is this a good name? Kernel? Window?
struct FiniteBasisFunction{F<:Function,T<:Real}
    _f::F # callback x::Real -> y::Real, x ∈ [0, Inf), 
    # y(0) should probably be 1 (current node at x=0) 
    # y(1) should probably be 0 (next neighbor at x=1)
    width::T # ∈ (0, Inf), width up to which the basis function will be evaluated (inclusively)
end

function basis_value(basis::FiniteBasisFunction, x)
    x < 0 && error("undefined for negative values")
    x > basis.width && return 0.0
    basis._f(x)
end

RectangularBasis(width = 0.5) = FiniteBasisFunction(x -> 1.0, width)

HannBasis(width = 1.0) = FiniteBasisFunction(
    x -> cos(pi * x / 2 / width)^2, 
    width
)

KaiserBasis(width = 1.2; alpha = 3.5) = FiniteBasisFunction(
    x -> besseli(0, pi * alpha * sqrt(1 - (x / width)^2)) / besseli(0, pi * alpha),
    width
)

TriangularBasis(width=1.) = FiniteBasisFunction(
    x -> 1 - x/width,
    width
)

function LanczosBasis(lobes::Int = 3)
    @assert lobes > 0 "Lanczos Basis is only defined for at least one lobe"

    FiniteBasisFunction(
    x -> begin
        x == 1 && return 1.
        lobes*sin(pi*x)*sin(pi*x/lobes)/(pi^2*x^2)
    end,
    lobes
    )
end

# TODO adjust naming to Julia conventions

