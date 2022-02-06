export FiniteBasisFunction, HannBasis, RectBasis

struct FiniteBasisFunction{F<:Function, T<:Real}
    _f::F # callback x::Real -> y::Real, x ∈ [0, Inf), 
    # y(0) should probably be 1 (current node at x=0) 
    # y(1) should probably be 0 (next neighbor at x=1)
    width::T # ∈ (0, Inf), width up to which the basis function will be evaluated (inclusively)
end

function val(basis::FiniteBasisFunction, x)
    x < 0 && error("undefined for negative values")
    x > basis.width && return 0.
    basis._f(x)
end

HannBasis(width=1.) = FiniteBasisFunction(x -> cos(pi*x/2/width)^2, width)

RectBasis(width=.5) = FiniteBasisFunction(x -> 1., width)

# TODO put Kaiser back in, even if it performs badly, it's nice to have the choice for experimentation

# TODO triangle_basis

# TODO adjust naming to Julia conventions

