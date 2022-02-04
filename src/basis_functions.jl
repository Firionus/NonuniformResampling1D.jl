struct FiniteBasisFunction{F<:Function, T<:Real}
    f::F # callback x::Real -> y::Real, x ∈ [0, Inf), 
    # y(0) should probably be 1 (current node at x=0) 
    # y(1) should probably be 0 (next neighbor at x=1)
    width::T # ∈ (0, Inf), width up to which the basis function will be evaluated (inclusively)
end

HannBasis(width=1.) = FiniteBasisFunction(x -> begin
    x < 0 && error("undefined for negative values")
    x > width && return 0.
    # normalized integral would require pre-factor 2/width
    cos(pi*x/2/width)^2
end, width)

RectBasis(width=.5) = FiniteBasisFunction(x -> begin
    x < 0 && error("undefined for negative values")
    x > width && return 0.
    # normalized integral would require pre-factor 2/width
    1.
end, width)

# TODO put Kaiser back in, even if it performs badly, it's nice to have the choice for experimentation

# TODO triangle_basis

# TODO adjust naming to Julia conventions

