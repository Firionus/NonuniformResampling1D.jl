module Regrid1D

export regrid1d, HannBasis

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

function regrid1d(xin, yin, xout; 
    min_upsampling = 4, 
    min_downsampling = 4,
    upsampling_function::FiniteBasisFunction = HannBasis(1.),
    downsampling_funciton::FiniteBasisFunction = HannBasis(1.),
    )
    # allocate
    yout = Array{Float64, 1}(undef, length(xout))

    # TODO upsample input by min_upsampling with upsampling_function (own function)

    # TODO iterate over the spaces between output points, 
    # taking the average between the left and right space to each output point

    # in each inter-space check that sample rate is higher than min_downsampling
    # otherwise upsample again with upsampling_function
    # then take mean of points weighted with downsampling_function

    println("hello world - this is WIP")
end

end # module
