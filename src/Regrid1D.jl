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

# Rough layout of algorithm
# =========================

# determine slice (a slice is a space between two output nodes)
# needed for slice: xstart, xstop (with xstop > xstart (check this condition holds!)), out_ind_start, out_ind_stop
# at the beginning and end: assume slice is as big as next one

# determine involved input nodes (throw if non-existent, don't assume any boundary condition)

# determine linear grid needed for slice
# slice_grid_step = min(slice_width/min_downsampling, min(involved_input_nodes_steps)/min_upsampling)

# high resolution linear grid x-positions
# starting 1/2 step next to xstart

# resample to high resolution linear grid with upsampling_function

# for left node of slice
# weight high resolution points with downsampling_function and accumulate while also accumulating the sum of the weights
# normalize sum such that sum of weight is .5
# add to previous value of left node

# for right node of slice: as left node, but assign instead of adding to previous value

# repeat until we're through



# TODO document requirement for xin/xout to be ascending
function regrid1d(xin, yin, xout; 
    min_upsampling = 4, 
    min_downsampling = 4,
    upsampling_function::FiniteBasisFunction = HannBasis(1.),
    downsampling_funciton::FiniteBasisFunction = HannBasis(1.),
    )
    # allocate
    yout = Array{Float64, 1}(undef, length(xout))

    # virtual first slice
    slice_width = xout[2] - xout[1]
    xstart = xout[1] - slice_width
    xstop = xout[1]
    out_ind_start = 0 # doesn't exist
    out_ind_stop = 1

    while true
        # determine input nodes involved in slice 
        
        # damn, this is gonna be really hard because how do we know there isn't
        # some point out really far that has a large distance to its neighbor
        # and therefore creates a far-reaching effect all the way to us? (for
        # basis width > 1)


    end

    println("hello world - this is WIP")
end

end # module
