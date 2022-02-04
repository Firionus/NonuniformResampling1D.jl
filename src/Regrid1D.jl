module Regrid1D

using SpecialFunctions, Interpolations

export regrid, HannBasis, KaiserBasis

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

function kaiser(x, alpha, width) 
    besseli(0, pi*alpha*sqrt(1-(x/width)^2))/besseli(0, pi*alpha)
end

function KaiserBasis(width = 1.2; alpha = 3.5, samples=64) 
    if samples < 2 # no approximation, e.g. with samples = 0
        f = x -> kaiser(x, alpha, width)
    else
        # approximate for higher speed with linearly interpolated LUT
        # the default of 64 samples results in errors of approximately < 0.3e-3
        x = range(0, width, length=samples)
        y = kaiser.(x, alpha, width)
        itp = interpolate(y, BSpline(Linear()))
        sitp = scale(itp, x)
        f = x -> sitp(x) 
        # evaluating these interpolations is somehow slow, as kaiser takes 10x longer than Hann. 
        # But overall, that's not a big issue probably, since the difference between Hann and Kaiser is negligible when looking at the results. 
        # So, we can probably take out the Kaiser window without much problems for my application
    end

    FiniteBasisFunction(x -> begin
    x < 0 && error("undefined for negative values")
    x > width && return 0.
    f(x)
    end, width)
end

function find_first_above_or_equal(cutoff, x::StepRangeLen)
    # this doesn't work for all cases - assert the right ones
    @assert Float64(x.step) > 0
    @assert x.offset == 1
    # basic math:
    # y = x.ref + (i-x.offset)*x.step
    # i = (y - x.ref)/x.step + x.offset
    ind_float = (cutoff - x.ref)/x.step
    ind = Int(ceil(Float64(ind_float))) + x.offset
    @assert ind <= x.len # none of the elements are above or equal the cutoff
    return max(1, ind)
end

function find_last_below_or_equal(cutoff, x::StepRangeLen)
    # this doesn't work for all cases - assert the right ones
    @assert Float64(x.step) > 0
    @assert x.offset == 1
    # basic math:
    # y = x.ref + (i-x.offset)*x.step
    # i = (y - x.ref)/x.step + x.offset
    ind_float = (cutoff - x.ref)/x.step
    ind = Int(floor(Float64(ind_float))) + x.offset
    @assert ind >= 1 # none of the elements are below or equal the cutoff
    return min(x.len, ind)
end

function regrid(xin::StepRangeLen, yin, xout,
    smoothing_function::FiniteBasisFunction = HannBasis(1.1),
    )

    # extent of smoothing function
    extent = smoothing_function.width

    step = xin.step|>Float64

    # allocate
    yout = Array{Float64, 1}(undef, length(xout))


    out_ind = 1


    # virtual first slice
    # TODO document that this results in nearest neighbor-like behavior when in_resolution < out_resolution
    slice_width = max(step, xout[2] - xout[1]) 
    slice_start = xout[1] - slice_width
    slice_stop = xout[1]


    while true
        # start with left slice
        input_start = slice_stop - slice_width*extent
        input_stop = slice_stop
        @assert input_start >= xin[1]
        @assert input_stop <= xin[end]

        #find relevant input points
        input_start_ind = find_first_above_or_equal(input_start, xin) # TODO if we always select or equal, we'll select values exactly on the output point twice and weigh them too heavily
        input_stop_ind = find_last_below_or_equal(input_stop, xin)
        # require at least one point in the window
        @assert input_stop_ind >= input_start_ind

        val_acc = 0.
        win_acc = 0.
        for input_ind in input_start_ind:input_stop_ind
            x = xin[input_ind]
            rel_pos = (input_stop - x)/slice_width
            win_value::Float64 = smoothing_function.f(rel_pos)
            win_acc += win_value
            val_acc += yin[input_ind]*win_value
        end
        # assing on left slice
        yout[out_ind] = val_acc / win_acc / 2 # normalization

        # left slice done
        # now right slice
        if out_ind < xout|>length # keep previous slice width on last element
            slice_width = max(step, xout[out_ind + 1] - xout[out_ind])
        end
        slice_start = xout[out_ind]
        slice_stop = xout[out_ind] + slice_width

        input_start = slice_start
        input_stop = slice_start + slice_width*extent
        @assert input_start >= xin[1]
        @assert input_stop <= xin[end]

        input_start_ind = find_first_above_or_equal(input_start, xin)
        input_stop_ind = find_last_below_or_equal(input_stop, xin)
        # require at least one point in the window
        @assert input_stop_ind >= input_start_ind

        val_acc::Float64 = 0.
        win_acc::Float64 = 0.
        for input_ind in input_start_ind:input_stop_ind
            x = xin[input_ind]
            rel_pos = (x - input_start)/slice_width
            win_value::Float64 = smoothing_function.f(rel_pos)
            win_acc += win_value
            val_acc += yin[input_ind]*win_value
        end
        # add on right slice
        yout[out_ind] += val_acc / win_acc / 2 # normalization

        # right slice done
        # prepare next left slice

        out_ind += 1
        out_ind > length(xout) && break

        slice_width = max(step, xout[out_ind] - xout[out_ind - 1])
        slice_start = xout[out_ind] - slice_width
        slice_stop = xout[out_ind]
    end

    yout
end

end # module
