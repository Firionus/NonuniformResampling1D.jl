module Regrid1D

include("basis_functions.jl")
export HannBasis, RectBasis

export regrid

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

function find_last_below(cutoff, x::StepRangeLen)
    # this doesn't work for all cases - assert the right ones
    @assert Float64(x.step) > 0
    @assert x.offset == 1
    # basic math:
    # y = x.ref + (i-x.offset)*x.step
    # i = (y - x.ref)/x.step + x.offset
    ind_float = (cutoff - x.ref)/x.step
    ind = Int(floor(Float64(ind_float)|>prevfloat)) + x.offset
    @assert ind >= 1 # none of the elements are below or equal the cutoff
    return min(x.len, ind)
end

function regrid(xin::StepRangeLen, yin, xout,
    smoothing_function::FiniteBasisFunction = RectBasis(.5),
    )
    in_step = Float64(xin.step)
    # allocate
    yout = Array{Float64, 1}(undef, length(xout))

    # first slice
    right_slice_width = max(in_step, xout[2] - xout[1])
    left_slice_width = right_slice_width
    
    out_ind = 1
    
    while true
        yout[out_ind] = interpolate_point(xin, yin, xout[out_ind], left_slice_width, right_slice_width, smoothing_function)

        # prepare next point
        out_ind += 1
        out_ind > length(xout) && break

        if out_ind < length(xout) # keep previous slice width on last element
            right_slice_width = max(in_step, xout[out_ind + 1] - xout[out_ind])
        end

        left_slice_width = max(in_step, xout[out_ind] - xout[out_ind - 1])
    end

    yout
end

function interpolate_point(xin, yin, xpoint, left_unit_width, right_unit_width, basis::FiniteBasisFunction)
    # start with left slice
    slice_width = left_unit_width

    input_width = slice_width*basis.width
    input_start = xpoint - input_width
    input_stop = xpoint
    @assert input_start >= xin[1]
    @assert input_stop <= xin[end] "Not enough points at the end of the input"

    #find relevant input points
    input_start_ind = find_first_above_or_equal(input_start, xin) 
    input_stop_ind = find_last_below(input_stop, xin)
    # require at least one point in the window
    @assert input_stop_ind >= input_start_ind

    input_x = xin[input_start_ind:input_stop_ind]
    input_y = @view yin[input_start_ind:input_stop_ind]
    left_slice_contribution = .5*slice_weighted_mean(xpoint, slice_width, input_x, input_y, basis)

    # right slice
    slice_width = right_unit_width

    input_width = slice_width*basis.width
    input_start = xpoint
    input_stop = xpoint + input_width
    @assert input_start >= xin[1]
    @assert input_stop <= xin[end] "Not enough points at the end of the input"

    #find relevant input points
    input_start_ind = find_first_above_or_equal(input_start, xin)
    input_stop_ind = find_last_below_or_equal(input_stop, xin)
    # require at least one point in the window
    @assert input_stop_ind >= input_start_ind "input range $(input_start) to $(input_stop) does not contain at least one point"

    input_x = xin[input_start_ind:input_stop_ind]
    input_y = @view yin[input_start_ind:input_stop_ind]
    right_slice_contribution = .5*slice_weighted_mean(xpoint, slice_width, input_x, input_y, basis)
    
    return left_slice_contribution + right_slice_contribution
end

function slice_weighted_mean(xpoint, slice_width, xin_points, yin_points, basis)
    val_acc::Float64 = 0.
    win_acc::Float64 = 0.
    for i in eachindex(xin_points)
        x = xin_points[i]
        y = yin_points[i]
        rel_pos = abs(x - xpoint)/slice_width
        win_value::Float64 = basis.f(rel_pos)
        win_acc += win_value
        val_acc += y*win_value
    end
    return val_acc / win_acc # normalization
end

end # module
