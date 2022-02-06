module Regrid1D

include("basis_functions.jl")
include("range_utilities.jl")

export regrid

function regrid(xin::StepRangeLen, yin, xout,
    smoothing_function::FiniteBasisFunction = RectBasis(.5);
    required_input_points=4, upsampling_basis=HannBasis(2.) # TODO change upsampling default to Lanczos Basis
    )
    # allocate
    yout = Array{Float64, 1}(undef, length(xout))

    # first slice
    right_slice_width = xout[2] - xout[1]
    left_slice_width = right_slice_width
    
    out_ind = 1
    
    while true
        yout[out_ind] = interpolate_point(xin, yin, xout[out_ind], left_slice_width, right_slice_width, 
        smoothing_function, 
        required_input_points=required_input_points, 
        upsampling_basis=upsampling_basis)

        # prepare next point
        out_ind += 1
        out_ind > length(xout) && break

        if out_ind < length(xout) # keep previous slice width on last element
            right_slice_width = xout[out_ind + 1] - xout[out_ind]
        end

        left_slice_width = xout[out_ind] - xout[out_ind - 1]
    end

    yout
end

function interpolate_point(xin, yin, xpoint, left_unit_width, right_unit_width, basis::FiniteBasisFunction;
    required_input_points=1, upsampling_basis=missing,
    )
    # start with left slice
    slice_width = left_unit_width

    input_width = slice_width*basis.width
    input_start = xpoint - input_width
    input_stop = xpoint
    @assert input_start >= xin[1] "Not enough points at the beginning of the input"
    @assert input_stop <= xin[end] "Not enough points at the end of the input"

    #find relevant input points
    input_start_ind = find_first_above_or_equal(input_start, xin) 
    input_stop_ind = find_last_below(input_stop, xin)

    input_points = input_stop_ind - input_start_ind + 1
    if input_points < required_input_points
        #upsample
        @assert !ismissing(upsampling_basis) "Not enough input points and no upsampling basis"
        upsample_step = input_width/required_input_points
        input_x = range(input_start + upsample_step/2, step=upsample_step, length=required_input_points)
        in_step = Float64(xin.step)
        input_y = [interpolate_point(xin, yin, up_x, in_step, in_step, upsampling_basis) for up_x in input_x] 
        # TODO use buffer for upsampled values (length == required_input_values) that is allocated at a call to `regrid` instead of allocating every time
    else # don't upsample
        input_x = xin[input_start_ind:input_stop_ind]
        input_y = @view yin[input_start_ind:input_stop_ind]
    end

    left_slice_contribution = .5*slice_weighted_mean(xpoint, slice_width, input_x, input_y, basis)

    # right slice
    slice_width = right_unit_width

    input_width = slice_width*basis.width
    input_start = xpoint
    input_stop = xpoint + input_width
    @assert input_start >= xin[1] "Not enough points at the beginning of the input"
    @assert input_stop <= xin[end] "Not enough points at the end of the input"

    #find relevant input points
    input_start_ind = find_first_above_or_equal(input_start, xin)
    input_stop_ind = find_last_below_or_equal(input_stop, xin)

    input_points = input_stop_ind - input_start_ind + 1
    if input_points < required_input_points
        #upsample
        @assert !ismissing(upsampling_basis) "Not enough input points and no upsampling basis"
        upsample_step = input_width/required_input_points
        input_x = range(input_start + upsample_step/2, step=upsample_step, length=required_input_points)
        in_step = Float64(xin.step)
        input_y = [interpolate_point(xin, yin, xpoint, in_step, in_step, upsampling_basis) for xpoint in input_x] 
    else # don't upsample
        input_x = xin[input_start_ind:input_stop_ind]
        input_y = @view yin[input_start_ind:input_stop_ind]
    end

    right_slice_contribution = .5*slice_weighted_mean(xpoint, slice_width, input_x, input_y, basis)
    
    return left_slice_contribution + right_slice_contribution
end

function prepare_input_or_upsample()
    input_points = input_stop_ind - input_start_ind + 1
    if input_points < required_input_points
        #upsample
        @assert !ismissing(upsampling_basis) "Not enough input points and no upsampling basis"
        upsample_step = input_width/required_input_points
        input_x = range(input_start + upsample_step/2, step=upsample_step, length=required_input_points)
        in_step = Float64(xin.step)
        input_y = [interpolate_point(xin, yin, xpoint, in_step, in_step, HannBasis(3.)) for xpoint in input_x] 
        # TODO change to Lanczos Basis
    else # don't upsample
        input_x = xin[input_start_ind:input_stop_ind]
        input_y = @view yin[input_start_ind:input_stop_ind]
    end
    return (input_x, input_y)
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
    if win_acc <= 0
        error("shouldn't happen")
    end
    return val_acc / win_acc # normalization
end

end # module
