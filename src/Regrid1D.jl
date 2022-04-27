module Regrid1D

include("basis_functions.jl")
include("range_utilities.jl")

export regrid

# TODO change from required_input_points per slice to required_input_points per unit width
function regrid(xin::StepRangeLen, yin, xout,
    smoothing_function::FiniteBasisFunction = RectangularBasis();
    required_input_points=4, upsampling_basis=LanczosBasis()
    )
    # allocate
    yout = Array{Float64, 1}(undef, length(xout))

    # first slice
    out_ind = 1
    right_slice_width = calculate_right_slice_width(xout, out_ind)
    left_slice_width = right_slice_width
    
    while true
        yout[out_ind] = interpolate_point(xin, yin, xout[out_ind], left_slice_width, right_slice_width, 
        smoothing_function, 
        required_input_points=required_input_points, 
        upsampling_basis=upsampling_basis)

        # prepare next point
        out_ind += 1
        out_ind > length(xout) && break

        if out_ind < length(xout) # keep previous slice width on last element
            right_slice_width = calculate_right_slice_width(xout, out_ind)
        end

        left_slice_width = xout[out_ind] - xout[out_ind - 1]
    end

    yout
end

function calculate_right_slice_width(xout, out_ind)
    right_slice_width = xout[out_ind + 1] - xout[out_ind]
    right_slice_width > 0 || throw(ArgumentError(
        "xout must be increasing everywhere. Violated between index $out_ind "*
        "and $(out_ind + 1): $(xout[out_ind]) >= $(xout[out_ind+1])"
    ))
    right_slice_width
end

struct SliceContribution
        slice_width
        input_width
        input_start
        input_stop
        input_start_ind
        input_stop_ind
        upsampling_required
end

function SliceContribution(slice_width, basis::FiniteBasisFunction, xpoint, xin, left::Bool, required_input_points)
    input_width = slice_width*basis.width
    if left
        input_start = xpoint - input_width
        input_stop = xpoint
    else
        input_start = xpoint
        input_stop = xpoint + input_width
    end
    @assert input_start > xin[1] - Float64(xin.step) "Not enough points at the beginning of the input"
    @assert input_stop < xin[end] + Float64(xin.step) "Not enough points at the end of the input"
    # TODO Float64(xin.step) is calculated again and again in different places -> bad?

    #find relevant input indices
    input_start_ind = find_first_above_or_equal(input_start, xin) 
    input_stop_ind = if left find_last_below(input_stop, xin) else find_last_below_or_equal(input_stop, xin) end

    input_points = input_stop_ind - input_start_ind + 1
    upsampling_required = input_points < required_input_points

    SliceContribution(
        slice_width,
        input_width,
        input_start,
        input_stop,
        input_start_ind,
        input_stop_ind,
        upsampling_required,
    )
end

# TODO is it called unit or slice? What's what? Make consistent and document!
function interpolate_point(xin, yin, xpoint, left_unit_width, right_unit_width, basis::FiniteBasisFunction;
    required_input_points=1, upsampling_basis=missing,
    )

    left = SliceContribution(left_unit_width, basis, xpoint, xin, true, required_input_points)
    right = SliceContribution(right_unit_width, basis, xpoint, xin, false, required_input_points)

    if left.upsampling_required || right.upsampling_required
        # upsample both slices, even if only required for one of them
        @assert !ismissing(upsampling_basis) "Interpolation required but no upsampling basis given"
        upsample_step = min(left.input_width, right.input_width)/required_input_points
        left_x, left_y = upsample_prepare_input(left, xin, yin, upsample_step, upsampling_basis)
        right_x, right_y = upsample_prepare_input(right, xin, yin, upsample_step, upsampling_basis)
    else
        # neither of the slices needs upsampling
        left_x, left_y = prepare_input(left, xin, yin)
        right_x, right_y = prepare_input(right, xin, yin)
    end

    left_val_acc, left_win_acc = slice_weighted_mean(xpoint, left_unit_width, left_x, left_y, basis)
    right_val_acc, right_win_acc = slice_weighted_mean(xpoint, right_unit_width, right_x, right_y, basis)

    return (left_val_acc + right_val_acc)/(left_win_acc + right_win_acc)
end

function prepare_input(sc::SliceContribution, xin, yin)
    input_x = xin[sc.input_start_ind:sc.input_stop_ind]
    input_y = @view yin[sc.input_start_ind:sc.input_stop_ind]
    return (input_x, input_y)
end

function upsample_prepare_input(sc::SliceContribution, xin, yin, upsample_step, upsampling_basis)
    input_x = range(sc.input_start + upsample_step/2, step=upsample_step, stop=sc.input_stop)
    in_step = Float64(xin.step)
    input_y = [interpolate_point(xin, yin, up_x, in_step, in_step, upsampling_basis) for up_x in input_x] 
    # TODO use buffer for upsampled values (length == required_input_values) that is allocated at a call to `regrid` instead of allocating every time
    return (input_x, input_y)
end

function slice_weighted_mean(xpoint, slice_width, xin_points, yin_points, basis)
    val_acc::Float64 = 0.
    win_acc::Float64 = 0.
    for i in eachindex(xin_points)
        x = xin_points[i]
        y = yin_points[i]
        rel_pos = abs(x - xpoint)/slice_width
        win_value::Float64 = basis_value(basis, rel_pos)
        win_acc += win_value
        val_acc += y*win_value
    end
    @assert win_acc > 0 "Slice weight should be bigger zero"
    return (val_acc, win_acc)
end

end # module
