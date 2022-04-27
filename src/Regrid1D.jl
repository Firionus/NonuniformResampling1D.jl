module Regrid1D

include("basis_functions.jl")
include("range_utilities.jl")

export regrid

# TODO change from required_input_points per slice to required_input_points per unit width to normalize
# the parameter against different width basis function
function regrid(xin::StepRangeLen, yin, xout,
    smoothing_function::FiniteBasisFunction = RectangularBasis();
    required_input_points=4, upsampling_basis=LanczosBasis()
    )
    # allocate
    yout = Array{Float64, 1}(undef, length(xout))

    # first unit
    out_ind = 1
    right_unit_width = calculate_right_unit_width(xout, out_ind)
    left_unit_width = right_unit_width
    
    while true
        yout[out_ind] = interpolate_point(xin, yin, xout[out_ind], left_unit_width, right_unit_width, 
        smoothing_function, 
        required_input_points=required_input_points, 
        upsampling_basis=upsampling_basis)

        # prepare next point
        left_unit_width = right_unit_width # next left unit width is the last right unit width
        out_ind += 1
        out_ind > length(xout) && break

        if out_ind < length(xout) # keep previous slice width on last element
            right_unit_width = calculate_right_unit_width(xout, out_ind)
        end
    end

    yout
end

function calculate_right_unit_width(xout, out_ind)
    right_unit_width = xout[out_ind + 1] - xout[out_ind]
    right_unit_width > 0 || throw(ArgumentError(
        "xout must be increasing everywhere. Violated between index $out_ind "*
        "and $(out_ind + 1): $(xout[out_ind]) >= $(xout[out_ind+1])"
    ))
    right_unit_width
end

struct Slice
        width
        start
        stop
        start_ind
        stop_ind
        upsampling_required
end

function Slice(unit_width, basis::FiniteBasisFunction, xpoint, xin, left::Bool, required_input_points)
    # the unit width is the width to the neighboring point
    # now calculate slice width, which is the width in which points will be considered according to the finite basis
    width = unit_width*basis.width
    if left
        start = xpoint - width
        stop = xpoint
    else
        start = xpoint
        stop = xpoint + width
    end
    @assert start > xin[1] - Float64(xin.step) "Not enough points at the beginning of the input"
    @assert stop < xin[end] + Float64(xin.step) "Not enough points at the end of the input"
    # TODO Float64(xin.step) is calculated again and again in different places -> bad?

    #find relevant input indices
    start_ind = find_first_above_or_equal(start, xin) 
    stop_ind = if left find_last_below(stop, xin) else find_last_below_or_equal(stop, xin) end

    input_points = stop_ind - start_ind + 1
    upsampling_required = input_points < required_input_points

    Slice(
        width,
        start,
        stop,
        start_ind,
        stop_ind,
        upsampling_required,
    )
end

function interpolate_point(xin, yin, xpoint, left_unit_width, right_unit_width, basis::FiniteBasisFunction;
    required_input_points=1, upsampling_basis=missing,
    )

    left = Slice(left_unit_width, basis, xpoint, xin, true, required_input_points)
    right = Slice(right_unit_width, basis, xpoint, xin, false, required_input_points)

    if left.upsampling_required || right.upsampling_required
        # upsample both slices, even if only required for one of them, to keep balance between left and right
        # while still maintaining continuity when an output point is moved slightly over an input point
        @assert !ismissing(upsampling_basis) "interpolation required but no upsampling basis given"
        upsample_step = min(left.width, right.width)/required_input_points
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

function prepare_input(sc::Slice, xin, yin)
    input_x = xin[sc.start_ind:sc.stop_ind]
    input_y = @view yin[sc.start_ind:sc.stop_ind]
    return (input_x, input_y)
end

function upsample_prepare_input(slice::Slice, xin, yin, upsample_step, upsampling_basis)
    input_x = range(slice.start + upsample_step/2, step=upsample_step, stop=slice.stop)
    in_step = Float64(xin.step)
    input_y = [interpolate_point(xin, yin, up_x, in_step, in_step, upsampling_basis) for up_x in input_x] 
    # TODO use buffer for upsampled values (length == required_input_values) that is allocated at a call to `regrid` instead of allocating every time
    return (input_x, input_y)
end

function slice_weighted_mean(xpoint, unit_width, xin_points, yin_points, basis)
    val_acc::Float64 = 0.
    win_acc::Float64 = 0.
    for i in eachindex(xin_points)
        x = xin_points[i]
        y = yin_points[i]
        rel_pos = abs(x - xpoint)/unit_width
        win_value::Float64 = basis_value(basis, rel_pos)
        win_acc += win_value
        val_acc += y*win_value
    end
    @assert win_acc > 0 "Slice weight should be bigger zero"
    return (val_acc, win_acc)
end

end # module
