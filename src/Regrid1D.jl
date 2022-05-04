module Regrid1D

include("WindowFunction.jl")
include("window_functions.jl")
include("range_utilities.jl")

export regrid

function regrid(xin::AbstractRange, yin, xout, smoothing_function=rect_window(); kwargs...)
    regrid(StepRangeLen(xin), yin, xout, smoothing_function; kwargs...)
end

function regrid(xin::StepRangeLen, yin, xout,
    smoothing_function::WindowFunction = rect_window();
    required_points_per_slice::Integer=Int(round(4 * smoothing_function.width)), 
    upsampling_function::WindowFunction=lanczos_window()
    )
    # validate inputs
    @assert required_points_per_slice >= 1 "required_points_per_slice must at least be 1"
    @assert smoothing_function.width > 0 "smoothing_function must have a width bigger than 0. "*
        "Maybe use upsampling_function as smoothing_function."
    @assert upsampling_function.width >= 1 "upsampling_function width must at least be 1 to ensure "*
        "there's at least one input point to each side of the upsampled point"
    # allocate
    yout = Array{Float64, 1}(undef, length(xout))

    # first unit
    out_ind = 1
    right_unit_width = calculate_right_unit_width(xout, out_ind)
    left_unit_width = right_unit_width
    
    while true
        yout[out_ind] = interpolate_point(xin, yin, xout[out_ind], left_unit_width, right_unit_width, 
        smoothing_function, 
        required_points_per_slice=required_points_per_slice, 
        upsampling_function=upsampling_function)

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

function Slice(unit_width, window::WindowFunction, xpoint, xin, left::Bool, required_points_per_slice)
    # the unit width is the width to the neighboring point  
    # now calculate slice width, which is the width in which points will be
    # considered according to the window funciton width
    width = unit_width*window.width
    if left
        start = xpoint - width
        stop = xpoint
    else
        start = xpoint
        stop = xpoint + width
    end
    @assert start > xin[1] - Float64(xin.step) "not enough points at the beginning of the input for slice from $start to $stop"
    @assert stop < xin[end] + Float64(xin.step) "not enough points at the end of the input for slice from $start to $stop"

    #find relevant input indices
    start_ind = find_first_above_or_equal(start, xin) 
    stop_ind = if left find_last_below(stop, xin) else find_last_below_or_equal(stop, xin) end

    input_points = stop_ind - start_ind + 1
    upsampling_required = input_points < required_points_per_slice

    Slice(
        width,
        start,
        stop,
        start_ind,
        stop_ind,
        upsampling_required,
    )
end

function interpolate_point(xin, yin, xpoint, left_unit_width, right_unit_width, window::WindowFunction;
    required_points_per_slice=1, upsampling_function=missing)

    left = Slice(left_unit_width, window, xpoint, xin, true, required_points_per_slice)
    right = Slice(right_unit_width, window, xpoint, xin, false, required_points_per_slice)

    if left.upsampling_required || right.upsampling_required
        # upsample both slices, even if only required for one of them, to keep balance between left and right
        # while still maintaining continuity when an output point is moved slightly over an input point
        @assert !ismissing(upsampling_function) "interpolation required but no upsampling_function given"
        upsample_step = min(left.width, right.width)/required_points_per_slice
        left_x, left_y = upsample_prepare_input(left, xin, yin, upsample_step, upsampling_function)
        right_x, right_y = upsample_prepare_input(right, xin, yin, upsample_step, upsampling_function)
    else
        # neither of the slices needs upsampling
        left_x, left_y = prepare_input(left, xin, yin)
        right_x, right_y = prepare_input(right, xin, yin)
    end

    left_val_acc, left_win_acc = slice_weighted_mean(xpoint, left_unit_width, left_x, left_y, window)
    right_val_acc, right_win_acc = slice_weighted_mean(xpoint, right_unit_width, right_x, right_y, window)

    return (left_val_acc + right_val_acc)/(left_win_acc + right_win_acc)
end

function prepare_input(sc::Slice, xin, yin)
    input_x = xin[sc.start_ind:sc.stop_ind]
    input_y = @view yin[sc.start_ind:sc.stop_ind]
    return (input_x, input_y)
end

function upsample_prepare_input(slice::Slice, xin, yin, upsample_step, upsampling_function)
    input_x = range(slice.start + upsample_step/2, step=upsample_step, stop=slice.stop)
    in_step = Float64(xin.step)
    input_y = [interpolate_point(xin, yin, up_x, in_step, in_step, upsampling_function) for up_x in input_x] 
    return (input_x, input_y)
end

function slice_weighted_mean(xpoint, unit_width, xin_points, yin_points, window)
    val_acc::Float64 = 0.
    win_acc::Float64 = 0.
    for i in eachindex(xin_points)
        x = xin_points[i]
        y = yin_points[i]
        rel_pos = abs(x - xpoint)/unit_width
        win_value::Float64 = window(rel_pos)
        win_acc += win_value
        val_acc += y*win_value
    end
    @assert win_acc > 0 "Slice weight should be bigger zero"
    return (val_acc, win_acc)
end

end # module
