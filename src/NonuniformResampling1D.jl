module NonuniformResampling1D

using ArgCheck

include("WindowFunction.jl")
include("window_functions.jl")
include("range_utilities.jl")

export nuresample

"""
    nuresample(xin, yin, xout, smoothing_function; kwargs...)

Take data `yin` sampled at uniformly spaced locations `xin` and resample at
nonuniform locations `xout`. 

# Positional Arguments

- `xin` is an `AbstractRange` of input locations
- `yin` is an `AbstractArray` with the input values
- `xout` is an `AbstractArray` with  the the increasing output locations
- `smoothing_function` is a [`WindowFunction`](@ref). The width of the window
  function is scaled proportional to the distance to the nearest neighbor on the
  left and right side of the output point. At the start and end of `xout` a
  symmetrical window is used. The default is a rectangular window without
  overlap. 
    
# Keyword Arguments

- `upsampling_function`: a [`WindowFunction`](@ref) that is used to perform
  upsampling when the density of input points is not high enough compared to the
  density of output points. The default is Lanczos3. 
- `required_points_per_slice`: the number of input points that must be in a
  slice (left or right half of a window) before smoothing. If the number of
  input points is lower than this value, upsampling is performed to create
  `required_points_per_slice` many points before applying the smoothing
  function.  
  The default is `round(4 * smoothing_function.width)`, which means 4 points are
  required per unit width to the nearest neighbor.  
  The minimum value is 1, in which case upsampling is only performed if no input
  point falls in a slice. 

# Examples

```jldoctest
julia> nuresample(1:9, 1:9, [4.2, 6.2])
2-element Vector{Float64}:
 4.201106036510037
 6.201106036510037
``` 

```jldoctest
julia> nuresample(1:9, 1:9, [4.2, 6.2], 
        rect_window(.5), # moving average without overlap
        required_points_per_slice = 1, # no upsampling
        upsampling_function = lanczos_window(2)) # Lanczos2 interpolation if a slice were empty
2-element Vector{Float64}:
 4.5
 6.5
``` 

# Output Type

Returns an `Array{Float64, 1}`. Other output types are currently unsupported. 
"""
function nuresample(xin::AbstractRange, yin::AbstractArray, xout::AbstractArray,
    smoothing_function::WindowFunction = rect_window();
    required_points_per_slice::Integer=Int(round(4 * smoothing_function.width)), 
    upsampling_function::WindowFunction=lanczos_window()
    )
    # validate inputs
    Base.require_one_based_indexing(yin, xout)
    if step(xin) < 0
        xin, yin = (xin, yin).|>reverse
    end
    @assert step(xin) > 0 "xin should be increasing at this point"
    @argcheck required_points_per_slice >= 1 "required_points_per_slice must at least be 1"
    @argcheck smoothing_function.width > 0 "smoothing_function must have a width bigger than 0. "*
        "Maybe use upsampling_function as smoothing_function."
    @argcheck upsampling_function.width >= 1 "upsampling_function width must at least be 1 to ensure "*
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
    start > first(xin) - step(xin) || "not enough points at the beginning of " +
    "the input for slice from $start to $stop"|>error
    stop < last(xin) + step(xin) || "not enough points at the end of the input " +
    "for slice from $start to $stop"|>error

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
    # TODO is this safe with non-convetional indexing?
    input_x = xin[sc.start_ind:sc.stop_ind]
    input_y = @view yin[sc.start_ind:sc.stop_ind]
    return (input_x, input_y)
end

function upsample_prepare_input(slice::Slice, xin, yin, upsample_step, upsampling_function)
    input_x = range(slice.start + upsample_step/2, step=upsample_step, stop=slice.stop)
    in_step = step(xin)
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
