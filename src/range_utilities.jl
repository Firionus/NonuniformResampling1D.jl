"""
    find_float_offset(r::AbstractRange, x::Real)

Return the imagined floating point offset, i.e. the index of the range if its
index were to start at 0, such that the value of `r` is `x`.

Indices must be spaced by +1, i.e. valid indices must form a UnitRange.

The output is not restricted to the set of possible offsets within the range. 

The basic underlying math is:  
`y = first + step*(i-firstindex)` for all i in firstindex:1:lastindex,  
therefore:  
`i = (y-first)/step + firstindex`.
"""
function find_float_offset(r::AbstractRange, x::Real)
    (x - first(r))/step(r)
end

"""
    find_first_above_or_equal(cutoff::Real, r::AbstractRange)

Returns the index of the first element where the value is above or equal to the
cutoff. 

Assumes `step(r) > 0`.
"""
function find_first_above_or_equal(cutoff::Real, r::AbstractRange)
    i = (find_float_offset(r, cutoff)|>ceil|>Int) + firstindex(r)
    @assert i <= lastindex(r) # none of the elements are above or equal the cutoff
    max(firstindex(r), i)
end

"""
    find_last_below_or_equal(cutoff::Real, r::AbstractRange)

Returns the index of the last element where the value is below or equal to the
cutoff. 

Assumes `step(r) > 0`.
"""
function find_last_below_or_equal(cutoff::Real, r::AbstractRange)
    i = (find_float_offset(r, cutoff)|>floor|>Int) + firstindex(r)
    @assert i >= firstindex(r) # none of the elements are below or equal the cutoff
    return min(lastindex(r), i)
end

"""
    find_last_below_or_equal(cutoff::Real, r::AbstractRange)

Returns the index of the last element where the value is below the cutoff. 

Assumes `step(r) > 0`.
"""
function find_last_below(cutoff::Real, r::AbstractRange)
    i = (find_float_offset(r, cutoff)|>prevfloat|>floor|>Int) + firstindex(r)
    @assert i >= firstindex(r) # none of the elements are below or equal the cutoff
    return min(lastindex(r), i)
end