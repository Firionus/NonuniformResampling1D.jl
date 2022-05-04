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

# TODO test
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