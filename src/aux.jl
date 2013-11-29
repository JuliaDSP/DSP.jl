export arraysplit

# Split an array into subarrays of length N, with overlapping regions
# of length M.
function arraysplit(s, n::Integer, m::Integer)
    # n = m is a problem - the algorithm will not terminate.
    if !(0 <= m < n)
        error("m must be between zero and n.")
    end

    # the length of the non-overlapping array stride
    l = n - m
    
    # total number of strides is the total length of the signal divided
    # by the unique number of elements per stride.  extra elements at the
    # end of of the signal are dropped.
    k = int(length(s)/l - n/l + 1)
    [s[(a*l + 1):(a*l + n)] for a=0:(k-1)]
end