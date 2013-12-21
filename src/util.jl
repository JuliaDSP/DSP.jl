module Util

export unwrap!, unwrap

function unwrap!(m::Array, dim::Integer=ndims(m); range::Number=2pi)
    if size(m, dim) < 2
        return m
    end
    for i = 2:size(m, dim)
        d = slicedim(m, dim, i) - slicedim(m, dim, i-1)
        slice_tuple = ntuple(ndims(m), n->(n==dim ? (i:i) : (1:size(m,n))))
        offset = floor((d+pi) / (2pi)) * 2pi
#        println("offset: ", offset)
#        println("typeof(offset): ", typeof(offset))
#        println("typeof(m[slice_tuple...]): ", typeof(m[slice_tuple...]))
#        println("slice_tuple: ", slice_tuple)
#        println("m[slice_tuple...]: ", m[slice_tuple...])
        m[slice_tuple...] = m[slice_tuple...] - offset
    end
    return m
end

function unwrap(m::Array, args...; kwargs...)
    unwrap!(copy(m), args...; kwargs...)
end

#function unwrap(m::Array; args...)
#    unwrap!(copy(m); args...)
#end

end # end module definition
