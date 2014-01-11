module Util

export unwrap!, unwrap, hilbert

function unwrap!{T <: FloatingPoint}(m::Array{T}, dim::Integer=ndims(m);
                                     range::Number=2pi)
    thresh = range / 2
    if size(m, dim) < 2
        return m
    end
    for i = 2:size(m, dim)
        d = slicedim(m, dim, i) - slicedim(m, dim, i-1)
        slice_tuple = ntuple(ndims(m), n->(n==dim ? (i:i) : (1:size(m,n))))
        offset = floor((d+thresh) / (range)) * range
#        println("offset: ", offset)
#        println("typeof(offset): ", typeof(offset))
#        println("typeof(m[slice_tuple...]): ", typeof(m[slice_tuple...]))
#        println("slice_tuple: ", slice_tuple)
#        println("m[slice_tuple...]: ", m[slice_tuple...])
        m[slice_tuple...] = m[slice_tuple...] - offset
    end
    return m
end

function unwrap{T <: FloatingPoint}(m::Array{T}, args...; kwargs...)
    unwrap!(copy(m), args...; kwargs...)
end

function hilbert{T <: Real}(x::Array{T})
# Return the Hilbert transform of x (a real signal).
# Code inspired by Scipy's implementation, which is under BSD license.
    X = vcat(rfft(x), zeros(int(floor(length(x)/2))-1))
    N = length(X)
    h = zeros(N)
    if N % 2 == 0
        h[1] = h[N/2+1] = 1
        h[2:N/2] = 2
    else
        h[1] = 1
        h[2:(N+1)/2+1] = 2
    end
    return ifft(X.*h)
end

end # end module definition
