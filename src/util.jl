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
        offset = floor((d.+thresh) / (range)) * range
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

function hilbert{T <: FFTW.fftwReal}(x::Array{T})
# Return the Hilbert transform of x (a real signal).
# Code inspired by Scipy's implementation, which is under BSD license.
    N = length(x)
    X = zeros(Complex{T}, N)
    p = FFTW.Plan(x, X, 1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT)
    FFTW.execute(T, p.plan)
    for i = 2:div(N, 2)+isodd(N)
        @inbounds X[i] *= 2.0
    end
    return ifft!(X)
end
hilbert{T <: Real}(x::Array{T}) = hilbert(float(x))

end # end module definition
