module Util

export unwrap!, unwrap, hilbert, Frequencies, fftintype, fftouttype,
       fftabs2type, fftfreq, rfftfreq, nextfastfft

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

function hilbert{T<:FFTW.fftwReal}(x::StridedVector{T})
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
hilbert{T<:Real}(x::AbstractVector{T}) = hilbert(convert(Vector{fftintype(T)}, x))

function hilbert{T<:Real}(x::AbstractArray{T})
    N = size(x, 1)
    xc = Array(fftintype(T), N)
    X = Array(fftouttype(T), N)
    out = similar(x, fftouttype(T))

    p1 = FFTW.Plan(xc, X, 1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT)
    p2 = FFTW.Plan(X, X, 1, FFTW.BACKWARD, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT)

    normalization = 1/N
    off = 1
    for i = 1:Base.trailingsize(x, 2)
        copy!(xc, 1, x, off, N)

        # fft
        fill!(X, 0)
        FFTW.execute(T, p1.plan)

        # scale real part
        for i = 2:div(N, 2)+isodd(N)
            @inbounds X[i] *= 2.0
        end

        # ifft
        FFTW.execute(T, p2.plan)

        # scale and copy to output
        for j = 1:N
            @inbounds out[off+j-1] = X[j]*normalization
        end

        off += N
    end

    out
end

## FFT TYPES

# Get the input element type of FFT for a given type
fftintype{T<:Base.FFTW.fftwNumber}(::Type{T}) = T
fftintype{T<:Real}(::Type{T}) = Float64
fftintype{T<:Complex}(::Type{T}) = Complex128

# Get the return element type of FFT for a given type
fftouttype{T<:Base.FFTW.fftwComplex}(::Type{T}) = T
fftouttype{T<:Base.FFTW.fftwReal}(::Type{T}) = Complex{T}
fftouttype{T<:Union(Real,Complex)}(::Type{T}) = Complex128

# Get the real part of the return element type of FFT for a given type
fftabs2type{T<:Base.FFTW.fftwReal}(::Type{Complex{T}}) = T
fftabs2type{T<:Base.FFTW.fftwReal}(::Type{T}) = T
fftabs2type{T<:Union(Real,Complex)}(::Type{T}) = Float64

## FREQUENCY VECTOR

immutable Frequencies <: AbstractVector{Float64}
    nreal::Int
    n::Int
    multiplier::Float64
end

unsafe_getindex(x::Frequencies, i::Int) =
    (i-1+ifelse(i <= x.nreal, 0, -x.n))*x.multiplier
function Base.getindex(x::Frequencies, i::Int)
    (i >= 1 && i <= x.n) || throw(BoundsError())
    unsafe_getindex(x, i)
end
Base.start(x::Frequencies) = 1
Base.next(x::Frequencies, i::Int) = (unsafe_getindex(x, i), i+1)
Base.done(x::Frequencies, i::Int) = i > x.n
Base.size(x::Frequencies) = (x.n,)
Base.similar(x::Frequencies, T::Type, args...) = Array(T, args...)
Base.step(x::Frequencies) = x.multiplier

fftfreq(n::Int, fs::Real=1) = Frequencies(((n-1) >> 1)+1, n, fs/n)
rfftfreq(n::Int, fs::Real=1) = Frequencies((n >> 1)+1, (n >> 1)+1, fs/n)
Base.fftshift(x::Frequencies) = (x.nreal-x.n:x.nreal-1)*x.multiplier

# Get next fast FFT size for a given signal length
const FAST_FFT_SIZES = [2, 3, 5, 7]
nextfastfft(n) = nextprod(FAST_FFT_SIZES, n)

end # end module definition
