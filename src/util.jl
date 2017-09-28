module Util
using ..DSP: @importffts
import Base: *
@importffts

export  unwrap!,
        unwrap,
        hilbert,
        Frequencies,
        fftintype,
        fftouttype,
        fftabs2type,
        fftfreq,
        rfftfreq,
        nextfastfft,
        dB,
        dBa,
        pow2db,
        amp2db,
        db2pow,
        db2amp,
        rms,
        rmsfft,
        unsafe_dot,
        polyfit,
        shiftin!

function unwrap!(m::Array{T}, dim::Integer=ndims(m); range::Number=2pi) where T<:AbstractFloat
    thresh = range / 2
    if size(m, dim) < 2
        return m
    end
    for i = 2:size(m, dim)
        d = slicedim(m, dim, i:i) - slicedim(m, dim, i-1:i-1)
        slice_tuple = ntuple(n->(n==dim ? (i:i) : (1:size(m,n))), ndims(m))
        offset = floor.((d.+thresh) / (range)) * range
#        println("offset: ", offset)
#        println("typeof(offset): ", typeof(offset))
#        println("typeof(m[slice_tuple...]): ", typeof(m[slice_tuple...]))
#        println("slice_tuple: ", slice_tuple)
#        println("m[slice_tuple...]: ", m[slice_tuple...])
        m[slice_tuple...] = m[slice_tuple...] - offset
    end
    return m
end

function unwrap(m::Array{T}, args...; kwargs...) where T<:AbstractFloat
    unwrap!(copy(m), args...; kwargs...)
end

function hilbert(x::StridedVector{T}) where T<:FFTW.fftwReal
# Return the Hilbert transform of x (a real signal).
# Code inspired by Scipy's implementation, which is under BSD license.
    N = length(x)
    X = zeros(Complex{T}, N)
    p = plan_rfft(x)
    A_mul_B!(view(X, 1:(N >> 1)+1), p, x)
    for i = 2:div(N, 2)+isodd(N)
        @inbounds X[i] *= 2.0
    end
    return ifft!(X)
end
hilbert(x::AbstractVector{T}) where {T<:Real} = hilbert(convert(Vector{fftintype(T)}, x))

function hilbert(x::AbstractArray{T}) where T<:Real
    N = size(x, 1)
    xc = Array{fftintype(T)}(N)
    X = Array{fftouttype(T)}(N)
    out = similar(x, fftouttype(T))

    p1 = plan_rfft(xc)
    Xsub = view(X, 1:(N >> 1)+1)
    p2 = plan_bfft!(X)

    normalization = 1/N
    off = 1
    for i = 1:Base.trailingsize(x, 2)
        copy!(xc, 1, x, off, N)

        # fft
        fill!(X, 0)
        A_mul_B!(Xsub, p1, xc)

        # scale real part
        for i = 2:div(N, 2)+isodd(N)
            @inbounds X[i] *= 2.0
        end

        # ifft
        A_mul_B!(X, p2, X)

        # scale and copy to output
        @simd for j = 1:N
            @inbounds out[off+j-1] = X[j]*normalization
        end

        off += N
    end

    out
end

## FFT TYPES

# Get the input element type of FFT for a given type
fftintype(::Type{T}) where {T<:FFTW.fftwNumber} = T
fftintype(::Type{T}) where {T<:Real} = Float64
fftintype(::Type{T}) where {T<:Complex} = Complex128

# Get the return element type of FFT for a given type
fftouttype(::Type{T}) where {T<:FFTW.fftwComplex} = T
fftouttype(::Type{T}) where {T<:FFTW.fftwReal} = Complex{T}
fftouttype(::Type{T}) where {T<:Union{Real,Complex}} = Complex128

# Get the real part of the return element type of FFT for a given type
fftabs2type(::Type{Complex{T}}) where {T<:FFTW.fftwReal} = T
fftabs2type(::Type{T}) where {T<:FFTW.fftwReal} = T
fftabs2type(::Type{T}) where {T<:Union{Real,Complex}} = Float64

## FREQUENCY VECTOR

struct Frequencies <: AbstractVector{Float64}
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
Base.step(x::Frequencies) = x.multiplier

fftfreq(n::Int, fs::Real=1) = Frequencies(((n-1) >> 1)+1, n, fs/n)
rfftfreq(n::Int, fs::Real=1) = Frequencies((n >> 1)+1, (n >> 1)+1, fs/n)
fftshift(x::Frequencies) = (x.nreal-x.n:x.nreal-1)*x.multiplier

# Get next fast FFT size for a given signal length
const FAST_FFT_SIZES = [2, 3, 5, 7]
nextfastfft(n) = nextprod(FAST_FFT_SIZES, n)
nextfastfft(n1, n2...) = tuple(nextfastfft(n1), nextfastfft(n2...)...)
nextfastfft(n::Tuple) = nextfastfft(n...)


## COMMON DSP TOOLS

struct dBconvert end
struct dBaconvert end
const dB = dBconvert()
const dBa = dBaconvert()
# for using e.g. 3dB or -3dBa
*(a::Real, ::dBconvert) = db2pow(a)
*(a::Real, ::dBaconvert) = db2amp(a)
# convert dB to power ratio
db2pow(a::Real) = 10^(a/10)
# convert dB to amplitude ratio
db2amp(a::Real) = 10^(a/20)
# convert power ratio to dB
pow2db(a::Real) = 10*log10(a)
# convert amplitude ratio to dB
amp2db(a::Real) = 20*log10(a)

# root mean square
rms(s::AbstractArray{T}) where {T<:Number} = sqrt(sum(abs2, s)/length(s))
# root mean square of fft of signal
rmsfft(f::AbstractArray{T}) where {T<:Complex} = sqrt(sum(abs2, f))/length(f)



# Computes the dot product of a single column of a, specified by aColumnIdx, with the vector b.
# The number of elements used in the dot product determined by the size(A)[1].
# Note: bIdx is the last element of b used in the dot product.
function unsafe_dot(a::AbstractMatrix, aColIdx::Integer, b::AbstractVector, bLastIdx::Integer)
    aLen     = size(a, 1)
    bBaseIdx = bLastIdx - aLen
    dotprod  = a[1, aColIdx] * b[ bBaseIdx + 1]
    @simd for i in 2:aLen
        @inbounds dotprod += a[i, aColIdx] * b[bBaseIdx + i]
    end

    return dotprod
end

@inline function unsafe_dot(a::Matrix{T}, aColIdx::Integer, b::Vector{T}, bLastIdx::Integer) where T<:Base.LinAlg.BlasReal
    BLAS.dot(size(a, 1), pointer(a, size(a, 1)*(aColIdx-1) + 1), 1, pointer(b, bLastIdx - size(a, 1) + 1), 1)
end

function unsafe_dot(a::AbstractMatrix, aColIdx::Integer, b::AbstractVector{T}, c::AbstractVector{T}, cLastIdx::Integer) where T
    aLen = size(a, 1)
    bLen = length(b)
    bLen == aLen-1  || error( "length(b) must equal to length(a)[1] - 1" )
    cLastIdx < aLen || error( "cLastIdx but be < length(a)")

    dotprod = a[1, aColIdx] * b[cLastIdx]
    @simd for i in 2:aLen-cLastIdx
        @inbounds dotprod += a[i, aColIdx] * b[i+cLastIdx-1]
    end
    @simd for i in 1:cLastIdx
        @inbounds dotprod += a[aLen-cLastIdx+i, aColIdx] * c[i]
    end

    return dotprod
end

function unsafe_dot(a::T, b::AbstractArray, bLastIdx::Integer) where T
    aLen     = length(a)
    bBaseIdx = bLastIdx - aLen
    @inbounds dotprod  = a[1] * b[bBaseIdx + 1]
    @simd for i in 2:aLen
        @inbounds dotprod += a[i] * b[bBaseIdx + i]
    end

    return dotprod
end

@inline function unsafe_dot(a::Vector{T}, b::Array{T}, bLastIdx::Integer) where T<:Base.LinAlg.BlasReal
    BLAS.dot(length(a), pointer(a), 1, pointer(b, bLastIdx - length(a) + 1), 1)
end

function unsafe_dot(a::AbstractVector, b::AbstractVector{T}, c::AbstractVector{T}, cLastIdx::Integer) where T
    aLen    = length(a)
    dotprod = zero(a[1]*b[1])
    @simd for i in 1:aLen-cLastIdx
        @inbounds dotprod += a[i] * b[i+cLastIdx-1]
    end
    @simd for i in 1:cLastIdx
        @inbounds dotprod += a[aLen-cLastIdx+i] * c[i]
    end

    return dotprod
end



# Shifts b into the end a.
# julia> DSP.Util.shiftin!( [1,2,3,4], [5, 6])
# 4-element Array{Int64,1}:
#  3
#  4
#  5
#  6
function shiftin!(a::AbstractVector{T}, b::AbstractVector{T}) where T
    aLen = length(a)
    bLen = length(b)

    if bLen >= aLen
        copy!(a, 1, b, bLen - aLen + 1, aLen)
    else

        for i in 1:aLen-bLen
            @inbounds a[i] = a[i+bLen]
        end
        bIdx = 1
        for i in aLen-bLen+1:aLen
            @inbounds a[i] = b[bIdx]
            bIdx += 1
        end
    end

    return a
end


end # end module definition
