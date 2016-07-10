module Util
using Compat
import Compat.view
import Base.Operators: *

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
        shiftin!,
        finddelay(),
        shiftsignals(),
        shiftsignals!(),
        alignsignals(),
        alignsignals!(),
        @julia_newer_than

macro julia_newer_than(version, iftrue, iffalse)
    isa(version, Expr) && version.head === :macrocall && length(version.args) == 2 && version.args[1] === @compat(Symbol("@v_str")) ||
        throw(ArgumentError("invalid syntax"))
    VERSION >= convert(VersionNumber, version.args[2]) ? esc(iftrue) : esc(iffalse)
end

function unwrap!{T <: AbstractFloat}(m::Array{T}, dim::Integer=ndims(m);
                                     range::Number=2pi)
    thresh = range / 2
    if size(m, dim) < 2
        return m
    end
    for i = 2:size(m, dim)
        d = slicedim(m, dim, i:i) - slicedim(m, dim, i-1:i-1)
        slice_tuple = ntuple(n->(n==dim ? (i:i) : (1:size(m,n))), ndims(m))
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

function unwrap{T <: AbstractFloat}(m::Array{T}, args...; kwargs...)
    unwrap!(copy(m), args...; kwargs...)
end

function hilbert{T<:FFTW.fftwReal}(x::StridedVector{T})
# Return the Hilbert transform of x (a real signal).
# Code inspired by Scipy's implementation, which is under BSD license.
    N = length(x)
    X = zeros(Complex{T}, N)
    @julia_newer_than v"0.4.0-dev+6068" begin
        p = plan_rfft(x)
        A_mul_B!(view(X, 1:(N >> 1)+1), p, x)
    end begin
        p = FFTW.Plan(x, X, 1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT)
        FFTW.execute(T, p.plan)
    end
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

    @julia_newer_than v"0.4.0-dev+6068" begin
        p1 = plan_rfft(xc)
        Xsub = view(X, 1:(N >> 1)+1)
        p2 = plan_bfft!(X)
    end begin
        p1 = FFTW.Plan(xc, X, 1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT)
        p2 = FFTW.Plan(X, X, 1, FFTW.BACKWARD, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT)
    end

    normalization = 1/N
    off = 1
    for i = 1:Base.trailingsize(x, 2)
        copy!(xc, 1, x, off, N)

        # fft
        fill!(X, 0)
        @julia_newer_than(v"0.4.0-dev+6068",
                          A_mul_B!(Xsub, p1, xc),
                          FFTW.execute(T, p1.plan))

        # scale real part
        for i = 2:div(N, 2)+isodd(N)
            @inbounds X[i] *= 2.0
        end

        # ifft
        @julia_newer_than(v"0.4.0-dev+6068",
                          A_mul_B!(X, p2, X),
                          FFTW.execute(T, p2.plan))

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
fftintype{T<:Base.FFTW.fftwNumber}(::Type{T}) = T
fftintype{T<:Real}(::Type{T}) = Float64
fftintype{T<:Complex}(::Type{T}) = Complex128

# Get the return element type of FFT for a given type
fftouttype{T<:Base.FFTW.fftwComplex}(::Type{T}) = T
fftouttype{T<:Base.FFTW.fftwReal}(::Type{T}) = Complex{T}
fftouttype{T<:@compat(Union{Real,Complex})}(::Type{T}) = Complex128

# Get the real part of the return element type of FFT for a given type
fftabs2type{T<:Base.FFTW.fftwReal}(::Type{Complex{T}}) = T
fftabs2type{T<:Base.FFTW.fftwReal}(::Type{T}) = T
fftabs2type{T<:@compat(Union{Real,Complex})}(::Type{T}) = Float64

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
Base.step(x::Frequencies) = x.multiplier

fftfreq(n::Int, fs::Real=1) = Frequencies(((n-1) >> 1)+1, n, fs/n)
rfftfreq(n::Int, fs::Real=1) = Frequencies((n >> 1)+1, (n >> 1)+1, fs/n)
Base.fftshift(x::Frequencies) = (x.nreal-x.n:x.nreal-1)*x.multiplier

# Get next fast FFT size for a given signal length
const FAST_FFT_SIZES = [2, 3, 5, 7]
nextfastfft(n) = nextprod(FAST_FFT_SIZES, n)
nextfastfft(n1, n2...) = tuple(nextfastfft(n1), nextfastfft(n2...)...)
nextfastfft(n::Tuple) = nextfastfft(n...)


## COMMON DSP TOOLS

immutable dBconvert end
immutable dBaconvert end
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
rms{T<:Number}(s::AbstractArray{T}) = sqrt(sumabs2(s)/length(s))
# root mean square of fft of signal
rmsfft{T<:Complex}(f::AbstractArray{T}) = sqrt(sumabs2(f))/length(f)



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

@inline function unsafe_dot{T<:Base.LinAlg.BlasReal}(a::Matrix{T}, aColIdx::Integer, b::Vector{T}, bLastIdx::Integer)
    BLAS.dot(size(a, 1), pointer(a, size(a, 1)*(aColIdx-1) + 1), 1, pointer(b, bLastIdx - size(a, 1) + 1), 1)
end

function unsafe_dot{T}(a::AbstractMatrix, aColIdx::Integer, b::AbstractVector{T}, c::AbstractVector{T}, cLastIdx::Integer)
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

function unsafe_dot{T}(a::T, b::AbstractArray, bLastIdx::Integer)
    aLen     = length(a)
    bBaseIdx = bLastIdx - aLen
    @inbounds dotprod  = a[1] * b[bBaseIdx + 1]
    @simd for i in 2:aLen
        @inbounds dotprod += a[i] * b[bBaseIdx + i]
    end

    return dotprod
end

@inline function unsafe_dot{T<:Base.LinAlg.BlasReal}(a::Vector{T}, b::Array{T}, bLastIdx::Integer)
    BLAS.dot(length(a), pointer(a), 1, pointer(b, bLastIdx - length(a) + 1), 1)
end

function unsafe_dot{T}(a::AbstractVector, b::AbstractVector{T}, c::AbstractVector{T}, cLastIdx::Integer)
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
function shiftin!{T}(a::Vector{T}, b::Vector{T})
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



# Delay Utility Functions:

# finddelay(): use peak of cross correlation from origin of time to calculate
# time lag between two Arrays (signals)

function finddelay{T <: Real}(x::AbstractArray{T, 1}, u::AbstractArray{T, 1})

sₓᵤ = xcorr(x, u)

ct_idx = cld(length(sₓᵤ), 2)

_, pk_idx = findmax(sₓᵤ, 1)

δ = ct_idx - pk_idx[1]

return δ

end

# shiftsignals(): shift elements of an Array (signal) of a given amount of
# samples

function shiftsignals{T <: Real}(u::AbstractArray{T, 1}, δ::Int)

lᵤ = length(u)

y = zeros(T, lᵤ)

if δ > 0
    y[1:(lᵤ - δ)] = u[(δ + 1):lᵤ]
else
    y[(-δ + 1):lᵤ] = u[1:(lᵤ - -δ)]
end

return y

end

function shiftsignals!{T <: Real}(u::AbstractArray{T, 1}, δ::Int)

lᵤ = length(u)

if δ > 0

    deleteat!(u, 1:δ)
    
    # append!() could be used, but this is faster and prevents allocation.
    for d = 1:δ
        insert!(u, lᵤ - δ + d, 0)
    end
    
elseif δ < 0

    deleteat!(u, (lᵤ - -δ + 1):lᵤ)
    
    # prepend!() could be used, but this is faster and prevents allocation.
    for d = 1:(-δ)
        insert!(u, d, 0)
    end

end

end

# alignsignals(): use finddelay() and shiftsignals() to time align Arrays
# (signals)

function alignsignals{T <: Real}(x::AbstractArray{T, 1}, u::AbstractArray{T, 1})

δ = finddelay(x, u)

y = shiftsignals(u, δ)

return y, δ

end

function alignsignals!{T <: Real}(x::AbstractArray{T, 1}, u::AbstractArray{T, 1})

δ = finddelay(x, u)

shiftsignals!(u, δ)

return δ

end

end # end module definition
