# Implementations of filt, filtfilt, and fftfilt

#
# filt and filt!
#


## PolynomialRatio
_zerosi(f::PolynomialRatio{:z,T}, ::AbstractArray{S}) where {T,S} =
    zeros(promote_type(T, S), max(-firstindex(f.a), -firstindex(f.b)))

"""
    filt!(out, f, x[, si])

Same as [`filt()`](@ref) but writes the result into the `out`
argument. Output array `out` may not be an alias of `x`, i.e. filtering may
not be done in place.
"""
filt!(out, f::PolynomialRatio{:z}, x::AbstractArray, si=_zerosi(f, x)) =
    filt!(out, coefb(f), coefa(f), x, si)

"""
    filt(f, x[, si])

Apply filter or filter coefficients `f` along the first dimension
of array `x`. If `f` is a filter coefficient object, `si`
is an optional array representing the initial filter state (defaults
to zeros). If `f` is a `PolynomialRatio`, `Biquad`, or
`SecondOrderSections`, filtering is implemented directly. If
`f` is a `ZeroPoleGain` object, it is first converted to a
`SecondOrderSections` object.  If `f` is a Vector, it is
interpreted as an FIR filter, and a naïve or FFT-based algorithm is
selected based on the data and filter length.
"""
filt(f::PolynomialRatio{:z}, x, si=_zerosi(f, x)) = filt(coefb(f), coefa(f), x, si)

## SecondOrderSections
_zerosi(f::SecondOrderSections{:z,T,G}, ::AbstractArray{S}) where {T,G,S} =
    zeros(promote_type(T, G, S), 2, length(f.biquads))

# filt! algorithm (no checking, returns si)
function _filt!(out::AbstractArray, si::AbstractArray{S,N}, f::SecondOrderSections{:z},
                x::AbstractArray, col::Int) where {S,N}
    g = f.g
    biquads = f.biquads
    n = length(biquads)
    @inbounds for i in axes(x, 1)
        yi = x[i, col]
        for fi = 1:n
            biquad = biquads[fi]
            xi = yi
            yi = muladd(biquad.b0, xi, si[1, fi])
            si[1, fi] = muladd(biquad.a1, -yi, muladd(biquad.b1, xi, si[2, fi]))
            si[2, fi] = muladd(biquad.b2, xi, -biquad.a2 * yi)
        end
        out[i, col] = yi*g
    end
    si
end

function filt!(out::AbstractArray, f::SecondOrderSections{:z}, x::AbstractArray,
                    si::AbstractArray{S,N}=_zerosi(f, x)) where {S,N}
    biquads = f.biquads
    ncols = Base.trailingsize(x, 2)

    size(x) != size(out) && error("out size must match x")
    (size(si, 1) != 2 || size(si, 2) != length(biquads) || (N > 2 && Base.trailingsize(si, 3) != ncols)) &&
        error("si must be 2 x nbiquads or 2 x nbiquads x nsignals")

    initial_si = si
    g = f.g
    n = length(biquads)
    for col = 1:ncols
        _filt!(out, initial_si[:, :, N > 2 ? col : 1], f, x, col)
    end
    out
end

filt(f::SecondOrderSections{:z,T,G}, x::AbstractArray{S}, si=_zerosi(f, x)) where {T,G,S<:Number} =
    filt!(Array{promote_type(T, G, S)}(undef, size(x)), f, x, si)

## Biquad
_zerosi(::Biquad{:z,T}, ::AbstractArray{S}) where {T,S} =
    zeros(promote_type(T, S), 2)

# filt! algorithm (no checking, returns si)
function _filt!(out::AbstractArray, si1::Number, si2::Number, f::Biquad{:z},
                x::AbstractArray, col::Int)
    @inbounds for i in axes(x, 1)
        xi = x[i, col]
        yi = muladd(f.b0, xi, si1)
        si1 = muladd(f.a1, -yi, muladd(f.b1, xi, si2))
        si2 = muladd(f.b2, xi, -f.a2 * yi)
        out[i, col] = yi
    end
    (si1, si2)
end

# filt! variant that preserves si
function filt!(out::AbstractArray, f::Biquad{:z}, x::AbstractArray,
                    si::AbstractArray{S,N}=_zerosi(f, x)) where {S,N}
    ncols = Base.trailingsize(x, 2)

    size(x) != size(out) && error("out size must match x")
    (size(si, 1) != 2 || (N > 1 && Base.trailingsize(si, 2) != ncols)) &&
        error("si must have two rows and 1 or nsignals columns")

    for col = 1:ncols
        _filt!(out, si[1, N > 1 ? col : 1], si[2, N > 1 ? col : 1], f, x, col)
    end
    out
end

filt(f::Biquad{:z,T}, x::AbstractArray{S}, si=_zerosi(f, x)) where {T,S<:Number} =
    filt!(Array{promote_type(T, S)}(undef, size(x)), f, x, si)

## For arbitrary filters, convert to SecondOrderSections
filt(f::FilterCoefficients{:z}, x) = filt(convert(SecondOrderSections, f), x)
filt!(out, f::FilterCoefficients{:z}, x) = filt!(out, convert(SecondOrderSections, f), x)

"""
    DF2TFilter(coef[, si])

Construct a stateful direct form II transposed filter with
coefficients `coef`. `si` is an optional array representing the
initial filter state (defaults to zeros). If `f` is a
`PolynomialRatio`, `Biquad`, or `SecondOrderSections`,
filtering is implemented directly. If `f` is a `ZeroPoleGain`
object, it is first converted to a `SecondOrderSections` object.
"""
struct DF2TFilter{T<:FilterCoefficients{:z},S<:Array}
    coef::T
    state::S

    function DF2TFilter{Ti,Si}(coef::PolynomialRatio{:z}, state::Vector) where {Ti,Si}
        length(state) == max(length(coefa(coef)), length(coefb(coef)))-1 ||
            throw(ArgumentError("length of state vector must match filter order"))
        new{Ti,Si}(coef, state)
    end
    function DF2TFilter{Ti,Si}(coef::SecondOrderSections{:z}, state::Matrix) where {Ti,Si}
        (size(state, 1) == 2 && size(state, 2) == length(coef.biquads)) ||
            throw(ArgumentError("state must be 2 x nbiquads"))
        new{Ti,Si}(coef, state)
    end
    function DF2TFilter{Ti,Si}(coef::Biquad{:z}, state::Vector) where {Ti,Si}
        length(state) == 2 || throw(ArgumentError("length of state must be 2"))
        new{Ti,Si}(coef, state)
    end
end

## PolynomialRatio
DF2TFilter(coef::PolynomialRatio{:z,T},
           state::Vector{S}=zeros(T, max(length(coefa(coef)), length(coefb(coef)))-1)) where {T,S} =
    DF2TFilter{PolynomialRatio{:z,T}, Vector{S}}(coef, state)

function filt!(out::AbstractVector, f::DF2TFilter{<:PolynomialRatio,<:Vector}, x::AbstractVector)
    length(x) != length(out) && throw(ArgumentError("out size must match x"))

    si = f.state
    n = length(si) + 1
    b = coefb(f.coef)
    if length(b) < n
        append!(b, zeros(n-length(b)))
    end
    a = coefa(f.coef)
    if length(a) < n
        append!(a, zeros(n-length(a)))
    end
    if n == 1
        mul!(out, x, b[1])
    else
        @inbounds for i in eachindex(x, out)
            xi = x[i]
            val = muladd(b[1], xi, si[1])
            for j=2:n-1
                si[j-1] = muladd(a[j], -val, muladd(b[j], xi, si[j]))
            end
            si[n-1] = muladd(b[n], xi, -a[n] * val)
            out[i] = val
        end
    end
    out
end

## SecondOrderSections
DF2TFilter(coef::SecondOrderSections{:z,T,G},
           state::Matrix{S}=zeros(promote_type(T, G), 2, length(coef.biquads))) where {T,G,S} =
    DF2TFilter{SecondOrderSections{:z,T,G}, Matrix{S}}(coef, state)

function filt!(out::AbstractVector, f::DF2TFilter{<:SecondOrderSections,<:Matrix}, x::AbstractVector)
    length(x) != length(out) && throw(ArgumentError("out size must match x"))
    _filt!(out, f.state, f.coef, x, 1)
    out
end

## Biquad
DF2TFilter(coef::Biquad{:z,T}, state::Vector{S}=zeros(T, 2)) where {T,S} =
    DF2TFilter{Biquad{:z,T}, Vector{S}}(coef, state)
function filt!(out::AbstractVector, f::DF2TFilter{<:Biquad,<:Vector}, x::AbstractVector)
    length(x) != length(out) && throw(ArgumentError("out size must match x"))
    si = f.state
    (si[1], si[2]) = _filt!(out, si[1], si[2], f.coef, x, 1)
    out
end

# Variant that allocates the output
filt(f::DF2TFilter{<:FilterCoefficients{:z},<:Array{T}}, x::AbstractVector) where {T} =
    filt!(Vector{T}(undef, length(x)), f, x)

# Fall back to SecondOrderSections
DF2TFilter(coef::FilterCoefficients{:z}) = DF2TFilter(convert(SecondOrderSections, coef))

#
# filtfilt
#

# Extrapolate the beginning of a signal for use by filtfilt. This
# computes:
#
# [(2 * x[1]) .- x[pad_length+1:-1:2],
#  x,
#  (2 * x[end]) .- x[end-1:-1:end-pad_length]]
#
# in place in output. The istart and n parameters determine the portion
# of the input signal x to extrapolate.
function extrapolate_signal!(out, ostart, sig, istart, n, pad_length)
    length(out) >= n+2*pad_length || error("output is incorrectly sized")
    x = 2*sig[istart]
    for i = 1:pad_length
        out[ostart+i-1] = x - sig[istart+pad_length+1-i]
    end
    copyto!(out, ostart+pad_length, sig, istart, n)
    x = 2*sig[istart+n-1]
    for i = 1:pad_length
        out[ostart+n+pad_length+i-1] = x - sig[istart+n-1-i]
    end
    out
end

# Zero phase digital filtering by processing data in forward and reverse direction
function iir_filtfilt(b::AbstractVector, a::AbstractVector, x::AbstractArray)
    zi = filt_stepstate(b, a)
    pad_length = 3 * (max(length(a), length(b)) - 1)
    t = Base.promote_eltype(b, a, x)
    zitmp = similar(zi, t)
    extrapolated = Vector{t}(undef, size(x, 1)+pad_length*2)
    out = similar(x, t)

    istart = 1
    for i = 1:Base.trailingsize(x, 2)
        extrapolate_signal!(extrapolated, 1, x, istart, size(x, 1), pad_length)
        reverse!(filt!(extrapolated, b, a, extrapolated, mul!(zitmp, zi, extrapolated[1])))
        filt!(extrapolated, b, a, extrapolated, mul!(zitmp, zi, extrapolated[1]))
        for j = 1:size(x, 1)
            @inbounds out[j, i] = extrapolated[end-pad_length+1-j]
        end
        istart += size(x, 1)
    end

    out
end

"""
    filtfilt(coef, x)

Filter `x` in the forward and reverse directions using filter
coefficients `coef`. The initial state of the filter is computed so
that its response to a step function is steady state. Before
filtering, the data is extrapolated at both ends with an
odd-symmetric extension of length
`3*(max(length(b), length(a))-1)`.

Because `filtfilt` applies the given filter twice, the effective
filter order is twice the order of `coef`. The resulting signal has
zero phase distortion.
"""
function filtfilt(b::AbstractVector, x::AbstractArray)
    nb = length(b)
    # Only need as much padding as the order of the filter
    t = Base.promote_eltype(b, x)
    extrapolated = similar(x, t, size(x, 1)+2nb-2, Base.trailingsize(x, 2))

    # Convolve b with its reverse
    newb = filt(b, reverse(b))
    resize!(newb, 2nb-1)
    for i = 1:nb-1
        newb[nb+i] = newb[nb-i]
    end

    # Extrapolate each column
    for i = 1:size(extrapolated, 2)
        extrapolate_signal!(extrapolated, (i-1)*size(extrapolated, 1)+1, x, (i-1)*size(x, 1)+1, size(x, 1), nb-1)
    end

    # Filter
    out = filt(newb, extrapolated)

    # Drop garbage at start
    reshape(out[2nb-1:end, :], size(x))
end

# Choose whether to use FIR or iir_filtfilt depending on
# length of a
function filtfilt(b::AbstractVector, a::AbstractVector, x::AbstractArray)
    if length(a) == 1
        if a[1] != 1
            b /= a[1]
        end
        filtfilt(b, x)
    else
        iir_filtfilt(b, a, x)
    end
end

# Zero phase digital filtering for second order sections
function filtfilt(f::SecondOrderSections{:z,T,G}, x::AbstractArray{S}) where {T,G,S}
    zi = filt_stepstate(f)
    pad_length = 6 * length(f.biquads)
    t = Base.promote_type(T, G, S)
    zitmp = similar(zi, t)
    extrapolated = Vector{t}(undef, size(x, 1)+pad_length*2)
    out = similar(x, t)

    istart = 1
    for i = 1:Base.trailingsize(x, 2)
        extrapolate_signal!(extrapolated, 1, x, istart, size(x, 1), pad_length)
        reverse!(filt!(extrapolated, f, extrapolated, mul!(zitmp, zi, extrapolated[1])))
        filt!(extrapolated, f, extrapolated, mul!(zitmp, zi, extrapolated[1]))
        for j = 1:size(x, 1)
            @inbounds out[j, i] = extrapolated[end-pad_length+1-j]
        end
        istart += size(x, 1)
    end

    out
end

# Support for other filter types
filtfilt(f::FilterCoefficients{:z}, x) = filtfilt(convert(SecondOrderSections, f), x)
filtfilt(f::PolynomialRatio{:z}, x) = filtfilt(coefb(f), coefa(f), x)

## Initial filter state

# Compute an initial state for filt with coefficients (b,a) such that its
# response to a step function is steady state.
function filt_stepstate(b::Union{AbstractVector{T}, T}, a::Union{AbstractVector{T}, T}) where T<:Number
    scale_factor = a[1]
    if !isone(scale_factor)
        a = a ./ scale_factor
        b = b ./ scale_factor
    end

    bs = length(b)
    as = length(a)
    sz = max(bs, as)
    sz > 0 || error("a and b must have at least one element each")
    sz == 1 && return T[]

    # Pad the coefficients with zeros if needed
    bs<sz && (b = copyto!(zeros(eltype(b), sz), b))
    as<sz && (a = copyto!(zeros(eltype(a), sz), a))

    # construct the companion matrix A and vector B:
    A = [-a[2:end] Matrix{T}(I, sz-1, sz-2)]
    B = @. b[2:end] - a[2:end] * b[1]
    # Solve si = A*si + B
    # (I - A)*si = B
    scale_factor \ (I - A) \ B
 end

function filt_stepstate(f::SecondOrderSections{:z,T}) where T
    biquads = f.biquads
    si = Matrix{T}(undef, 2, length(biquads))
    y = one(T)
    for i = 1:length(biquads)
        biquad = biquads[i]
        a1, a2, b0, b1, b2 = biquad.a1, biquad.a2, biquad.b0, biquad.b1, biquad.b2

        # At steady state, we have:
        #  y = s1 + b0*x
        # s1 = s2 + b1*x - a1*y
        # s2 = b2*x - a2*y
        # where x is the input and y is the output. Solving these
        # equations yields the following.
        den = (1 + a1 + a2)
        si[1, i] = muladd((a1 + a2), -b0, b1 + b2) / den * y
        si[2, i] = muladd(a1, b2, muladd(-a2, (b0 + b1), b2)) / den * y
        y *= (b0 + b1 + b2) / den
    end
    si
end

"""
    tdfilt(h, x)

Apply filter or filter coefficients `h` along the first dimension
of array `x` using a naïve time-domain algorithm
"""
function tdfilt(h::AbstractVector{H}, x::AbstractArray{T}) where {H,T<:Real}
    filt(h, one(H), x)
end

"""
    tdfilt!(out, h, x)

Like `tdfilt`, but writes the result into array `out`. Output array `out` may
not be an alias of `x`, i.e. filtering may not be done in place.
"""
function tdfilt!(out::AbstractArray, h::AbstractVector{H}, x::AbstractArray) where H
    filt!(out, h, one(H), x)
end

filt(h::AbstractVector{H}, x::AbstractArray{T}) where {H,T} =
    filt!(Array{promote_type(H, T)}(undef, size(x)), h, x)

#
# fftfilt and filt
#

"""
    fftfilt(h, x)

Apply FIR filter taps `h` along the first dimension of array `x`
using an FFT-based overlap-save algorithm.
"""
function fftfilt(b::AbstractVector{T}, x::AbstractArray{T},
                 nfft::Integer=optimalfftfiltlength(length(b), length(x))) where T<:Real
    _fftfilt!(Array{T}(undef, size(x)), b, x, nfft)
end

"""
    fftfilt!(out, h, x)

Like `fftfilt` but writes result into out array.
"""
function fftfilt!(
    out::AbstractArray{<:Real},
    b::AbstractVector{<:Real},
    x::AbstractArray{<:Real},
    nfft::Integer=optimalfftfiltlength(length(b), length(x))
)
    size(out) == size(x) || throw(ArgumentError("out and x must be the same size"))
    _fftfilt!(out, b, x, nfft)
end

# Like fftfilt! but does not check if out and x are the same size
function _fftfilt!(
    out::AbstractArray{<:Real},
    b::AbstractVector{H},
    x::AbstractArray{T},
    nfft::Integer
) where {T<:Real,H<:Real}
    nb = length(b)
    nx = size(x, 1)
    normfactor = nfft
    W = promote_type(H, T)

    L = min(nx, nfft - (nb - 1))
    tmp1 = Vector{W}(undef, nfft)
    tmp2 = Vector{Complex{W}}(undef, nfft >> 1 + 1)

    p1 = plan_rfft(tmp1)
    p2 = plan_brfft(tmp2, nfft)

    # FFT of filter
    filterft = similar(tmp2)
    tmp1[1:nb] .= b ./ normfactor
    tmp1[nb+1:end] .= zero(W)
    mul!(filterft, p1, tmp1)

    # FFT of chunks
    for colstart = 0:nx:length(x)-1, off = 1:L:nx
        npadbefore = max(0, nb - off)
        xstart = off - nb + npadbefore + 1
        n = min(nfft - npadbefore, nx - xstart + 1)

        tmp1[1:npadbefore] .= zero(W)
        tmp1[npadbefore+n+1:end] .= zero(W)

        copyto!(tmp1, npadbefore+1, x, colstart+xstart, n)
        mul!(tmp2, p1, tmp1)
        broadcast!(*, tmp2, tmp2, filterft)
        mul!(tmp1, p2, tmp2)

        # Copy to output
        copyto!(out, colstart+off, tmp1, nb, min(L, nx - off + 1))
    end

    out
end

# Filter x using FIR filter b, heuristically choosing to perform convolution in
# the time domain using tdfilt or in the frequency domain using fftfilt
function filt(b::AbstractVector{T}, x::AbstractArray{T}) where T<:Number
    filt_choose_alg!(Array{T}(undef, size(x)), b, x)
end

# Like filt but mutates output array
function filt!(out::AbstractArray, b::AbstractVector, x::AbstractArray)
    size(out) == size(x) || throw(ArgumentError("out must be the same size as x"))
    filt_choose_alg!(out, b, x)
end

# Perform FIR filtering with either time domain or fft based algorithms, writing
# output to `out` argument. Does not check that x and out are the same size.
function filt_choose_alg!(
    out::AbstractArray{<:Real},
    b::AbstractVector{<:Real},
    x::AbstractArray{<:Real}
)
    nb = length(b)

    if nb > SMALL_FILT_CUTOFF
        nx = size(x, 1)
        nfft = optimalfftfiltlength(nb, nx)
        _fftfilt!(out, b, x, nfft)
    else
        tdfilt!(out, b, x)
    end
end

function filt_choose_alg!(out::AbstractArray, b::AbstractVector, x::AbstractArray)
    tdfilt!(out, b, x)
end
