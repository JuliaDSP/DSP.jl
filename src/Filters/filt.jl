# Implementations of filt, filtfilt, and fftfilt

#
# filt and filt!
#
using ..DSP: _filt_fir!, _filt_iir!

## PolynomialRatio

"""
    filt!(out, f, x)

Same as [`filt()`](@ref) but writes the result into the `out`
argument. Output array `out` may not be an alias of `x`, i.e. filtering may
not be done in place.
"""
filt!(out, f::PolynomialRatio{:z}, x::AbstractArray) = filt!(out, coefb(f), coefa(f), x)

"""
    filt(f::FilterCoefficients{:z}, x::AbstractArray)

Apply filter or filter coefficients `f` along the first dimension
of array `x`. If `f` is a `PolynomialRatio`, `Biquad`, or
`SecondOrderSections`, filtering is implemented directly. If
`f` is a `ZeroPoleGain` object, it is first converted to a
`SecondOrderSections` object.  If `f` is a Vector, it is
interpreted as an FIR filter, and a naïve or FFT-based algorithm is
selected based on the data and filter length.
"""
filt(f::PolynomialRatio{:z}, x) = filt(coefb(f), coefa(f), x)

## SecondOrderSections

# filt! algorithm (no checking, returns si)
function _filt!(out::AbstractArray, si::AbstractArray{S,N}, f::SecondOrderSections{:z},
                x::AbstractArray, col::Union{Int,CartesianIndex}) where {S,N}
    g = f.g
    biquads = f.biquads
    @inbounds for i in axes(x, 1)
        yi = x[i, col]
        for fi in eachindex(biquads)
            biquad = biquads[fi]
            xi = yi
            yi = muladd(biquad.b0, xi, si[1, fi])
            si[1, fi] = muladd(biquad.a1, -yi, muladd(biquad.b1, xi, si[2, fi]))
            si[2, fi] = muladd(biquad.b2, xi, -biquad.a2 * yi)
        end
        out[i, col] = g * yi
    end
    si
end

function filt!(out::AbstractArray, f::SecondOrderSections{:z,T,G}, x::AbstractArray{S}) where {T,G,S}
    size(x) != size(out) && throw(DimensionMismatch("out size must match x"))

    si = Matrix{promote_type(T, G, S)}(undef, 2, length(f.biquads))
    for col in CartesianIndices(axes(x)[2:end])
        fill!(si, zero(eltype(si)))
        _filt!(out, si, f, x, col)
    end
    out
end

filt(f::SecondOrderSections{:z,T,G}, x::AbstractArray{S}) where {T,G,S<:Number} =
    filt!(similar(x, promote_type(T, G, S)), f, x)

## Biquad

# filt! algorithm (no checking, returns si)
function _filt!(out::AbstractArray, si1::Number, si2::Number, f::Biquad{:z},
                x::AbstractArray, col::Union{Int,CartesianIndex})
    @inbounds for i in axes(x, 1)
        xi = x[i, col]
        yi = muladd(f.b0, xi, si1)
        si1 = muladd(f.a1, -yi, muladd(f.b1, xi, si2))
        si2 = muladd(f.b2, xi, -f.a2 * yi)
        out[i, col] = yi
    end
    (si1, si2)
end

function filt!(out::AbstractArray, f::Biquad{:z,T}, x::AbstractArray{S}) where {T,S}
    size(x) != size(out) && throw(DimensionMismatch("out size must match x"))

    for col in CartesianIndices(axes(x)[2:end])
        _filt!(out, zero(promote_type(T, S)), zero(promote_type(T, S)), f, x, col)
    end
    out
end

filt(f::Biquad{:z,T}, x::AbstractArray{S}) where {T,S<:Number} =
    filt!(similar(x, promote_type(T, S)), f, x)

## For arbitrary filters, convert to SecondOrderSections
filt(f::FilterCoefficients{:z}, x) = filt(convert(SecondOrderSections, f), x)
filt!(out, f::FilterCoefficients{:z}, x) = filt!(out, convert(SecondOrderSections, f), x)

"""
    DF2TFilter(coef::FilterCoefficients{:z})
    DF2TFilter(coef::FilterCoefficients{:z}, coldims::Tuple)
    DF2TFilter(coef::FilterCoefficients{:z}, sitype::Type, coldims::Tuple = ())
    DF2TFilter(coef::FilterCoefficients{:z}, si)

Construct a stateful direct form II transposed filter with
coefficients `coef`.

The initial filter state defaults to zeros (of a type derived from `coef`)
suitable for vector input. Another element type of the state can be specified
with `sitype`.

To allow column-wise filtering of higher-dimensional input, the size of the
extra dimensions have to be given in `coldims`. To e.g. column-wise filter an
input with size `(L, N1, N2)`, set `coldims` to `(N1, N2)`.

Alternatively, an array representing the initial filter state can be passed
as `si`.

If `coef` is a `PolynomialRatio`, `Biquad`, or `SecondOrderSections`,
filtering is implemented directly. If `coef` is a `ZeroPoleGain`
object, it is first converted to a `SecondOrderSections` object.
"""
struct DF2TFilter{T<:FilterCoefficients{:z},S<:Array}
    coef::T
    state::S

    function DF2TFilter{Ti,Si}(coef::PolynomialRatio{:z}, state::Array) where {Ti,Si}
        size(state, 1) == max(length(coefa(coef)), length(coefb(coef)))-1 ||
            throw(ArgumentError("length of state vector must match filter order"))
        new{Ti,Si}(coef, state)
    end
    function DF2TFilter{Ti,Si}(coef::SecondOrderSections{:z}, state::Array) where {Ti,Si}
        size(state, 1) == 2 && size(state, 2) == length(coef.biquads) ||
            throw(ArgumentError("state must be 2 x nbiquads"))
        new{Ti,Si}(coef, state)
    end
    function DF2TFilter{Ti,Si}(coef::Biquad{:z}, state::Array) where {Ti,Si}
        size(state, 1) == 2 || throw(ArgumentError("length of state must be 2"))
        new{Ti,Si}(coef, state)
    end
end

DF2TFilter(coef::Union{PolynomialRatio{:z,T},Biquad{:z,T}}, state::S) where {T,S<:Array} =
    DF2TFilter{typeof(coef),S}(coef, state)

DF2TFilter(coef::SecondOrderSections{:z,T,G}, state::S) where {T,G,S<:Array} =
    DF2TFilter{SecondOrderSections{:z,T,G},S}(coef, state)

## PolynomialRatio
DF2TFilter(coef::PolynomialRatio{:z,T}, coldims::Tuple{Vararg{Integer}}) where {T} = DF2TFilter(coef, T, coldims)
DF2TFilter(coef::PolynomialRatio{:z,T}, ::Type{V}=T, coldims::Tuple{Vararg{Integer}}=()) where {T,V} =
    DF2TFilter(coef, zeros(promote_type(T, V), max(length(coefa(coef)), length(coefb(coef))) - 1, coldims...))

function filt!(out::AbstractArray{<:Any,N}, f::DF2TFilter{<:PolynomialRatio,Array{T,N}} where T, x::AbstractArray{<:Any,N}) where N
    size(x) != size(out) && throw(ArgumentError("out size must match x"))

    si = f.state
    size(x)[2:end] != size(si)[2:end] && throw(ArgumentError("state size must match x"))

    n = size(si, 1) + 1
    b = coefb(f.coef)
    if n == 1
        mul!(out, x, b[1])
    else
        a = coefa(f.coef)
        if length(a) != 1
            for col in CartesianIndices(axes(x)[2:end])
                _filt_iir!(out, b, a, x, view(si, :, col), col)
            end
        elseif n <= SMALL_FILT_CUTOFF
            vtup = ntuple(j -> b[j], Val(length(b)))
            for col in CartesianIndices(axes(x)[2:end])
                _filt_fir!(out, vtup, x, view(si, :, col), col, Val(true))
            end
        else
            for col in CartesianIndices(axes(x)[2:end])
                _filt_fir!(out, b, x, view(si, :, col), col)
            end
        end
    end
    return out
end

## SecondOrderSections
DF2TFilter(coef::SecondOrderSections{:z,T}, coldims::Tuple{Vararg{Integer}}) where {T} = DF2TFilter(coef, T, coldims)
DF2TFilter(coef::SecondOrderSections{:z,T,G}, ::Type{V}=T, coldims::Tuple{Vararg{Integer}}=()) where {T,G,V} =
    DF2TFilter(coef, zeros(promote_type(T, G, V), 2, length(coef.biquads), coldims...))

function filt!(out::AbstractArray{<:Any,N}, f::DF2TFilter{<:SecondOrderSections,<:Array}, x::AbstractArray{<:Any,N}) where N
    size(x) != size(out) && throw(ArgumentError("out size must match x"))
    size(x)[2:end] != size(f.state)[3:end] && throw(ArgumentError("state size must match x"))
    for col in CartesianIndices(axes(x)[2:end])
        _filt!(out, view(f.state, :, :, col), f.coef, x, col)
    end
    return out
end

## Biquad
DF2TFilter(coef::Biquad{:z,T}, coldims::Tuple{Vararg{Integer}}) where {T} = DF2TFilter(coef, T, coldims)
DF2TFilter(coef::Biquad{:z,T}, ::Type{V}=T, coldims::Tuple{Vararg{Integer}}=()) where {T,V} =
    DF2TFilter(coef, zeros(promote_type(T, V), 2, coldims...))

function filt!(out::AbstractArray{<:Any,N}, f::DF2TFilter{<:Biquad,Array{T,N}} where T, x::AbstractArray{<:Any,N}) where N
    size(x) != size(out) && throw(ArgumentError("out size must match x"))
    si = f.state
    size(x)[2:end] != size(si)[2:end] && throw(ArgumentError("state size must match x"))
    for col in CartesianIndices(axes(x)[2:end])
        (si[1, col], si[2, col]) = _filt!(out, si[1, col], si[2, col], f.coef, x, col)
    end
    return out
end

"""
    filt(f::DF2TFilter{<:FilterCoefficients{:z},<:Array{T}}, x::AbstractArray{V}) where {T,V}

Apply the [stateful filter](@ref stateful-filter-objects) `f` on `x`.

!!! warning
    The output array has eltype `promote_type(T, V)`, where
    `T` is the eltype of the filter state.\n
    For more control over the output type, provide a preallocated
    output array `out` to `filt!(out, f, x)`.
"""
filt(f::DF2TFilter{<:FilterCoefficients{:z},<:Array{T}}, x::AbstractArray{V}) where {T,V} =
    filt!(similar(x, promote_type(T, V)), f, x)

# Fall back to SecondOrderSections
DF2TFilter(coef::FilterCoefficients{:z}, coldims::Tuple{Vararg{Integer}}=()) =
    DF2TFilter(convert(SecondOrderSections, coef), coldims)
DF2TFilter(coef::FilterCoefficients{:z}, ::Type{V}, coldims::Tuple{Vararg{Integer}}=()) where {V} =
    DF2TFilter(convert(SecondOrderSections, coef), V, coldims)

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
function extrapolate_signal!(out, sig, pad_length)
    n = length(sig)
    length(out) >= n + 2 * pad_length || throw(ArgumentError("output is incorrectly sized"))
    x = 2 * sig[1]
    for i = 1:pad_length
        out[i] = x - sig[2+pad_length-i]
    end
    copyto!(out, 1+pad_length, sig, 1, n)
    x = 2 * sig[n]
    for i = 1:pad_length
        out[n+pad_length+i] = x - sig[n-i]
    end
    out
end

# Zero phase digital filtering by processing data in forward and reverse direction
function iir_filtfilt(b::AbstractVector, a::AbstractVector, x::AbstractArray)
    pad_length = min(3 * (max(length(a), length(b)) - 1), size(x, 1) - 1)
    zi, bn, an = filt_stepstate(b, a)
    t = Base.promote_eltype(bn, an, x)
    zitmp = similar(zi, t)
    extrapolated = Vector{t}(undef, size(x, 1) + 2 * pad_length)
    out = similar(x, t)

    for col in CartesianIndices(axes(x)[2:end])
        extrapolate_signal!(extrapolated, view(x, :, col), pad_length)
        _filt_iir!(extrapolated, bn, an, extrapolated, mul!(zitmp, zi, extrapolated[1]), CartesianIndex())
        reverse!(extrapolated)
        _filt_iir!(extrapolated, bn, an, extrapolated, mul!(zitmp, zi, extrapolated[1]), CartesianIndex())
        for j in axes(x, 1)
            out[j, col] = extrapolated[end-pad_length+1-j]
        end
    end

    out
end

"""
    filtfilt(coef::FilterCoefficients, x::AbstractArray)
    filtfilt(b::AbstractVector, x::AbstractArray)
    filtfilt(b::AbstractVector, a::AbstractVector, x::AbstractArray)

Filter `x` in the forward and reverse directions using either a
`FilterCoefficients` object `coef`, or the coefficients `b` and optionally `a`
as in [`filt`](@ref). The initial state of the filter is computed so
that its response to a step function is steady state. Before
filtering, the data is extrapolated at both ends with an
odd-symmetric extension of length
`min(3*(max(length(b), length(a))-1), size(x, 1) - 1)`.

Because `filtfilt` applies the given filter twice, the effective
filter order is twice the order of `coef`. The resulting signal has
zero phase distortion.
"""
function filtfilt end

function filtfilt(b::AbstractVector, x::AbstractArray)
    nb = length(b)
    # Only need as much padding as the order of the filter
    T = Base.promote_eltype(b, x)
    extrapolated = similar(x, T, (size(x, 1)+2*(nb-1), size(x)[2:end]...))

    # Convolve b with its reverse
    newb = reverse(b, 1)
    filt!(newb, b, newb)
    resize!(newb, 2nb-1)
    for i = 1:nb-1
        newb[nb+i] = newb[nb-i]
    end

    # Extrapolate each column
    for col in CartesianIndices(axes(x)[2:end])
        @views extrapolate_signal!(extrapolated[:, col], x[:, col], nb-1)
    end

    # Filter
    filt!(extrapolated, newb, extrapolated)

    # Drop garbage at start
    return extrapolated[2nb-1:end, axes(extrapolated)[2:end]...]
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
    pad_length = min(6 * length(f.biquads), size(x, 1) - 1)
    t = promote_type(T, G, S)
    zitmp = similar(zi, t)
    extrapolated = Vector{t}(undef, size(x, 1) + pad_length * 2)
    out = similar(x, t)

    for col in CartesianIndices(axes(x)[2:end])
        extrapolate_signal!(extrapolated, view(x, :, col), pad_length)
        _filt!(extrapolated, mul!(zitmp, zi, extrapolated[1]), f, extrapolated, CartesianIndex())
        reverse!(extrapolated)
        _filt!(extrapolated, mul!(zitmp, zi, extrapolated[1]), f, extrapolated, CartesianIndex())
        for j in axes(x, 1)
            out[j, col] = extrapolated[end-pad_length+1-j]
        end
    end

    out
end

# Support for other filter types
filtfilt(f::FilterCoefficients{:z}, x) = filtfilt(convert(SecondOrderSections, f), x)
filtfilt(f::PolynomialRatio{:z}, x) = filtfilt(coefb(f), coefa(f), x)

## Initial filter state

# Compute an initial state for filt with coefficients (b,a) such that its
# response to a step function is steady state. Also returns padded (b, a).
function filt_stepstate(b::AbstractVector{V}, a::AbstractVector{V}) where V<:Number
    T = typeof(one(V) / one(V))
    scale_factor = a[1]
    if !isone(scale_factor)
        a = a ./ scale_factor
        b = b ./ scale_factor
    elseif T !== V
        a = convert.(T, a)
        b = convert.(T, b)
    end

    bs = length(b)
    as = length(a)
    sz = max(bs, as)
    sz > 0 || throw(ArgumentError("a and b must have at least one element each"))

    # Pad the coefficients with zeros if needed
    bs < sz && (b = copyto!(zeros(T, sz), b))
    as < sz && (a = copyto!(zeros(T, sz), a))
    sz == 1 && return (T[], b, a)

    # construct the companion matrix A and vector B:
    A = [-a[2:end] Matrix{T}(I, sz-1, sz-2)]
    B = @views @. muladd(a[2:end], -b[1], b[2:end])
    # Solve si = A*si + B
    # (I - A)*si = B
    si = (((I - A) \ B) .*= scale_factor)
    return (si, b, a)
end

filt_stepstate(b::AbstractVector{<:Number}, a::AbstractVector{<:Number}) =
    filt_stepstate(promote(b, a)...)

function filt_stepstate(f::SecondOrderSections{:z,T}) where T
    biquads = f.biquads
    si = Matrix{T}(undef, 2, length(biquads))
    y = one(T)
    for i in eachindex(biquads)
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
    tdfilt(h::AbstractVector, x::AbstractArray)

Apply filter or filter coefficients `h` along the first dimension
of array `x` using a naïve time-domain algorithm
"""
function tdfilt(h::AbstractVector{H}, x::AbstractArray{T}) where {H,T<:Real}
    filt(h, one(H), x)
end

"""
    tdfilt!(out::AbstractArray, h::AbstractVector, x::AbstractArray)

Like `tdfilt`, but writes the result into array `out`. Output array `out` may
not be an alias of `x`, i.e. filtering may not be done in place.
"""
function tdfilt!(out::AbstractArray, h::AbstractVector{H}, x::AbstractArray) where H
    filt!(out, h, one(H), x)
end

filt(h::AbstractVector{H}, x::AbstractArray{T}) where {H,T} =
    filt!(similar(x, promote_type(H, T)), h, x)

#
# fftfilt and filt
#

"""
    fftfilt(h::AbstractVector{<:Real}, x::AbstractArray{<:Real})

Apply FIR filter taps `h` along the first dimension of array `x`
using an FFT-based overlap-save algorithm.
"""
function fftfilt(b::AbstractVector{H}, x::AbstractArray{T},
                 nfft::Integer=optimalfftfiltlength(length(b), length(x))) where {H<:Real,T<:Real}
    _fftfilt!(similar(x, promote_type(H, T)), b, x, nfft)
end

"""
    fftfilt!(out::AbstractArray, h::AbstractVector, x::AbstractArray)

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

        fill!(tmp1, zero(W))
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
    filt_choose_alg!(similar(x, T), b, x)
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
