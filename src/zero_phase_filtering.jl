# Zero phase digital filtering for Julia
#
# Contributors:
#     Robert Luke (Robert.Luke@med.kuleuven.be)
#     Matt Bauman (mbauman@gmail.com)
#     Simon Kornblith (simon@simonster.com)

module ZeroPhaseFiltering
using ..FilterDesign
export filtfilt

## filtfilt

# Extrapolate the beginning of a signal for use by filtfilt. This
# computes:
#
# [(2 * x[1]) .- x[pad_length+1:-1:2],
#  x,
#  (2 * x[end]) .- x[end-1:-1:end-pad_length]]
#
# in place in output. The istart and n parameters determine the portion
# of the input signal x to extrapolate.
function extrapolate_signal!(out, sig, istart, n, pad_length)
    length(out) == n+2*pad_length || error("output is incorrectly sized")
    x = 2*sig[istart]
    for i = 1:pad_length
        out[i] = x - sig[istart+pad_length+1-i]
    end
    copy!(out, pad_length+1, sig, istart, n)
    x = 2*sig[istart+n-1]
    for i = 1:pad_length
        out[n+pad_length+i] = x - sig[istart+n-1-i]
    end
    out
end

# Zero phase digital filtering by processing data in forward and reverse direction
function filtfilt(b::AbstractVector, a::AbstractVector, x::AbstractArray)
    zi = filt_stepstate(b, a)
    zitmp = copy(zi)
    pad_length = 3 * (max(length(a), length(b)) - 1)
    t = Base.promote_eltype(b, a, x)
    extrapolated = Array(t, size(x, 1)+pad_length*2)
    out = similar(x, t)

    istart = 1
    for i = 1:size(x, 2)
        extrapolate_signal!(extrapolated, x, istart, size(x, 1), pad_length)
        reverse!(filt!(extrapolated, b, a, extrapolated, scale!(zitmp, zi, extrapolated[1])))
        filt!(extrapolated, b, a, extrapolated, scale!(zitmp, zi, extrapolated[1]))
        for j = 1:size(x, 1)
            @inbounds out[j, i] = extrapolated[end-pad_length+1-j]
        end
        istart += size(x, 1)
    end

    out
end

# Extract si for a biquad, multiplied by a scaling factor
function biquad_si!(zitmp, zi, i, scal)
    zitmp[1] = zi[1, i]*scal
    zitmp[2] = zi[2, i]*scal
    zitmp
end

# Zero phase digital filtering for second order sections
function filtfilt{T,G,S}(f::SOSFilter{T,G}, x::AbstractArray{S})
    zi = filt_stepstate(f)
    zi2 = zeros(2)
    zitmp = zeros(2)
    pad_length = 3*size(zi, 1)
    t = Base.promote_type(T, G, S)
    extrapolated = Array(t, size(x, 1)+pad_length*2)
    out = similar(x, t)

    istart = 1
    for i = 1:size(x, 2)
        # First biquad
        extrapolate_signal!(extrapolated, x, istart, size(x, 1), pad_length)
        f2 = f.biquads[1]*f.g
        reverse!(filt!(extrapolated, f2, extrapolated, biquad_si!(zitmp, zi, 1, extrapolated[1])))
        reverse!(filt!(extrapolated, f2, extrapolated, biquad_si!(zitmp, zi, 1, extrapolated[1])))

        # Subsequent biquads
        for j = 2:length(f.biquads)
            f2 = f.biquads[j]
            extrapolate_signal!(extrapolated, extrapolated, pad_length+1, size(x, 1), pad_length)
            reverse!(filt!(extrapolated, f2, extrapolated, biquad_si!(zitmp, zi, j, extrapolated[1])))
            reverse!(filt!(extrapolated, f2, extrapolated, biquad_si!(zitmp, zi, j, extrapolated[1])))
        end

        # Copy to output
        copy!(out, istart, extrapolated, pad_length+1, size(x, 1))
        istart += size(x, 1)
    end

    out
end

# Support for other filter types
filtfilt(f::Filter, x) = filtfilt(convert(TFFilter, f), x)
filtfilt(f::TFFilter, x) = filtfilt(coefb(f), coefa(f), x)

## Initial filter state

# Compute an initial state for filt with coefficients (b,a) such that its
# response to a step function is steady state.
function filt_stepstate{T<:Number}(b::Union(AbstractVector{T}, T), a::Union(AbstractVector{T}, T))
    scale_factor = a[1]
    if scale_factor != 1.0
        a = a ./ scale_factor
        b = b ./ scale_factor
    end

    bs = length(b)
    as = length(a)
    sz = max(bs, as)
    sz > 0 || error("a and b must have at least one element each")
    sz == 1 && return T[]

    # Pad the coefficients with zeros if needed
    bs<sz && (b = copy!(zeros(eltype(b), sz), b))
    as<sz && (a = copy!(zeros(eltype(a), sz), a))

    # construct the companion matrix A and vector B:
    A = [-a[2:end] [eye(T, sz-2); zeros(T, 1, sz-2)]]
    B = b[2:end] - a[2:end] * b[1]
    # Solve si = A*si + B
    # (I - A)*si = B
    scale_factor \ (I - A) \ B
 end

function filt_stepstate{T}(f::SOSFilter{T})
    biquads = f.biquads
    si = Array(T, 2, length(biquads))
    for i = 1:length(biquads)
        biquad = biquads[i]
        A = [one(T)+biquad.a1 -one(T)
                   +biquad.a2  one(T)]
        B = [biquad.b1 - biquad.a1*biquad.b0,
             biquad.b2 - biquad.a2*biquad.b0]
        si[:, i] = A \ B
    end
    si[1, 1] *= f.g
    si[2, 1] *= f.g
    si
end

end
