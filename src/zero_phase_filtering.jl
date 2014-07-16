# Zero phase digital filtering for Julia
#
# Contributors:
#     Robert Luke (Robert.Luke@med.kuleuven.be)
#     Matt Bauman (mbauman@gmail.com)
#     Simon Kornblith (simon@simonster.com)

module ZeroPhaseFiltering

using ..FilterDesign

export filtfilt

##############
#
# filtfilt
#
##############

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
function filtfilt{T}(b::AbstractVector, a::AbstractVector, x::AbstractArray{T})
    zi = filt_stepstate(b, a)
    zitmp = copy(zi)
    pad_length = 3 * (max(length(a), length(b)) - 1)
    extrapolated = Array(T, size(x, 1)+pad_length*2)
    out = similar(x)

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

# Support for filter types
filtfilt(f::Filter, x) = filtfilt(convert(TFFilter, f), x)
filtfilt(f::TFFilter, x) = filtfilt(coeffs(f.b), coeffs(f.a), x)


##############
#
# Determine initial state from Matt Bauman
# https://github.com/mbauman
#
##############

# Compute an initial state for filt with coefficients (b,a) such that its
# response to a step function is steady state.
function filt_stepstate{T<:Number}(b::Union(AbstractVector{T}, T), a::Union(AbstractVector{T}, T))

    scale_factor = a[1]
    if scale_factor != 1.0
        a = a ./ scale_factor
        b = b ./ scale_factor
    end

    sz = length(a)
    sz == length(b) || error("a and b must be the same length")
    sz > 0 || error("a and b must have at least one element each")
    a[1] == 1 || error("a and b must be normalized such that a[1] == 1")

    sz == 1 && return T[]

    # construct the companion matrix A and vector B:
    A = [-a[2:end] [eye(T, sz-2); zeros(T, 1, sz-2)]]
    B = b[2:end] - a[2:end] * b[1]
    # Solve si = A*si + B
    # (I - A)*si = B
    scale_factor \ (eye(size(A)[1]) - A) \ B
 end


end
