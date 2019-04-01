# This file was formerly a part of Julia. License is MIT: https://julialang.org/license

import Base.trailingsize
import Compat.LinearAlgebra.BLAS

_zerosi(b,a,T) = zeros(promote_type(eltype(b), eltype(a), T), max(length(a), length(b))-1)

"""
    filt(b, a, x, [si])

Apply filter described by vectors `a` and `b` to vector `x`, with an optional initial filter
state vector `si` (defaults to zeros).
"""
function filt(b::Union{AbstractVector, Number}, a::Union{AbstractVector, Number},
              x::AbstractArray{T}, si::AbstractArray{S} = _zerosi(b,a,T)) where {T,S}
    filt!(Array{promote_type(eltype(b), eltype(a), T, S)}(undef, size(x)), b, a, x, si)
end

# in-place filtering: returns results in the out argument, which may shadow x
# (and does so by default)

"""
    filt!(out, b, a, x, [si])

Same as [`filt`](@ref) but writes the result into the `out` argument, which may
alias the input `x` to modify it in-place.
"""
function filt!(out::AbstractArray, b::Union{AbstractVector, Number}, a::Union{AbstractVector, Number},
               x::AbstractArray{T}, si::AbstractArray{S,N} = _zerosi(b,a,T)) where {T,S,N}
    isempty(b) && throw(ArgumentError("filter vector b must be non-empty"))
    isempty(a) && throw(ArgumentError("filter vector a must be non-empty"))
    a[1] == 0  && throw(ArgumentError("filter vector a[1] must be nonzero"))
    if size(x) != size(out)
        throw(ArgumentError("output size $(size(out)) must match input size $(size(x))"))
    end

    as = length(a)
    bs = length(b)
    sz = max(as, bs)
    silen = sz - 1
    ncols = trailingsize(x,2)

    if size(si, 1) != silen
        throw(ArgumentError("initial state vector si must have max(length(a),length(b))-1 rows"))
    end
    if N > 1 && trailingsize(si,2) != ncols
        throw(ArgumentError("initial state vector si must be a vector or have the same number of columns as x"))
    end

    size(x,1) == 0 && return out
    sz == 1 && return mul!(out, x, b[1]/a[1]) # Simple scaling without memory

    # Filter coefficient normalization
    if a[1] != 1
        norml = a[1]
        a = a ./ norml
        b = b ./ norml
    end
    # Pad the coefficients with zeros if needed
    bs<sz   && (b = copyto!(zeros(eltype(b), sz), b))
    1<as<sz && (a = copyto!(zeros(eltype(a), sz), a))

    initial_si = si
    for col = 1:ncols
        # Reset the filter state
        si = initial_si[:, N > 1 ? col : 1]
        if as > 1
            _filt_iir!(out, b, a, x, si, col)
        else
            _filt_fir!(out, b, x, si, col)
        end
    end
    return out
end

function _filt_iir!(out, b, a, x, si, col)
    silen = length(si)
    @inbounds for i=1:size(x, 1)
        xi = x[i,col]
        val = si[1] + b[1]*xi
        for j=1:(silen-1)
            si[j] = si[j+1] + b[j+1]*xi - a[j+1]*val
        end
        si[silen] = b[silen+1]*xi - a[silen+1]*val
        out[i,col] = val
    end
end

function _filt_fir!(out, b, x, si, col)
    silen = length(si)
    @inbounds for i=1:size(x, 1)
        xi = x[i,col]
        val = si[1] + b[1]*xi
        for j=1:(silen-1)
            si[j] = si[j+1] + b[j+1]*xi
        end
        si[silen] = b[silen+1]*xi
        out[i,col] = val
    end
end

"""
    deconv(b,a) -> c

Construct vector `c` such that `b = conv(a,c) + r`.
Equivalent to polynomial division.
"""
function deconv(b::StridedVector{T}, a::StridedVector{T}) where T
    lb = size(b,1)
    la = size(a,1)
    if lb < la
        return [zero(T)]
    end
    lx = lb-la+1
    x = zeros(T, lx)
    x[1] = 1
    filt(b, a, x)
end

function _zeropad!(padded::AbstractVector, u::AbstractVector)
    ulen = length(u)
    padlen = length(padded)
    copyto!(padded, 1, u, 1, ulen)
    padded[ulen + 1:padlen] .= 0
    padded
end
function _zeropad!(padded::AbstractArray{<:Any, N}, u::AbstractArray{<:Any, N}) where N
    datainds = CartesianIndices(u)
    copyto!(padded, datainds, u, datainds)
    pad_ranges = Vector{UnitRange{Int}}(undef, N)
    pad_ranges .= axes(padded)
    for i = N:-1:1
        nu = size(u, i)
        pad_ranges[i] = nu + 1:size(padded, i)
        padinds = CartesianIndices(NTuple{N, UnitRange{Int}}(pad_ranges))
        padded[padinds] .= 0
        pad_ranges[i] = 1:nu
    end
    padded
end
function _zeropad(u, padded_size)
    _zeropad!(similar(u, padded_size), u)
end

function _conv(
    u::StridedArray{T, N}, v::StridedArray{T, N}, paddims
) where {T<:Real, N}
    padded = _zeropad(u, paddims)
    p = plan_rfft(padded)
    uf = p * padded
    _zeropad!(padded, v)
    vf = p * padded
    uf .*= vf
    irfft(uf, paddims[1])
end
function _conv(u, v, paddims)
    upad = _zeropad(u, paddims)
    vpad = _zeropad(v, paddims)
    p! = plan_fft!(upad)
    p! * upad # Operates in place on upad
    p! * vpad
    upad .*= vpad
    ifft!(upad)
end

_conv_clip!(y::AbstractVector, minpad) = resize!(y, minpad[1])
_conv_clip!(y::AbstractArray, minpad) = y[CartesianIndices(minpad)]

"""
    conv(u,v)

Convolution of two arrays. Uses FFT algorithm.
"""
function conv(u::StridedArray{T, N}, v::StridedArray{T, N}) where {T<:BLAS.BlasFloat, N}
    su = size(u)
    sv = size(v)
    minpad = su .+ sv .- 1
    nfft = map(n -> n > 1024 ? nextprod([2,3,5], n) : nextpow(2, n), minpad)
    padsize = NTuple{N, Int}(nfft)
    y = _conv(u, v, padsize)
    _conv_clip!(y, minpad)
end

conv(u::StridedArray{T, N}, v::StridedArray{T, N}) where {T<:Integer, N} =
    round.(Int, conv(float(u), float(v)))
function conv(
    u::StridedArray{<:Integer, N}, v::StridedArray{<:BLAS.BlasFloat, N}
) where N
    conv(float(u), v)
end
function conv(
    u::StridedArray{<:BLAS.BlasFloat, N}, v::StridedArray{<:Integer, N}
) where N
    conv(u, float(v))
end

function conv(A::StridedArray{T}, B::StridedArray{T}) where T
    maxnd = max(ndims(A), ndims(B))
    return conv(cat(A, dims=maxnd), cat(B, dims=maxnd))
end

"""
    conv(u,v,A)

2-D convolution of the matrix `A` with the 2-D separable kernel generated by
the vectors `u` and `v`.
Uses 2-D FFT algorithm.
"""
function conv(u::StridedVector{T}, v::StridedVector{T}, A::StridedMatrix{T}) where T
    m = length(u)+size(A,1)-1
    n = length(v)+size(A,2)-1
    B = zeros(T, m, n)
    B[1:size(A,1),1:size(A,2)] = A
    u = fft([u;zeros(T,m-length(u))])
    v = fft([v;zeros(T,n-length(v))])
    C = ifft(fft(B) .* (u * transpose(v)))
    if T <: Real
        return real(C)
    end
    return C
end


function check_padmode_kwarg(padmode::Symbol, su::Integer, sv::Integer)
    if padmode == :default_longest
        if su != sv
            Base.depwarn(
            """
            The default value of `padmode` will be changing from `:longest` to
            `:none` in a future release of DSP. In preparation for this change,
            leaving `padmode` unspecified is currently deprecated. To keep
            current behavior specify `padmode=:longest`. To avoid this warning,
            specify padmode = :none or padmode = :longest where appropriate.
            """
                ,
                :xcorr
            )
        end
        :longest
    else
        padmode
    end
end

"""
    xcorr(u,v; padmode = :longest)

Compute the cross-correlation of two vectors. The size of the output depends on
the padmode keyword argument: with padmode = :none the length of the
result will be length(u) + length(v) - 1, as with conv. With
padmode = :longest the shorter of the arguments will be padded so they
are equal length. This gives a result with length 2*max(length(u), length(v))-1,
with the zero-lag condition at the center.

!!! warning
    The default value of `padmode` will be changing from `:longest` to `:none`
    in a future release of DSP. In preparation for this change, leaving
    `padmode` unspecified is currently deprecated.
"""
function xcorr(u, v; padmode::Symbol = :default_longest)
    su = size(u,1); sv = size(v,1)
    padmode = check_padmode_kwarg(padmode, su, sv)
    if padmode == :longest
        if su < sv
            u = _zeropad(u, sv)
        elseif sv < su
            v = _zeropad(v, su)
        end
        conv(u, Compat.reverse(conj(v), dims=1))
    elseif padmode == :none
        conv(u, Compat.reverse(conj(v), dims=1))
    else
        throw(ArgumentError("padmode keyword argument must be either :none or :longest"))
    end
end
