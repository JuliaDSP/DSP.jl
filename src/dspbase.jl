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

function _zeropad(u, ntot, nu = length(u))
    padded = similar(u, ntot)
    copyto!(padded, 1, u, 1, nu)
    padded[nu + 1:ntot] .= 0
    padded
end

function _circ_conv(
    upad::StridedVector{T}, vpad::StridedVector{T}, np2::Integer
) where T<:Real
    p = plan_rfft(upad)
    irfft((p*upad).*(p*vpad), np2)
end
function _circ_conv(upad, vpad, ::Integer)
    p = plan_fft!(upad)
    ifft!((p*upad).*(p*vpad))
end

"""
    conv(u,v)

Convolution of two vectors. Uses FFT algorithm.
"""
function conv(u::StridedVector{T}, v::StridedVector{T}) where T<:BLAS.BlasFloat
    nu = length(u)
    nv = length(v)
    n = nu + nv - 1
    np2 = n > 1024 ? nextprod([2,3,5], n) : nextpow(2, n)
    upad = _zeropad(u, np2, nu)
    vpad = _zeropad(v, np2, nv)
    y = _circ_conv(upad, vpad, np2)
    return y[1:n]
end
conv(u::StridedVector{T}, v::StridedVector{T}) where {T<:Integer} = round.(Int, conv(float(u), float(v)))
conv(u::StridedVector{<:Integer}, v::StridedVector{<:BLAS.BlasFloat}) = conv(float(u), v)
conv(u::StridedVector{<:BLAS.BlasFloat}, v::StridedVector{<:Integer}) = conv(u, float(v))

"""
    conv2(u,v,A)

2-D convolution of the matrix `A` with the 2-D separable kernel generated by
the vectors `u` and `v`.
Uses 2-D FFT algorithm.
"""
function conv2(u::StridedVector{T}, v::StridedVector{T}, A::StridedMatrix{T}) where T
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

"""
    conv2(B,A)

2-D convolution of the matrix `B` with the matrix `A`. Uses 2-D FFT algorithm.
"""
function conv2(A::StridedMatrix{T}, B::StridedMatrix{T}) where T
    sa, sb = size(A), size(B)
    At = zeros(T, sa[1]+sb[1]-1, sa[2]+sb[2]-1)
    Bt = zeros(T, sa[1]+sb[1]-1, sa[2]+sb[2]-1)
    At[1:sa[1], 1:sa[2]] = A
    Bt[1:sb[1], 1:sb[2]] = B
    p = plan_fft(At)
    C = ifft((p*At).*(p*Bt))
    if T <: Real
        return real(C)
    end
    return C
end
conv2(A::StridedMatrix{T}, B::StridedMatrix{T}) where {T<:Integer} =
    round.(Int, conv2(float(A), float(B)))
conv2(u::StridedVector{T}, v::StridedVector{T}, A::StridedMatrix{T}) where {T<:Integer} =
    round.(Int, conv2(float(u), float(v), float(A)))

"""
    xcorr(u,v)

Compute the cross-correlation of two vectors.
"""
function xcorr(u, v)
    su = size(u,1); sv = size(v,1)
    if su < sv
        u = [u;zeros(eltype(u),sv-su)]
    elseif sv < su
        v = [v;zeros(eltype(v),su-sv)]
    end
    conv(u, Compat.reverse(conj(v), dims=1))
end
