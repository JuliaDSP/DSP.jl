# This file was formerly a part of Julia. License is MIT: https://julialang.org/license

import Base.trailingsize
import LinearAlgebra.BLAS

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

# padded must start at index 1, but u can have arbitrary offset
function _zeropad!(padded::AbstractVector, u::AbstractVector)
    datasize = length(u)
    # Use axes to accommodate arrays that do not start at index 1
    padsize = length(padded)
    data_first_i = first(axes(u, 1))
    copyto!(padded, 1, u, data_first_i, datasize)
    padded[1 + datasize:padsize] .= 0
    padded
end

# padded must start at index 1, but u can have arbitrary offset
function _zeropad!(padded::AbstractArray{<:Any, N},
                   u::AbstractArray{<:Any, N}) where N
    # Copy the data to the beginning of the padded array
    padsize = size(padded)
    pad_ax = axes(padded)
    datasize = size(u)
    pad_data_ranges = range.(1, datasize)
    copyto!(padded, CartesianIndices(pad_data_ranges), u, CartesianIndices(u))

    # Step through each dimension, and fill the trailing indices in that
    # dimension with zeros.
    #
    # This is accomplished by going from the last to the first dimension, and
    # for each dimension, finding the biggest rectangular region that needs to
    # be zero padded. This region corresponds to the  trailing indices after the
    # data in the selected dimension, all the indices in lower dimensions, and
    # the same indices as the data for higher dimensions, corresponding to
    # regions that were missed in that dimension. This should allow column-major
    # friendly access to the array.
    #
    # Geometric intuition for two dimensions:
    #   1 1 1 0 0     1 1 1    0 0
    #   1 1 1 0 0     1 1 1    0 0
    #   1 1 1 0 0  =  1 1 1    0 0
    #   0 0 0 0 0              0 0
    #   0 0 0 0 0     0 0 0    0 0
    #                 0 0 0

    pad_ranges = Vector{UnitRange{Int}}(undef, N)
    pad_ranges .= pad_ax # Initially select the entire dimension
    for i = N:-1:1
        # Trailing indices for this dimension
        pad_ranges[i] = datasize[i] + 1 : padsize[i]

        # Make the rectangular region and set it to zero
        pad_region = CartesianIndices(NTuple{N, UnitRange{Int}}(pad_ranges))
        padded[pad_region] .= 0

        # Set this dimension to be just the portion that was missed (same as data)
        pad_ranges[i] = pad_data_ranges[i]
    end
    padded
end
function _zeropad(u, padded_size)
    _zeropad!(similar(u, padded_size), u)
end

function _conv(
    u::AbstractArray{T, N}, v::AbstractArray{T, N}, paddims
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

_conv_clip!(y::AbstractVector, minpad) = sizehint!(resize!(y, minpad[1]), minpad[1])
_conv_clip!(y::AbstractArray, minpad) = y[CartesianIndices(minpad)]

"""
    conv(u,v)

Convolution of two arrays. Uses FFT algorithm.
"""
function conv(u::AbstractArray{T, N},
              v::AbstractArray{T, N}) where {T<:BLAS.BlasFloat, N}
    su = size(u)
    sv = size(v)
    minpad = su .+ sv .- 1
    padsize = map(n -> n > 1024 ? nextprod([2,3,5], n) : nextpow(2, n), minpad)
    y = _conv(u, v, padsize)
    _conv_clip!(y, minpad)
end
function conv(u::AbstractArray{<:BLAS.BlasFloat, N},
              v::AbstractArray{<:BLAS.BlasFloat, N}) where N
    fu, fv = promote(u, v)
    conv(fu, fv)
end
conv(u::AbstractArray{T, N}, v::AbstractArray{T, N}) where {T<:Number, N} =
    conv(float(u), float(v))
conv(u::AbstractArray{T, N}, v::AbstractArray{T, N}) where {T<:Integer, N} =
    round.(Int, conv(float(u), float(v)))
function conv(u::AbstractArray{<:Number, N},
              v::AbstractArray{<:BLAS.BlasFloat, N}) where N
    conv(float(u), v)
end
function conv(u::AbstractArray{<:BLAS.BlasFloat, N},
              v::AbstractArray{<:Number, N}) where N
    conv(u, float(v))
end

function conv(A::AbstractArray{T}, B::AbstractArray{T}) where T
    maxnd = max(ndims(A), ndims(B))
    return conv(cat(A, dims=maxnd), cat(B, dims=maxnd))
end

"""
    conv(u,v,A)

2-D convolution of the matrix `A` with the 2-D separable kernel generated by
the vectors `u` and `v`.
Uses 2-D FFT algorithm.
"""
function conv(u::AbstractVector{T}, v::AbstractVector{T}, A::AbstractMatrix{T}) where T
    # Arbitrary indexing offsets not implemented
    @assert !Base.has_offset_axes(u, v, A)
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
function xcorr(
    u::AbstractVector, v::AbstractVector; padmode::Symbol = :default_longest
)
    su = size(u,1); sv = size(v,1)
    padmode = check_padmode_kwarg(padmode, su, sv)
    if padmode == :longest
        if su < sv
            u = _zeropad(u, sv)
        elseif sv < su
            v = _zeropad(v, su)
        end
        conv(u, reverse(conj(v), dims=1))
    elseif padmode == :none
        conv(u, reverse(conj(v), dims=1))
    else
        throw(ArgumentError("padmode keyword argument must be either :none or :longest"))
    end
end
