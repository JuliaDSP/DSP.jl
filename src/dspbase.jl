# This file was formerly a part of Julia. License is MIT: https://julialang.org/license

import Base.trailingsize
import LinearAlgebra.BLAS

const SMALL_FILT_CUTOFF = 58

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
        elseif bs <= SMALL_FILT_CUTOFF
            _small_filt_fir!(out, b, x, si, col)
        else
            _filt_fir!(out, b, x, si, col)
        end
    end
    return out
end

# Transposed direct form II
function _filt_iir!(out, b, a, x, si, col)
    silen = length(si)
    @inbounds for i=1:size(x, 1)
        xi = x[i,col]
        val = muladd(xi, b[1], si[1])
        for j=1:(silen-1)
            si[j] = muladd(val, -a[j+1], muladd(xi, b[j+1], si[j+1]))
        end
        si[silen] = muladd(xi, b[silen+1], -a[silen+1]*val)
        out[i,col] = val
    end
end

# Transposed direct form II
function _filt_fir!(out, b, x, si, col)
    silen = length(si)
    @inbounds for i=1:size(x, 1)
        xi = x[i,col]
        val = muladd(xi, b[1], si[1])
        for j=1:(silen-1)
            si[j] = muladd(xi, b[j+1], si[j+1])
        end
        si[silen] = b[silen+1]*xi
        out[i,col] = val
    end
end

#
# filt implementation for FIR filters (faster than Base)
#

for n = 2:SMALL_FILT_CUTOFF
    silen = n-1
    si = [Symbol("si$i") for i = 1:silen]
    # Transposed direct form II
    @eval function _filt_fir!(out, b::NTuple{$n,T}, x, siarr, col) where {T}
        offset = (col - 1) * size(x, 1)

        $(Expr(:block, [:(@inbounds $(si[i]) = siarr[$i]) for i = 1:silen]...))
        @inbounds for i=1:size(x, 1)
            xi = x[i+offset]
            val = muladd(xi, b[1], $(si[1]))
            $(Expr(:block, [:($(si[j]) = muladd(xi, b[$(j+1)], $(si[j+1]))) for j = 1:(silen-1)]...))
            $(si[silen]) = b[$(silen+1)]*xi
            out[i+offset] = val
        end
    end
end

# Convert array filter tap input to tuple for small-filtering
let chain = :(throw(ArgumentError("invalid tuple size")))
    for n = SMALL_FILT_CUTOFF:-1:2
        chain = quote
            if length(h) == $n
                _filt_fir!(
                    out,
                    ($([:(@inbounds(h[$i])) for i = 1:n]...),),
                    x,
                    si,
                    col
                )
            else
                $chain
            end
        end
    end

    @eval function _small_filt_fir!(
        out::AbstractArray, h::AbstractVector{T}, x::AbstractArray, si, col
    ) where T
        $chain
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
    fill!(padded, zero(eltype(padded)))
    pad_data_ranges = UnitRange.(1, size(u))
    copyto!(padded, CartesianIndices(pad_data_ranges), u, CartesianIndices(u))

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

function _conv_clip!(
    y::AbstractVector,
    minpad,
    ::NTuple{<:Any, Base.OneTo{Int}},
    ::NTuple{<:Any, Base.OneTo{Int}}
)
    sizehint!(resize!(y, minpad[1]), minpad[1])
end
function _conv_clip!(
    y::AbstractArray,
    minpad,
    ::NTuple{<:Any, Base.OneTo{Int}},
    ::NTuple{<:Any, Base.OneTo{Int}}
)
    y[CartesianIndices(minpad)]
end
# For arrays with weird offsets
function _conv_clip!(y::AbstractArray, minpad, axesu, axesv)
    out_offsets = first.(axesu) .+ first.(axesv)
    out_axes = UnitRange.(out_offsets, out_offsets .+ minpad .- 1)
    out = similar(y, out_axes)
    copyto!(out, CartesianIndices(out), y, CartesianIndices(UnitRange.(1, minpad)))
end

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
    _conv_clip!(y, minpad, axes(u), axes(v))
end
function conv(u::AbstractArray{<:BLAS.BlasFloat, N},
              v::AbstractArray{<:BLAS.BlasFloat, N}) where N
    fu, fv = promote(u, v)
    conv(fu, fv)
end
conv(u::AbstractArray{T, N}, v::AbstractArray{T, N}) where {T<:Number, N} =
    conv(float(u), float(v))
conv(u::AbstractArray{<:Integer, N}, v::AbstractArray{<:Integer, N}) where {N} =
    round.(Int, conv(float(u), float(v)))
function conv(u::AbstractArray{<:Number, N},
              v::AbstractArray{<:BLAS.BlasFloat, N}) where N
    conv(float(u), v)
end
function conv(u::AbstractArray{<:BLAS.BlasFloat, N},
              v::AbstractArray{<:Number, N}) where N
    conv(u, float(v))
end

function conv(A::AbstractArray, B::AbstractArray)
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

dsp_reverse(v, ::NTuple{<:Any, Base.OneTo{Int}}) = reverse(v, dims = 1)
function dsp_reverse(v, vaxes)
    vsize = length(v)
    reflected_start = - first(vaxes[1]) - vsize + 1
    reflected_axes = (reflected_start : reflected_start + vsize - 1,)
    out = similar(v, reflected_axes)
    copyto!(out, reflected_start, Iterators.reverse(v), 1, vsize)
end


"""
    xcorr(u,v; padmode = :longest)

Compute the cross-correlation of two vectors, by calculating the similarity
between `u` and `v` with various offsets of `v`. Delaying `u` relative to `v`
will shift the result to the right.

The size of the output depends on the padmode keyword argument: with padmode =
:none the length of the result will be length(u) + length(v) - 1, as with conv.
With padmode = :longest the shorter of the arguments will be padded so they are
equal length. This gives a result with length 2*max(length(u), length(v))-1,
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
        conv(u, dsp_reverse(conj(v), axes(v)))
    elseif padmode == :none
        conv(u, dsp_reverse(conj(v), axes(v)))
    else
        throw(ArgumentError("padmode keyword argument must be either :none or :longest"))
    end
end
