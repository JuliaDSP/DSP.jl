# This file was formerly a part of Julia. License is MIT: https://julialang.org/license

const SMALL_FILT_CUTOFF = 66

"""
    filt(b::Union{AbstractVector,Number},
         a::Union{AbstractVector,Number},
         x::AbstractArray)

Apply filter described by vectors `a` and `b` to vector `x`.

Inputs that are `Number`s are treated as one-element `Vector`s.
"""
filt(b::Union{AbstractVector, Number}, a::Union{AbstractVector, Number}, x::AbstractArray{T}) where {T} =
    filt!(similar(x, promote_type(eltype(b), eltype(a), T)), b, a, x)

# in-place filtering: returns results in the out argument, which may shadow x
# (and does so by default)

"""
    filt!(out, b, a, x)

Same as [`filt`](@ref) but writes the result into the `out` argument, which may
alias the input `x` to modify it in-place.
"""
function filt!(out::AbstractArray, b::Union{AbstractVector, Number}, a::Union{AbstractVector, Number},
               x::AbstractArray{T}) where {T}
    isempty(b) && throw(ArgumentError("filter vector b must be non-empty"))
    isempty(a) && throw(ArgumentError("filter vector a must be non-empty"))
    a[1] == 0  && throw(ArgumentError("filter vector a[1] must be nonzero"))
    if size(x) != size(out)
        throw(ArgumentError(lazy"output size $(size(out)) must match input size $(size(x))"))
    end

    as = length(a)
    bs = length(b)
    sz = max(as, bs)

    iszero(size(x, 1)) && return out
    isone(sz) && return (k = b[1] / a[1]; @noinline mul!(out, x, k)) # Simple scaling without memory

    # Filter coefficient normalization
    if !isone(a[1])
        norml = a[1]
        a = @noinline broadcast(/, a, norml)
        b = @noinline broadcast(/, b, norml)
    end

    si = Vector{promote_type(eltype(b), eltype(a), T)}(undef, sz - 1)

    if as == 1 && bs <= SMALL_FILT_CUTOFF
        fill!(si, zero(eltype(si)))
        _small_filt_fir!(out, b, x, si, Val(bs))
    else
        for col in CartesianIndices(axes(x)[2:end])
            # Reset the filter state
            fill!(si, zero(eltype(si)))
            if as > 1
                _filt_iir!(out, b, a, x, si, col)
            else
                _filt_fir!(out, b, x, si, col)
            end
        end
    end
    return out
end

# Transposed direct form II
function _filt_iir!(out, b, a, x, si, col)
    silen = length(si)
    @inbounds for i in axes(x, 1)
        xi = x[i, col]
        val = muladd(xi, b[1], si[1])
        out[i, col] = val
        for j in 1:min(length(a), length(b), silen) - 1
            si[j] = muladd(val, -a[j+1], muladd(xi, b[j+1], si[j+1]))
        end
        if length(a) == length(b)
            si[silen] = muladd(xi, b[silen+1], -a[silen+1]*val)
        elseif length(a) > length(b)
            for j in length(b):silen-1
                si[j] = muladd(val, -a[j+1], si[j+1])
            end
            si[silen] = -a[silen+1]*val
        else
            for j in length(a):silen-1
                si[j] = muladd(xi, b[j+1], si[j+1])
            end
            si[silen] = xi*b[silen+1]
        end
    end
end

# Transposed direct form II
function _filt_fir!(out, b, x, si, col)
    silen = length(si)
    @inbounds for i in axes(x, 1)
        xi = x[i, col]
        out[i, col] = muladd(xi, b[1], si[1])
        for j=1:(silen-1)
            si[j] = muladd(xi, b[j+1], si[j+1])
        end
        si[silen] = b[silen+1] * xi
    end
end

#
# filt implementation for FIR filters
#

### NOTE ###
# Fragile because of the impact of @inbounds and checkbounds
# on the effects system

const SMALL_FILT_VECT_CUTOFF = 19

# Transposed direct form II
@generated function _filt_fir!(out, b::NTuple{N,T}, x, siarr, col, ::Val{StoreSI}=Val(false)) where {N,T,StoreSI}
    silen = N - 1
    si_end = Symbol(:si_, silen)

    quote
        N <= SMALL_FILT_VECT_CUTOFF && checkbounds(siarr, $silen)
        Base.@nextract $silen si siarr
        for i in axes(x, 1)
            xi = x[i, col]
            val = muladd(xi, b[1], si_1)
            Base.@nexprs $(silen - 1) j -> (si_j = muladd(xi, b[j+1], si_{j + 1}))
            $si_end = xi * b[N]
            if N > SMALL_FILT_VECT_CUTOFF
                @inbounds out[i, col] = val
            else
                out[i, col] = val
            end
        end
        if StoreSI
            Base.@nexprs $silen j -> siarr[j] = si_j
        end
        return nothing
    end
end

# Convert array filter tap input to tuple for small-filtering
function _small_filt_fir!(
    out::AbstractArray, h::AbstractVector, x::AbstractArray,
        si::AbstractVector, ::Val{bs}) where {bs}

    bs < 2 && throw(ArgumentError("invalid tuple size"))
    length(h) != bs && throw(ArgumentError("length(h) does not match bs"))
    b = ntuple(j -> h[j], Val(bs))
    for col in CartesianIndices(axes(x)[2:end])
        _filt_fir!(out, b, x, si, col)
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
