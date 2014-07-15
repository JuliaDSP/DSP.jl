# This file contains functions included in recent but not earlier julia versions


module BackwardsCompatability

export filt!


##############
#
# filt with initial conditions from Matt Bauman
# https://github.com/mbauman
#
# https://github.com/JuliaLang/julia/pull/7560
#
##############

import Base: DSP.filt
import Base.trailingsize



_zerosi(b,a,T) = zeros(promote_type(eltype(b), eltype(a), T), max(length(a), length(b))-1)

function filt{T,S}(b::Union(AbstractVector, Number), a::Union(AbstractVector, Number),
                   x::AbstractArray{T}, si::AbstractArray{S}=_zerosi(b,a,T))
    filt!(Array(promote_type(eltype(b), eltype(a), T, S), size(x)), b, a, x, si)
end

# in-place filtering: returns results in the out argument, which may shadow x
# (and does so by default)
function filt!{T,S,N}(out::AbstractArray, b::Union(AbstractVector, Number), a::Union(AbstractVector, Number),
                      x::AbstractArray{T}, si::AbstractArray{S,N}=_zerosi(b,a,T))
    isempty(b) && error("b must be non-empty")
    isempty(a) && error("a must be non-empty")
    a[1] == 0 && error("a[1] must be nonzero")
    size(x) != size(out) && error("out size must match x")

    as = length(a)
    bs = length(b)
    sz = max(as, bs)
    silen = sz - 1
    xs = size(x,1)
    ncols = trailingsize(x,2)

    size(si, 1) != silen && error("si must have max(length(a),length(b))-1 rows")
    N > 1 && trailingsize(si,2) != ncols && error("si must either be a vector or have the same number of columns as x")

    xs == 0 && return out
    sz == 1 && return scale!(out, x, b[1]/a[1]) # Simple scaling without memory

    # Filter coefficient normalization
    if a[1] != 1
        norml = a[1]
        a ./= norml
        b ./= norml
    end
    # Pad the coefficients with zeros if needed
    bs<sz && (b = copy!(zeros(eltype(b), sz), b))
    1<as<sz && (a = copy!(zeros(eltype(a), sz), a))

    initial_si = si
    for col = 1:ncols
        # Reset the filter state
        si = initial_si[:, N > 1 ? col : 1]
        if as > 1
            @inbounds for i=1:xs
                xi = x[i,col]
                val = si[1] + b[1]*xi
                for j=1:(silen-1)
                    si[j] = si[j+1] + b[j+1]*xi - a[j+1]*val
                end
                si[silen] = b[silen+1]*xi - a[silen+1]*val
                out[i,col] = val
            end
        else
            @inbounds for i=1:xs
                xi = x[i,col]
                val = si[1] + b[1]*xi
                for j=1:(silen-1)
                    si[j] = si[j+1] + b[j+1]*xi
                end
                si[silen] = b[silen+1]*xi
                out[i,col] = val
            end
        end
    end
    return out
end





end  # module
