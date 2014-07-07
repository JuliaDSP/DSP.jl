# Zero phase digital filtering for Julia
# Created by Robert Luke (Robert.Luke@med.kuleuven.be)
# Initial state and filter in place code by Matt Bauman https://github.com/mbauman

module ZeroPhaseFiltering

export filtfilt


##############
#
# filtfilt
#
##############

# zero phase digital filtering by processing data in forward and reverse direction
function filtfilt(b::AbstractVector, a::AbstractVector, x::AbstractVector)

    zi = filt_stepstate(b, a)

    pad_length = 3 * (max(length(a), length(b)) - 1)

    x = [vec(2*x[1] - x[pad_length+1:-1:2]) , x, vec(2 * x[end] - x[end-1:-1:end-pad_length])]

    x = flipud(filt!(b, a, x, zi*x[1]))
    x = flipud(filt!(b, a, x, zi*x[1]))

    # Return to original size by removing padded length
    x[pad_length+1: end-pad_length]
end


##############
#
# Determine initial state from Matt Bauman
# https://github.com/mbauman
#
# This will be removed when added to base
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
    (I - A) \ B
 end


##############
#
# filt with initial conditions from Matt Bauman
# https://github.com/mbauman
#
# This will be removed when added to base
#
##############

import Base: DSP.filt

# in-place filtering; both the input and filter state are modified in-place
function filt!{T<:Number}(b::Union(AbstractVector{T}, T), a::Union(AbstractVector{T}, T),
                          x::AbstractVector{T}, si::AbstractVector{T}=zeros(T, max(length(a), length(b))-1))
    if isempty(b); error("b must be non-empty"); end
    if isempty(a); error("a must be non-empty"); end
    if a[1]==0; error("a[1] must be nonzero"); end

    as = length(a)
    bs = length(b)
    sz = max(as, bs)

    if sz == 1
        # Simple scaling without memory; quick exit
        return scale!(x, b[1]/a[1])
    end

    if bs<sz
        # Ensure b has at least as many elements as a
        newb = zeros(T,sz)
        newb[1:bs] = b
        b = newb
    end

    xs = size(x,1)
    silen = sz-1
    size(si) == (silen,) || error("the vector of initial conditions must have exactly max(length(a),length(b))-1 elements")

    if a[1] != 1
        # Normalize the coefficients such that a[1] == 1
        norml = a[1]
        a ./= norml
        b ./= norml
    end

    if as > 1
        if as<sz
            # Pad a to be the same length as b
            newa = zeros(T,sz)
            newa[1:as] = a
            a = newa
        end

        @inbounds begin
            for i=1:xs
                val = si[1] + b[1]*x[i]
                for j=1:(silen-1)
                    si[j] = si[j+1] + b[j+1]*x[i] - a[j+1]*val
                end
                si[silen] = b[silen+1]*x[i] - a[silen+1]*val
                x[i] = val
            end
        end
    else
        @inbounds begin
            for i=1:xs
                val = si[1] + b[1]*x[i]
                for j=1:(silen-1)
                    si[j] = si[j+1] + b[j+1]*x[i]
                end
                si[silen] = b[silen+1]*x[i]
                x[i] = val
            end
        end
    end
    return x
end


end
