# Zero phase digital filtering for Julia
# Created by Robert Luke (Robert.Luke@med.kuleuven.be)
# Initial state and filter in place code by Matt Bauman https://github.com/mbauman

module ZeroPhaseFiltering

using ..BackwardsCompatability

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

    x = flipud(filt!(x, b, a, x, zi*x[1]))
    x = flipud(filt!(x, b, a, x, zi*x[1]))

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
    scale_factor \ (eye(size(A)[1]) - A) \ B
 end




end
