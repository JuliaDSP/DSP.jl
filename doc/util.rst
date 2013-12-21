:mod:`Util` - utility functions
=================================

.. function:: unwrap{T <: FloatingPoint}(m::Array{T}, dim::Integer=ndims(m);
                                         range::Number=2pi)

    Assumes m (along dimension dim) to be a sequences of values that have been
    wrapped to be inside the given range, and undoes the wrapping by
    identifying discontinuities. If dim is not given, the last dimension is
    used.

    A common usage is for a phase measurement over time, such as when comparing
    successive frames of a short-time-fourier-transform, as each frame is
    wrapped to stay within (-pi, pi].

.. function:: unwrap!{T <: FloatingPoint}(m::Array{T}, dim::Integer=ndims(m);
                                          range::Number=2pi)

    In-place version of unwrap(m, dim, range)
