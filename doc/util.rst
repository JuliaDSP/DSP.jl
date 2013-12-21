:mod:`Util` - utility functions
=================================

.. function:: unwrap(m::Array[, dim::Integer]; range::FloatingPoint=2pi)

    Assumes m (along dimension dim) to be a sequences of values that have been
    wrapped to be inside the given range, and undoes the wrapping by
    identifying discontinuities. If dim is not given, the last dimension is
    used.

    A common usage is for a phase measurement over time, such as when comparing
    successive frames of a short-time-fourier-transform, as each frame is
    wrapped to stay within (-pi, pi].

.. function:: unwrap!(m::Array[, dim::Integer]; range::FloatingPoint=2pi)

    In-place version of unwrap(v, dim, range)
