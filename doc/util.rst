:mod:`Util` - utility functions
=================================

.. function:: unwrap(m, dim=ndims(m); range=2pi)

    Assumes m (along dimension dim) to be a sequences of values that
    have been wrapped to be inside the given range (centered around
    zero), and undoes the wrapping by identifying discontinuities. If
    dim is not given, the last dimension is used.

    A common usage is for a phase measurement over time, such as when
    comparing successive frames of a short-time-fourier-transform, as
    each frame is wrapped to stay within (-pi, pi].

.. function:: unwrap!(m, dim=ndims(m); range=2pi)

    In-place version of unwrap(m, dim, range)

.. function:: hilbert(x)

    Computes the analytic representation of x, :math:`x_a = x + j
    \hat{x}`, where :math:`\hat{x}` is the Hilbert transform of x.

.. function:: fftfreq(n, fs=1)

    Return discrete fourier transform sample frequencies. The returned
    Frequencies object is an AbstractVector containing the frequency
    bin centers at every sample point. ``fs`` is the sample rate of the
    input signal.

.. function:: rfftfreq(n, fs=1)

    Return discrete fourier transform sample frequencies for use with
    ``rfft`. The returned Frequencies object is an AbstractVector
    containing the frequency bin centers at every sample point. ``fs``
    is the sample rate of the input signal.

.. function:: nextfastfft(n)

    Return the closest product of 2, 3, 5, and 7 greater than or equal
    to ``n``. FFTW contains optimized kernels for these sizes and
    computes Fourier transforms of input that is a product of these
    sizes faster than for input of other sizes.
