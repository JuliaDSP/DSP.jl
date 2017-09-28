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
    \hat{x}`, where :math:`\hat{x}` is the Hilbert transform of x,
    along the first dimension of x.

.. function:: fftfreq(n, fs=1)

    Return discrete fourier transform sample frequencies. The returned
    Frequencies object is an AbstractVector containing the frequency
    bin centers at every sample point. ``fs`` is the sample rate of the
    input signal.

.. function:: rfftfreq(n, fs=1)

    Return discrete fourier transform sample frequencies for use with
    ``rfft``. The returned Frequencies object is an AbstractVector
    containing the frequency bin centers at every sample point. ``fs``
    is the sample rate of the input signal.

.. function:: nextfastfft(n)

    Return the closest product of 2, 3, 5, and 7 greater than or equal
    to ``n``. FFTW contains optimized kernels for these sizes and
    computes Fourier transforms of input that is a product of these
    sizes faster than for input of other sizes.

.. function:: pow2db(a)

    Convert a power ratio to dB (decibel), or :math:`10\log_{10}(a)`. 
    The inverse of ``db2pow``.

.. function:: amp2db(a)

    Convert an amplitude ratio to dB (decibel), or :math:`20
    \log_{10}(a)=10\log_{10}(a^2)`. The inverse of ``db2amp``.

.. function:: db2pow(a)

    Convert dB to a power ratio. This function call also be called 
    using ``a*dB``, i.e. ``3dB == db2pow(3)``. The inverse of ``pow2db``.

.. function:: db2amp(a)

    Convert dB to an amplitude ratio. This function call also be called 
    using ``a*dBa``, i.e. ``3dBa == db2amp(3)``. The inverse of ``amp2db``.

.. function:: rms(s)

    Return the root mean square of signal ``s``.

.. function:: rmsfft(f)

    Return the root mean square of signal ``s`` given the FFT transform
    ``f = fft(s)``. Equivalent to ``rms(ifft(f))``.

.. function:: meanfreq(f, fs)

    Calculate the mean power frequency of `f` with a sampling frequency of `fs`, defined as:

    .. math:: \\textup{MPF} = \\frac{\\sum_{i=1}^{F} f_i X_i^2 }{\\sum_{i=0}^{F} X_i^2 } \\textup{Hz}

    where `F` is the Nyquist frequency, and `X` is the power spectral density.

