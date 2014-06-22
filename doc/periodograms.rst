:mod:`Periodograms` - periodogram estimation
===========================================

.. function:: arraysplit(s, n::Integer, m::Integer)

    Split an array into arrays of length n with overlapping regions of length m.

.. function:: periodogram(s; onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, window=nothing)

    Computes periodogram of a signal by FFT and returns a
    Periodogram object.

    For real signals, the two-sided periodogram is symmetric and this
    function returns a one-sided (real only) periodogram by default. A
    two-sided periodogram can be obtained by setting ``onesided=true``.

    ``nfft`` specifies the number of points to use for the Fourier
    transform. If ``length(s)`` < ``nfft``, then the input is padded
    with zeros. By default, ``nfft`` is the closest size for which the
    Fourier transform can be computed with maximal efficiency.

    ``fs`` is the sample rate of the original signal, and ``window`` is
    an optional window function or vector to be applied to the original
    signal before computing the Fourier transform. The computed
    periodogram is normalized so that the area under the periodogram is
    equal to the uncentered variance of the original signal.

.. function:: welch_pgram(s, n=div(length(s), 8), noverlap=div(n, 2); onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, window=nothing)

    Computes the Welch periodogram of a signal ``s`` based on ``n``
    segments with overlap ``noverlap`` and returns a Periodogram
    object. For a Bartlett periodogram, set ``noverlap=0``. See
    :func:`periodogram` for description of optional keyword arguments.

.. function:: spectrogram(s, n=div(length(s), 8), noverlap=div(n, 2); onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, window=nothing)

    Computes the spectrogram of a signal ``s`` based on ``n`` segments
    with overlap ``noverlap`` and returns a Spectrogram object. See
    :func:`periodogram` for description of optional keyword arguments.

.. function:: freq(p)

	Returns the frequency bin centers for a given Periodogram or
	Spectrogram object.

	See also: :func:`fftfreq`, :func:`rfftfreq`

.. function:: power(p)

    For a Periodogram, returns the computed power at each frequency as
    a Vector.

    For a Spectrogram, returns the computed power at each frequency and
    time bin as a Matrix. Dimensions are frequency × time.

.. function:: time(p)

    Returns the time bin centers for a given Spectrogram object.
 