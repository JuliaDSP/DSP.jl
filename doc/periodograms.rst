:mod:`Periodograms` - periodogram estimation
===========================================

.. function:: arraysplit(s, n, m)

    Split an array into arrays of length ``n`` with overlapping regions
    of length ``m``. Iterating or indexing the returned AbstractVector
    always yields the same Vector with different contents.

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

.. function:: spectrogram(s, n=div(length(s), 8), noverlap=div(n, 2); onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, window=nothing)

    Computes the STFT of a signal ``s`` based on ``n`` segments
    with overlap ``noverlap`` and returns a matrix containing the STFT coefficients. See :func:`periodogram` for description of optional keyword arguments.


.. function:: periodogram(s::AbstractMatrix; nfft=nextfastfft(size(s)), fs=1, radialsum=false, radialavg=false)

    Computes periodogram of a 2-d signal by FFT and returns a
    Periodogram2 object.

    Returns a 2-d periodogram by default. A radially summed or 
    averaged periodogram is returned as a Periodogram object 
    if ``radialsum`` or  ``radialavg`` is true, respectively.

    ``nfft`` specifies the number of points to use for the Fourier
    transform. If ``size(s)`` < ``nfft``, then the input is padded
    with zeros. By default, ``nfft`` is the closest size for which the
    Fourier transform can be computed with maximal efficiency. ``fs`` 
    is the sample rate of the original signal in both directions.
    
    For ``radialsum=true`` the value of ``power[k]`` is proportional to
    :math:`\frac{1}{N}\sum_{k\leq |k'|<k+1} |X[k']|^2`.
    For ``radialavg=true`` it is proportional to
    :math:`\frac{1}{N \#\{k\leq |k'|<k+1\}} \sum_{k\leq |k'|<k+1} |X[k']|^2`.
    The computation of ``|k'|`` takes into account non-square signals
    by scaling the coordinates of the wavevector accordingly.

.. function:: freq(p)

	Returns the frequency bin centers for a given Periodogram or
	Spectrogram object.
	
	Returns a tuple of frequency bin centers for a given Periodogram2 
	object.

	See also: :func:`fftfreq`, :func:`rfftfreq`

.. function:: power(p)

    For a Periodogram, returns the computed power at each frequency as
    a Vector.

    For a Spectrogram, returns the computed power at each frequency and
    time bin as a Matrix. Dimensions are frequency × time.

.. function:: time(p)

    Returns the time bin centers for a given Spectrogram object.
 
