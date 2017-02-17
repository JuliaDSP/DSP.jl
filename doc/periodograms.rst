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
    two-sided periodogram can be obtained by setting ``onesided=false``.

    ``nfft`` specifies the number of points to use for the Fourier
    transform. If ``length(s)`` < ``nfft``, then the input is padded
    with zeros. By default, ``nfft`` is the closest size for which the
    Fourier transform can be computed with maximal efficiency.

    ``fs`` is the sample rate of the original signal, and ``window`` is
    an optional window function or vector to be applied to the original
    signal before computing the Fourier transform. The computed
    periodogram is normalized so that the area under the periodogram is
    equal to the uncentered variance (or average power) of the original
    signal.

.. function:: welch_pgram(s, n=div(length(s), 8), noverlap=div(n, 2); onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, window=hanning)

    Computes the Welch periodogram of a signal ``s`` by dividing the signal
    into overlapping segments of length ``n``, averaging the periodograms of
    each segment, and returning a ``Periodogram`` object.  ``noverlap`` is the
    number of samples of overlap, with a new segment starting every
    ``n-noverlap`` samples of ``s``.  Each segment is padded with zeros to a
    length of ``nfft`` for efficiency.

    By default, the signal is windowed using the hanning window as a moderate
    tradeoff between high spectral resolution (the square window) and
    sensitivity to signals near the noise floor.  A window which falls to zero
    also reduces the spurious ringing which occurs when ``n!=nfft`` for signals
    with large mean.

    Compared to ``periodogram()``, ``welch_pgram()`` reduces noise in the
    estimated power spectral density as ``length(s)`` increases. In the basic
    periodogram, we increase the frequency resolution instead.  For a Bartlett
    periodogram, set ``noverlap=0``.

    See :func:`periodogram` for a description of the other keyword arguments.

.. function:: mt_pgram(s; onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, nw=4, ntapers=iceil(2nw)-1, window=dpss(length(s), nw, ntapers))

    Computes the multitaper periodogram of a signal ``s``.

    If ``window`` is not specified, the signal is tapered with
    ``ntapers`` discrete prolate spheroidal sequences with
    time-bandwidth product ``nw``. Each sequence is equally weighted;
    adaptive multitaper is not (yet) supported.

    If ``window`` is specified, each column is applied as a taper. The
    sum of periodograms is normalized by the total sum of squares of
    ``window``.

    See also: :func:`dpss`

.. function:: spectrogram(s, n=div(length(s), 8), noverlap=div(n, 2); onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, window=nothing)

    Computes the spectrogram of a signal ``s`` based on segments with ``n`` samples
    with overlap of ``noverlap`` samples, and returns a Spectrogram object. See
    :func:`periodogram` for description of optional keyword arguments.

.. function:: stft(s, n=div(length(s), 8), noverlap=div(n, 2); onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, window=nothing)

    Computes the STFT of a signal ``s`` based on segments with ``n`` samples
    with overlap of ``noverlap`` samples, and returns a matrix containing the STFT
    coefficients. See :func:`periodogram` for description of optional
    keyword arguments.

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
 
