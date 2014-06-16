:mod:`Periodogram` - periodogram estimation
===========================================

.. function:: periodogram(s, [window::Function])

    Compute periodogram of signal ``s`` by using the Fast Fourier Transform. The optional ``window`` argument defines a window function that accept an integer input argument ``n`` and returns an array of length ``n`` or a single scalar value that will multiplied the signal ``s``. This argument defaults to a scalar value of 1, which means no windowing is performed by default. The power spectral density of the signal will be normalised to take into account the window gain.

.. function:: welch_pgram(s, n, m, [window::Function])

    Compute Welch periodogram of signal ``s`` by averaging the power spectral density on ``n`` segments with overlap ``m``. The optional ``window`` argument defines a window function as for the ``periodogram`` function, that is applied on each data segment before the transform in Fourier domain is taken.

.. function:: bartlett_pgram(s, n, [window::Function])

    Compute Bartlett periodogram of signal ``s` based on averaging of the periodograms of segments of data of length ``n``. This is equivalent to ``welch_pgram(s, n, 0)``, 
    so that no overlap between the segment is used. Optional ``window`` argument defines a window function applied on each segment of data before the transform in Fourier domain is taken.

.. function:: spectrogram(s; n=length(s)/8, m=n/2, r=1, window::Function=n->one(eltype(s)))

    Compute Spectrogram of signal ``s`` based on ``n`` segments with overlap ``m``, sampling rate ``r``, and using the window function ``window``.

.. function:: arraysplit(s, n::Integer, m::Integer)

    Split an array into arrays of length n with overlapping regions of length m.