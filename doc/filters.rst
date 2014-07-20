:mod:`Filters` - filter design and filtering
============================================

Linear time-invariant filter representations
--------------------------------------------

DSP.jl supports common filter representations. Filters can be converted
from one type to another using ``convert``.

.. function:: ZPKFilter(z, p, k)

    Filter representation in terms of zeros (``z``), poles(``p``), and
    gain (``k``).

.. function:: TFFilter(b, a)

    Filter representation in terms of the coefficients of the numerator
    (``b``) and denominator ``a`` of the transfer function. ``b`` and
    ``a`` may be either ``Polynomial`` objects or vectors ordered from
    highest power to lowest.

.. function:: BiquadFilter(b0, b1, b2, a1, a2)

    Filter representation in terms of the transfer function of a single
    second-order section given by:

    .. math:: \frac{\verb!b0! s^2+\verb!b1! s+\verb!b2!}{s^2+\verb!a1! s + \verb!a2!}

    or equivalently:

    .. math:: \frac{\verb!b0!+\verb!b1! z^{-1}+\verb!b2! z^{-2}}{1+\verb!a1! z^{-1} + \verb!a2! z^{-2}}

.. function:: SOSFilter(biquads, gain)

    Filter representation in terms of a cascade of second-order
    sections and gain. ``biquads`` must be specified as a vector of
    ``BiquadFilters``.

Filter application
------------------

.. function:: filt(f, x[, si])

    Apply filter ``f`` along the first dimension of array ``x``. ``si``
    is an optional array representing the initial filter state
    (defaults to zeros).

.. function:: filt!(out, f, x[, si])

    Same as :func:`filt()` but writes the result into the ``out``
    argument, which may alias the input ``x`` to modify it in-place.

.. function:: filtfilt(f, x)

    Filter ``x`` in the forward and reverse directions. The initial
    state of the filter is computed so that its response to a step
    function is steady state. Before filtering, the data is
    extrapolated at both ends with an odd-symmetric extension of length
    ``3*(max(length(b), length(a))-1)``.

    Because ``filtfilt`` applies the given filter twice, the effective
    filter order is twice the order of ``f``. The resulting signal has
    zero phase distortion.

.. function:: fftfilt(b, x)

    Perform overlap-save filtering of ``x`` using filter ``b``.

.. function:: firfilt(b, x)

    Filter ``x`` using filter ``b``, using :func:`filt` or
    :func:`fftfilt` depending on the lengths of ``b`` and ``x``.

Filter design
-------------

.. function:: analogfilter(responsetype, family)

    Construct an analog filter.

.. function:: digitalfilter(responsetype, family)

    Construct a digital filter.

Filter response types
---------------------

.. function:: Lowpass(Wn)

    Low pass filter with normalized cutoff frequency ``Wn``.

.. function:: Highpass(Wn)

    High pass filter with normalized cutoff frequency ``Wn``.

.. function:: Bandpass(Wn1, Wn2)

    Band pass filter with normalized pass band (``Wn1``, ``Wn2``).

.. function:: Bandstop(Wn1, Wn2)

    Band stop filter with normalized stop band (``Wn1``, ``Wn2``).


Filter families
---------------

.. function:: Butterworth(n) 

    ``n`` pole Butterworth filter.

Filter response
-----------------------

.. function:: freqz(filter, w)

    Frequency response of a digital ``filter`` at normalised frequency
    or frequencies ``w`` in radians/sample.

.. function:: freqz(filter, hz, fs)

    Frequency response of a digital ``filter`` at frequency or
    frequencies ``hz`` with sampling rate ``fs``.

.. function:: freqs(filter, w)

    Frequency response of an analog ``filter`` at normalised frequency
    or frequencies ``w`` in radians/sample.

.. function:: freqs(filter, hz, fs)

    Frequency response of an analog ``filter`` at frequency or
    frequencies ``hz`` with sampling rate ``fs``.

Examples
--------

TODO
