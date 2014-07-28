:mod:`Filters` - filter design and filtering
============================================

Linear time-invariant filter representations
--------------------------------------------

DSP.jl supports common filter representations. Filters can be converted
from one type to another using ``convert``.

.. function:: ZPKFilter(z, p, k)

    Filter representation in terms of zeros (``z``), poles(``p``), and
    gain (``k``):

    .. math:: H(z) = k\frac{(s - \verb!z[1]!) \ldots (s - \verb!z[end]!)}{(s - \verb!p[1]!) \ldots (s - \verb!p[end]!)}

.. function:: TFFilter(b, a)

    Filter representation in terms of the coefficients of the numerator
    (``b``) and denominator ``a`` of the transfer function:

    .. math:: H(z) = \frac{\verb!b[1]! s^{n-1} + \ldots + \verb!b[n]!}{\verb!a[1]! s^{n-1} + \ldots + \verb!a[n]!}

    or equivalently:

    .. math:: H(z) = \frac{\verb!b[1]! + \ldots + \verb!b[n]! s^{-n+1}}{\verb!a[1]! + \ldots + \verb!b[n]! s^{-n+1}}

    ``b`` and ``a`` may be specified as ``Polynomial`` objects or
    vectors ordered from highest power to lowest.

.. function:: BiquadFilter(b0, b1, b2, a1, a2)

    Filter representation in terms of the transfer function of a single
    second-order section given by:

    .. math:: H(z) = \frac{\verb!b0! s^2+\verb!b1! s+\verb!b2!}{s^2+\verb!a1! s + \verb!a2!}

    or equivalently:

    .. math:: H(z) = \frac{\verb!b0!+\verb!b1! z^{-1}+\verb!b2! z^{-2}}{1+\verb!a1! z^{-1} + \verb!a2! z^{-2}}

.. function:: SOSFilter(biquads, gain)

    Filter representation in terms of a cascade of second-order
    sections and gain. ``biquads`` must be specified as a vector of
    ``BiquadFilters``.

Filter application
------------------

.. function:: filt(f, x[, si])

    Apply filter ``f`` along the first dimension of array ``x``. ``si``
    is an optional array representing the initial filter state
    (defaults to zeros). If ``f`` is a ``TFFilter``, ``BiquadFilter``,
    or ``SOSFilter``, filtering is implemented directly. If ``f`` is a
    ``ZPKFilter``, it is first converted to an ``SOSFilter``.

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

.. function:: analogfilter(responsetype, prototype)

    Construct an analog filter.

.. function:: digitalfilter(responsetype, prototype)

    Construct a digital filter.

Filter response types
---------------------

.. function:: Lowpass(Wn[; fs])

    Low pass filter with cutoff frequency ``Wn``. If ``fs`` is not
    specified, ``Wn`` is interpreted as a normalized frequency in
    half-cycles/sample.

.. function:: Highpass(Wn[; fs])

    High pass filter with cutoff frequency ``Wn``. If ``fs`` is not
    specified, ``Wn`` is interpreted as a normalized frequency in
    half-cycles/sample.

.. function:: Bandpass(Wn1, Wn2[; fs])

    Band pass filter with normalized pass band (``Wn1``, ``Wn2``). If
    ``fs`` is not specified, ``Wn`` is interpreted as a normalized
    frequency in half-cycles/sample.

.. function:: Bandstop(Wn1, Wn2[; fs])

    Band stop filter with normalized stop band (``Wn1``, ``Wn2``). If
    ``fs`` is not specified, ``Wn`` is interpreted as a normalized
    frequency in half-cycles/sample.


Filter prototypes
-----------------

.. function:: Butterworth(n)

    ``n`` pole Butterworth filter.

.. function:: Chebyshev1(n, ripple)

    ``n`` pole Chebyshev type I filter with ``ripple`` dB ripple in
    the passband.

.. function:: Chebyshev2(n, ripple)

    ``n`` pole Chebyshev type II filter with ``ripple`` dB ripple in
    the stopband.

.. function:: Elliptic(n, rp, rs) 

    ``n`` pole elliptic (Cauer) filter with ``rp`` dB ripple in the
    passband and ``rs`` dB attentuation in the stopband.


Filter response
---------------

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


Miscellaneous
-------------

 .. function:: coefb(f)

     Coefficients of the numerator of a TFFilter object, highest power
     first, i.e., the ``b`` passed to ``Base.filt()``

 .. function:: coefa(f)

     Coefficients of the denominator of a TFFilter object, highest power
     first, i.e., the ``a`` passed to ``Base.filt()``


Examples
--------

Construct a 4th order elliptic lowpass filter with normalized cutoff
frequency 0.2, 0.5 dB of passband ripple, and 30 dB attentuation in
the stopband and extract the coefficients of the numerator and
denominator of the transfer function::

  responsetype = Lowpass(0.2)
  prototype = Elliptic(4, 0.5, 30)
  tf = convert(TFFilter, digitalfilter(responsetype, prototype))
  numerator_coefs = coefb(tf)
  denominator_ceofs = coefa(tf)

Filter the data in ``x``, sampled at 1000 Hz, with a 4th order
Butterworth bandpass filter between 10 and 40 Hz::

  responsetype = Bandpass(10, 40; fs=1000)
  prototype = Butterworth(4)
  filt(digitalfilter(responsetype, prototype), x)
