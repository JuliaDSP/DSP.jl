:mod:`Filters` - filter design and filtering
============================================

DSP.jl differentiates between filter coefficients and filters. Filter
coefficient objects specify the response of the filter in one of
several standard forms. Filter objects carry the state of the filter
together with filter coefficients in an implementable form
(``PolynomialRatio``, ``Biquad``, or ``SecondOrderSections``).
When invoked on a filter coefficient object, ``filt`` preserves state
between invocations.

Linear time-invariant filter coefficient objects
------------------------------------------------

DSP.jl supports common filter representations. Filter coefficients can
be converted from one type to another using ``convert``.

.. function:: ZeroPoleGain(z, p, k)

    Filter representation in terms of zeros ``z``, poles ``p``, and
    gain ``k``:

    .. math:: H(x) = k\frac{(x - \verb!z[1]!) \ldots (x - \verb!z[end]!)}{(x - \verb!p[1]!) \ldots (x - \verb!p[end]!)}

.. function:: PolynomialRatio(b, a)

    Filter representation in terms of the coefficients of the numerator
    ``b`` and denominator ``a`` of the transfer function:

    .. math:: H(s) = \frac{\verb!b[1]! s^{n-1} + \ldots + \verb!b[n]!}{\verb!a[1]! s^{n-1} + \ldots + \verb!a[n]!}

    or equivalently:

    .. math:: H(z) = \frac{\verb!b[1]! + \ldots + \verb!b[n]! z^{-n+1}}{\verb!a[1]! + \ldots + \verb!a[n]! z^{-n+1}}

    ``b`` and ``a`` may be specified as ``Polynomial`` objects or
    vectors ordered from highest power to lowest.

.. function:: Biquad(b0, b1, b2, a1, a2)

    Filter representation in terms of the transfer function of a single
    second-order section given by:

    .. math:: H(s) = \frac{\verb!b0! s^2+\verb!b1! s+\verb!b2!}{s^2+\verb!a1! s + \verb!a2!}

    or equivalently:

    .. math:: H(z) = \frac{\verb!b0!+\verb!b1! z^{-1}+\verb!b2! z^{-2}}{1+\verb!a1! z^{-1} + \verb!a2! z^{-2}}

.. function:: SecondOrderSections(biquads, gain)

    Filter representation in terms of a cascade of second-order
    sections and gain. ``biquads`` must be specified as a vector of
    ``Biquads``.

Filter objects
----------------------------------

.. function:: DF2TFilter(coef[, si])

    Construct a stateful direct form II transposed filter with
    coefficients ``coef``. ``si`` is an optional array representing the
    initial filter state (defaults to zeros). If ``f`` is a
    ``PolynomialRatio``, ``Biquad``, or ``SecondOrderSections``,
    filtering is implemented directly. If ``f`` is a ``ZeroPoleGain``
    object, it is first converted to a ``SecondOrderSections`` object.

DSP.jl's ``FIRFilter`` type maintains state between calls to :func:`filt`, allowing
you to filter a signal of indefinite length in RAM-friendly chunks. ``FIRFilter``
contains nothing more that the state of the filter, and a ``FIRKernel``. There are
five different kinds of ``FIRKernel`` for single rate, up-sampling, down-sampling,
rational resampling, and arbitrary sample-rate conversion. You need not specify the
type of kernel. The ``FIRFilter`` constructor selects the correct kernel based on input
parameters.

.. function:: FIRFilter(h[, ratio])

    Construct a stateful FIRFilter object from the vector of filter taps ``h``.
    ``ratio`` is an optional rational integer which specifies
    the input to output sample rate relationship (e.g. ``147//160`` for
    converting recorded audio from 48 KHz to 44.1 KHz).

.. function:: FIRFilter(h, rate[, Nϕ])

    Returns a polyphase FIRFilter object from the vector of filter taps ``h``.
    ``rate`` is a floating point number that specifies the input to output
    sample-rate relationship :math:`\frac{fs_{out}}{fs_{in}}`. ``Nϕ`` is an
    optional parameter which specifies the number of *phases* created from
    ``h``. ``Nϕ`` defaults to 32.

Filter application
------------------

.. function:: filt(f, x[, si])

    Apply filter or filter coefficients ``f`` along the first dimension
    of array ``x``. If ``f`` is a filter coefficient object, ``si``
    is an optional array representing the initial filter state (defaults
    to zeros). If ``f`` is a ``PolynomialRatio``, ``Biquad``, or
    ``SecondOrderSections``, filtering is implemented directly. If
    ``f`` is a ``ZeroPoleGain`` object, it is first converted to a
    ``SecondOrderSections`` object.  If ``f`` is a Vector, it is
    interpreted as an FIR filter, and a naïve or FFT-based algorithm is
    selected based on the data and filter length.

.. function:: filt!(out, f, x[, si])

    Same as :func:`filt()` but writes the result into the ``out``
    argument, which may alias the input ``x`` to modify it in-place.

.. function:: filtfilt(coef, x)

    Filter ``x`` in the forward and reverse directions using filter
    coefficients ``f``. The initial state of the filter is computed so
    that its response to a step function is steady state. Before
    filtering, the data is extrapolated at both ends with an
    odd-symmetric extension of length
    ``3*(max(length(b), length(a))-1)``.

    Because ``filtfilt`` applies the given filter twice, the effective
    filter order is twice the order of ``f``. The resulting signal has
    zero phase distortion.

.. function:: fftfilt(b, x)

    Apply FIR filter ``b`` along the first dimension of array ``x``
    using an FFT-based overlap-save algorithm.


Filter design
-------------

Most analog and digital filters are constructed by composing
:ref:`response-types`, which determine the frequency response of the
filter, with :ref:`design-methods`, which determine how the filter is
constructed.

.. function:: analogfilter(responsetype, designmethod)

    Construct an analog filter. See below for possible response and
    filter types.

.. function:: digitalfilter(responsetype, designmethod)

    Construct a digital filter. See below for possible response and
    filter types.

For some filters, the design method inherently implies a response type.
Such filters are documented below.

.. function:: iirnotch(Wn, bandwidth[; fs])

    Second-order digital IIR notch filter at frequency ``Wn`` with
    bandwidth ``bandwidth``. If ``fs`` is not specified, ``Wn`` is
    interpreted as a normalized frequency in half-cycles/sample.

.. _response-types:

Filter response types
~~~~~~~~~~~~~~~~~~~~~

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
    ``fs`` is not specified, ``Wn1`` and ``Wn2`` are interpreted as
    normalized frequencies in half-cycles/sample.

.. function:: Bandstop(Wn1, Wn2[; fs])

    Band stop filter with normalized stop band (``Wn1``, ``Wn2``). If
    ``fs`` is not specified, ``Wn1`` and ``Wn2`` are interpreted as
    normalized frequencies in half-cycles/sample.

.. design-methods:

Filter design methods
~~~~~~~~~~~~~~~~~~~~~

IIR filter design methods
:::::::::::::::::::::::::

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


FIR filter design methods
:::::::::::::::::::::::::

.. function:: FIRWindow(window; scale=true)

    FIR filter design using window ``window``, a vector whose length
    matches the number of taps in the resulting filter.

    If ``scale`` is ``true`` (default), the designed FIR filter is
    scaled so that the following holds:

    - For :func:`Lowpass` and :func:`Bandstop` filters, the frequency
      response is unity at 0 (DC).
    - For :func:`Highpass` filters, the frequency response is unity
      at the Nyquist frequency.
    - For :func:`Bandpass` filters, the frequency response is unity
      in the center of the passband.

.. function:: FIRWindow(; transitionwidth, attenuation=60, scale=true)

    Kaiser window FIR filter design. The required number of taps is
    calculated based on ``transitionwidth`` (in half-cycles/sample)
    and stopband ``attenuation`` (in dB). ``attenuation`` defaults to
    60 dB.


Filter response
---------------

.. function:: freqz(filter, w = linspace(0, π, 250))

    Frequency response of a digital ``filter`` at normalised frequency
    or frequencies ``w`` in radians/sample.

.. function:: freqz(filter, hz, fs)

    Frequency response of a digital ``filter`` at frequency or
    frequencies ``hz`` with sampling rate ``fs``.

.. function:: phasez(filter, w = linspace(0, π, 250))

    Phase response of a digital ``filter`` at normalised frequency
    or frequencies ``w`` in radians/sample.

.. function:: impz(filter, n=100)

    Impulse response of a digital ``filter`` with ``n`` points.

.. function:: stepz(filter, n=100)

    Step response of a digital ``filter`` with ``n`` points.

.. function:: freqs(filter, w)

    Frequency response of an analog ``filter`` at normalised frequency
    or frequencies ``w`` in radians/sample.

.. function:: freqs(filter, hz, fs)

    Frequency response of an analog ``filter`` at frequency or
    frequencies ``hz`` with sampling rate ``fs``.


Miscellaneous
-------------

.. function:: coefb(f)

    Coefficients of the numerator of a PolynomialRatio object, highest power
    first, i.e., the ``b`` passed to ``Base.filt()``

.. function:: coefa(f)

    Coefficients of the denominator of a PolynomialRatio object, highest power
    first, i.e., the ``a`` passed to ``Base.filt()``


Examples
--------

Construct a 4th order elliptic lowpass filter with normalized cutoff
frequency 0.2, 0.5 dB of passband ripple, and 30 dB attentuation in
the stopband and extract the coefficients of the numerator and
denominator of the transfer function::

  responsetype = Lowpass(0.2)
  designmethod = Elliptic(4, 0.5, 30)
  tf = convert(PolynomialRatio, digitalfilter(responsetype, designmethod))
  numerator_coefs = coefb(tf)
  denominator_coefs = coefa(tf)

Filter the data in ``x``, sampled at 1000 Hz, with a 4th order
Butterworth bandpass filter between 10 and 40 Hz::

  responsetype = Bandpass(10, 40; fs=1000)
  designmethod = Butterworth(4)
  filt(digitalfilter(responsetype, designmethod), x)

Filter the data in ``x``, sampled at 50 Hz, with a 64 tap Hanning
window FIR lowpass filter at 5 Hz::

  responsetype = Lowpass(5; fs=50)
  designmethod = FIRWindow(hanning(64))
  filt(digitalfilter(responsetype, designmethod), x)
