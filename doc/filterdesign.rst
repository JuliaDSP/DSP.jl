:mod:`FilterDesign` - filter design functions
=============================================

.. function:: analogfilter(responsetype, designmethod)

    Construct an analog filter.

.. function:: digitalfilter(responsetype, designmethod)

    Construct a digital filter.

Filter response types
---------------------

.. function:: Lowpass(Wn)

    Low pass filter with normalized cutoff frequency Wn.

.. function:: Highpass(Wn)

    High pass filter with normalized cutoff frequency Wn.

.. function:: Bandpass(Wn1, Wn2)

    Band pass filter with normalized pass band (Wn1, Wn2).

.. function:: Bandstop(Wn1, Wn2)

    Band stop filter with normalized stop band (Wn1, Wn2).

Filter design methods
---------------------

.. function:: Butterworth(N::Integer) 

    N pole Butterworth filter.

Examples
--------

TODO
