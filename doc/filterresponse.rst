:mod:`FilterResponse` - Frequency response of filters
=====================================================

Digital filter response
-----------------------

.. function:: freqz(filter, w)

    Frequency response of a digital ``filter`` at normalised frequency ``w`` in radians/sample.

.. function:: freqz(filter, w)

    Frequency response of a digital ``filter`` at normalised frequencies ``w`` in radians/sample.

.. function:: freqz(filter, hz, fs)

    Frequency response of a digital ``filter`` at frequencies ``hz`` with sampling rate ``fs``.

Analog filter response
-----------------------

.. function:: freqs(filter, w)

    Frequency response of an analog ``filter`` at normalised frequency ``w`` in radians/sample.

.. function:: freqs(filter, w)

    Frequency response of an analog ``filter`` at normalised frequencies ``w`` in radians/sample.

.. function:: freqs(filter, hz, fs)

    Frequency response of an analog ``filter`` at frequencies ``hz`` with sampling rate ``fs``.
