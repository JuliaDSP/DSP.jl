:mod:`FilterResponse` - Frequency response of filters
=====================================================

Digital filter response
-----------------------

.. function:: freqz(filter::Filter, w::Number)

    Frequency response of a digital ``filter`` at normalised frequency ``w`` in radians/sample.

.. function:: freqz(filter::Filter, w::AbstractVector)

    Frequency response of a digital ``filter`` at normalised frequencies ``w`` in radians/sample.

.. function:: freqz(filter::Filter, hz::Union(Number, AbstractVector), fs::Integer)

    Frequency response of a digital ``filter`` at frequencies ``hz`` with sampling rate ``fs``.

Analog filter response
-----------------------

.. function:: freqs(filter::Filter, w::Number)

    Frequency response of an analog ``filter`` at normalised frequency ``w`` in radians/sample.

.. function:: freqs(filter::Filter, w::AbstractVector)

    Frequency response of an analog ``filter`` at normalised frequencies ``w`` in radians/sample.

.. function:: freqs(filter::Filter, hz::Union(Number, AbstractVector), fs::Integer)

    Frequency response of an analog ``filter`` at frequencies ``hz`` with sampling rate ``fs``.
