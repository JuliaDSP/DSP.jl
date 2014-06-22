:mod:`FFTFilt` - FFT-based FIR filtering
===========================================

.. function:: fftfilt(b, x)

    Perform overlap-save filtering of ``x`` using filter ``b``.

.. function:: firfilt(b, x)

    Filter ``x`` using filter ``b``, using ``filt`` or ``fftfilt`` depending on the lengths of ``b`` and ``x``.
