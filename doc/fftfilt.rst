:mod:`FFTFilt` - FFT-based FIR filtering
===========================================

.. function:: fftfilt{T<:Union(Float32,Float64)}(b::Vector{T}, x::Vector{T})

    Perform overlap-save filtering of ``x`` using filter ``b``.

.. function:: firfilt{T<:Union(Float32,Float64)}(b::Vector{T}, x::Vector{T})

    Filter ``x`` using filter ``b``, using ``filt`` or ``fftfilt`` depending on the lengths of ``b`` and ``x``.
