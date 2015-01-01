module Filters
using Polynomials, ..Util

include("types.jl")
export Filter, ZPKFilter, TFFilter, BiquadFilter, SOSFilter, coefa, coefb

include("filt.jl")
export filtfilt, fftfilt, firfilt

include("design.jl")
export FilterType, Butterworth, Chebyshev1, Chebyshev2, Elliptic,
       Lowpass, Highpass, Bandpass, Bandstop, analogfilter,
       digitalfilter

include("response.jl")
export freqs, freqz, phasez, impz, stepz
end
