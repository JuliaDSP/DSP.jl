module  Filters
using   Polynomials,
        ..Util

include("types.jl")
export  Filter,
        ZPKFilter,
        TFFilter,
        BiquadFilter,
        SOSFilter,
        coefa,
        coefb

include("filt.jl")
export  filtfilt,
        fftfilt,
        firfilt

include("design.jl")
export  FilterType,
        Butterworth,
        Chebyshev1,
        Chebyshev2,
        Elliptic,
        Lowpass,
        Highpass,
        Bandpass,
        Bandstop,
        analogfilter,
        digitalfilter,
        kaiserlength,
        firdes

include("response.jl")    
export  freqs,
        freqz

include("stream_filt.jl")
export  FIRFilter,
        outputlength,
        inputlength,
        reset

include("naive_resampling.jl")
export  naivefilt

end