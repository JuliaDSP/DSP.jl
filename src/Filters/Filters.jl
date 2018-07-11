module Filters
using ..DSP: @importffts, mul!, rmul!
using Polynomials, ..Util
import Base: *
using Compat: copyto!, undef
import Compat
using Compat.LinearAlgebra: I
@importffts
if VERSION >= v"0.7.0-DEV.602"
    import ..DSP: filt, filt!
else
    import Base: filt, filt!
end

include("coefficients.jl")
export FilterCoefficients,
        ZeroPoleGain,
        PolynomialRatio,
        Biquad,
        SecondOrderSections,
        coefa,
        coefb

include("filt.jl")
export  DF2TFilter,
        filtfilt,
        fftfilt

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
        iirnotch,
        kaiserord,
        FIRWindow,
        resample_filter

include("response.jl")
export  freqs,
        freqz,
        phasez,
        impz,
        stepz

include("stream_filt.jl")
export  FIRFilter,
        outputlength,
        inputlength,
        reset!,
        resample,
        setphase!,
        timedelay

include("remez_fir.jl")
export  remez

end
