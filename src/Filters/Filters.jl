module Filters
using ..DSP: @importffts, mul!, rmul!
using ..Unwrap
using Polynomials, ..Util
import Base: *
using Compat: copyto!, undef, argmin
import Compat
using Compat.LinearAlgebra: I
using Compat.Statistics: middle
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
        tdfilt,
        tdfilt!,
        fftfilt,
        fftfilt!

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
export  remez,
        RemezFilterType,
        filter_type_bandpass,
        filter_type_differentiator,
        filter_type_hilbert

end
