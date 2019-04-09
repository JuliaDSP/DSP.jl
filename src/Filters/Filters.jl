module Filters
using ..Unwrap
using Polynomials, ..Util
import Base: *
using LinearAlgebra: I, mul!, rmul!
using Statistics: middle
import ..DSP: filt, filt!
using FFTW

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
