module Filters
using ..Unwrap
using ..Util
using Polynomials: LaurentPolynomial, Polynomial, coeffs, derivative, fromroots, roots, indeterminate

import Base: *
using LinearAlgebra: I, mul!, rmul!
using Statistics: middle
using SpecialFunctions: ellipk
import ..DSP: filt, filt!, optimalfftfiltlength, os_fft_complexity, SMALL_FILT_CUTOFF
import Compat
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
export DF2TFilter,
        filtfilt,
        tdfilt,
        tdfilt!,
        fftfilt,
        fftfilt!

include("design.jl")
export FilterType,
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

include("filt_order.jl")
export buttord,
        ellipord,
        cheb1ord,
        cheb2ord,
        remezord


include("response.jl")
export freqresp,
        phaseresp,
        grpdelay,
        impresp,
        stepresp

include("stream_filt.jl")
export FIRFilter,
        outputlength,
        inputlength,
        reset!,
        resample,
        setphase!,
        timedelay

include("remez_fir.jl")
export remez,
        RemezFilterType,
        filter_type_bandpass,
        filter_type_differentiator,
        filter_type_hilbert

end
