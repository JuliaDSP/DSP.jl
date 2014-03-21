module DSP

include("util.jl")
include("windows.jl")
include("periodogram.jl")
include("fftfilt.jl")
include("filter_design.jl")

using DSP.Windows, DSP.Periodogram, DSP.FFTFilt, DSP.FilterDesign, DSP.Util

export
       # Util
       unwrap!, unwrap, hilbert,
       # Windows
       rect, hanning, hamming, tukey, cosine, lanczos,
       triang, bartlett, gaussian, bartlett_hann, blackman,
       kaiser, dpss,
       # Periodogram
       arraysplit, periodogram, welch_pgram, bartlett_pgram, spectrogram,
       # FFTFilt
       fftfilt, firfilt,
       # FilterDesign
       ZPKFilter, TFFilter, BiquadFilter, SOSFilter, Butterworth, Lowpass,
       Highpass, Bandpass, Bandstop, analogfilter, digitalfilter
end
