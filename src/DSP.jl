module DSP

include("windows.jl")
include("periodogram.jl")
include("fftfilt.jl")
include("filter_design.jl")

using DSP.Windows, DSP.Periodogram, DSP.FFTFilt, DSP.FilterDesign

export
       # Windows
       rect, hanning, hamming, tukey, cosine, lanczos, 
       triang, bartlett, gaussian, bartlett_hann, blackman, 
       kaiser,
       # Periodogram
       arraysplit, periodogram, welch_pgram, bartlett_pgram,
       # FFTFilt
       fftfilt, firfilt,
       # FilterDesign
       Butterworth, Lowpass, Highpass, Bandpass, Bandstop,
       analogfilter, digitalfilter
end
