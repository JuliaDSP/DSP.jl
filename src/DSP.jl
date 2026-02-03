module DSP

using LinearAlgebra: Transpose, mul!, rmul!
using IterTools: subsets

# FFT Backend system - must be included first
include("fft/FFTBackend.jl")
include("fft/FFTAImpl.jl")
using .FFTBackend
using .FFTAImpl

# Re-export FFT backend API
# Note: FFT functions (fft, rfft, etc.) are NOT exported to avoid conflicts with FFTW.jl
# Users should either:
# 1. Use FFTW.jl's exports: `using FFTW; fft(x)`
# 2. Qualify with DSP: `DSP.fft(x)` (uses the current backend)
# The backend selection (FFTABackend vs FFTWBackend) is automatic based on loaded packages.
export FFTABackend, set_fft_backend!, get_fft_backend
export fftfreq, rfftfreq, fftshift, ifftshift

export conv, conv!, deconv, filt, filt!, xcorr

# This function has methods added in `periodograms` but is not exported,
# so we define it here so one can do `DSP.allocate_output` instead of
# `DSP.Periodograms.allocate_output`.
function allocate_output end

include("dspbase.jl")

include("util.jl")
include("unwrap.jl")
include("windows.jl")
include("periodograms.jl")
include("Filters/Filters.jl")
include("lpc.jl")
include("estimation.jl")
include("diric.jl")

using Reexport
@reexport using .Util, .Windows, .Periodograms, .Filters, .LPC, .Unwrap, .Estimation

include("deprecated.jl")
end
