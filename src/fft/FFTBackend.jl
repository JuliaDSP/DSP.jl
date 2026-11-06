module FFTBackend

using AbstractFFTs: fftfreq, rfftfreq, fftshift, ifftshift
export fftfreq, rfftfreq, fftshift, ifftshift

# FFT-compatible types (replaces FFTW.fftwReal etc.)
const FFTReal = Union{Float32, Float64}
const FFTComplex = Union{ComplexF32, ComplexF64}
const FFTNumber = Union{FFTReal, FFTComplex}
export FFTReal, FFTComplex, FFTNumber

# FFT flags (for compatibility with FFTW flags interface)
const ESTIMATE = UInt32(64)   # FFTW.ESTIMATE value
const MEASURE = UInt32(0)     # FFTW.MEASURE value
const NO_FLAG = UInt32(64)    # Default (ignored by non-FFTW backends)
export ESTIMATE, MEASURE, NO_FLAG

# Abstract backend type
abstract type AbstractFFTBackend end
export AbstractFFTBackend

# Interface functions - backends must implement these
function fft end
function ifft end
function ifft! end
function bfft end
function rfft end
function irfft end
function brfft end

function plan_fft end
function plan_ifft end
function plan_bfft end
function plan_rfft end
function plan_irfft end
function plan_brfft end
function plan_fft! end
function plan_bfft! end

export fft, ifft, ifft!, bfft, rfft, irfft, brfft
export plan_fft, plan_ifft, plan_bfft, plan_rfft, plan_irfft, plan_brfft, plan_fft!, plan_bfft!

# Global default backend (mutable)
const _default_backend = Ref{AbstractFFTBackend}()

"""
    set_fft_backend!(b::AbstractFFTBackend)

Set the default FFT backend for DSP.jl.

# Example
```julia
using DSP
DSP.set_fft_backend!(DSP.FFTABackend())
```
"""
function set_fft_backend!(b::AbstractFFTBackend)
    _default_backend[] = b
end
export set_fft_backend!

"""
    get_fft_backend()

Get the current default FFT backend.
"""
function get_fft_backend()
    if !isassigned(_default_backend)
        error("No FFT backend set. This should not happen - please report this as a bug.")
    end
    return _default_backend[]
end
export get_fft_backend

# Convenience methods that use default backend
fft(x) = fft(get_fft_backend(), x)
ifft(x) = ifft(get_fft_backend(), x)
ifft!(x) = ifft!(get_fft_backend(), x)
bfft(x) = bfft(get_fft_backend(), x)
rfft(x) = rfft(get_fft_backend(), x)
irfft(x, d) = irfft(get_fft_backend(), x, d)
brfft(x, d) = brfft(get_fft_backend(), x, d)

plan_fft(x; kw...) = plan_fft(get_fft_backend(), x; kw...)
plan_ifft(x; kw...) = plan_ifft(get_fft_backend(), x; kw...)
plan_bfft(x; kw...) = plan_bfft(get_fft_backend(), x; kw...)
plan_rfft(x; kw...) = plan_rfft(get_fft_backend(), x; kw...)
plan_irfft(x, d; kw...) = plan_irfft(get_fft_backend(), x, d; kw...)
plan_brfft(x, d; kw...) = plan_brfft(get_fft_backend(), x, d; kw...)
plan_fft!(x; kw...) = plan_fft!(get_fft_backend(), x; kw...)
plan_bfft!(x; kw...) = plan_bfft!(get_fft_backend(), x; kw...)

end # module
