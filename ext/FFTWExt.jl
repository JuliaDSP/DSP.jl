module FFTWExt

using DSP: FFTBackend
using FFTW
import DSP.FFTBackend: AbstractFFTBackend, fft, ifft, ifft!, bfft, rfft, irfft, brfft,
                        plan_fft, plan_ifft, plan_bfft, plan_rfft, plan_irfft, plan_brfft,
                        plan_fft!, plan_bfft!

"""
    FFTWBackend <: AbstractFFTBackend

FFT backend using FFTW.jl. This provides optimized FFT implementations for Float32 and Float64.
Available when FFTW.jl is loaded.

# Example
```julia
using DSP
using FFTW
DSP.set_fft_backend!(DSP.FFTBackend.FFTWBackend())
```
"""
struct FFTWBackend <: AbstractFFTBackend end

# Transform implementations - delegate to FFTW
fft(::FFTWBackend, x) = FFTW.fft(x)
ifft(::FFTWBackend, x) = FFTW.ifft(x)
ifft!(::FFTWBackend, x) = FFTW.ifft!(x)
bfft(::FFTWBackend, x) = FFTW.bfft(x)
rfft(::FFTWBackend, x) = FFTW.rfft(x)
irfft(::FFTWBackend, x, d) = FFTW.irfft(x, d)
brfft(::FFTWBackend, x, d) = FFTW.brfft(x, d)

# Planning - delegate to FFTW (returns FFTW plan objects directly)
plan_fft(::FFTWBackend, x; flags=FFTW.ESTIMATE, kw...) = FFTW.plan_fft(x; flags=flags, kw...)
plan_ifft(::FFTWBackend, x; flags=FFTW.ESTIMATE, kw...) = FFTW.plan_ifft(x; flags=flags, kw...)
plan_bfft(::FFTWBackend, x; flags=FFTW.ESTIMATE, kw...) = FFTW.plan_bfft(x; flags=flags, kw...)
plan_rfft(::FFTWBackend, x; flags=FFTW.ESTIMATE, kw...) = FFTW.plan_rfft(x; flags=flags, kw...)
plan_irfft(::FFTWBackend, x, d; flags=FFTW.ESTIMATE, kw...) = FFTW.plan_irfft(x, d; flags=flags, kw...)
plan_brfft(::FFTWBackend, x, d; flags=FFTW.ESTIMATE, kw...) = FFTW.plan_brfft(x, d; flags=flags, kw...)
plan_fft!(::FFTWBackend, x; flags=FFTW.ESTIMATE, kw...) = FFTW.plan_fft!(x; flags=flags, kw...)
plan_bfft!(::FFTWBackend, x; flags=FFTW.ESTIMATE, kw...) = FFTW.plan_bfft!(x; flags=flags, kw...)

# Set FFTW as default when loaded (for better performance)
function __init__()
    FFTBackend.set_fft_backend!(FFTWBackend())
end

end # module
