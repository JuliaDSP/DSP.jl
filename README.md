DSP.jl provides a number of common DSP routines in Julia.  So far, the following functions are
implemented:

### Periodogram estimation
* `arraysplit(s, n::Integer, m::Integer)` - Split an array into arrays of length `n` with overlapping regions of length `m`.
* `periodogram(s)` - Compute periodogram of a signal by FFT.
* `welch_pgram(s, n, m)` - Compute Welch periodogram of a signal based on `n` segments with overlap `m`.
* `bartlett_pgram(s, n)` - Compute Bartlett periodogram. This is equivalent to `welch_pgram(s, n, 0)`.

### Window functions
* `rect(n::Integer)` - Rectangular window function of length `n`.
* `hanning(n::Integer)` - Hanning window of length `n`.
* `hamming(n::Integer)` - Hamming window of length `n`.
* `tukey(n::Integer, alpha::Real)` - Tukey window of length `n`, parameterized by `alpha`. For `alpha` = 0, the window is equivalent to a rectangular window. For `alpha` = 1, the window is a Hann window.
* `cosine(n::Integer)` - Cosine window of length N.  Also called the sine window for obvious reasons.
* `lanczos(n::Integer)` - Lanczos window of length `n`.
* `triang(n::Integer)` - Triangular window of length `n`.
* `bartlett(n::Integer)` - Bartlett window of length `n`.
* `gaussian(n::Integer, sigma::Real)` - Gaussian window of length N parameterized by the standard deviation `sigma`
* `bartlett_hann(n::Integer)` - Bartlett-Hann window of length `n`.
* `blackman(n::Integer)` - "Exact" Blackman window, alpha=0.16.
* `kaiser(n::Integer, alpha::Real)` - Kaiser window parameterized by `alpha`.

### FFT-based (overlap-save) FIR filtering
* `fftfilt(b, x)` - Perform FFT-based filtering of `x` using filter `b`.
* `firfilt(b, x)` - Filter `x` using filter `b`, using `filt` or `fftfilt` depending on the lengths of `b` and `x`.

### Filter design functions
* `analogfilter(responsetype, designmethod)` - Construct an analog filter.
* `digitalfilter(responsetype, designmethod)` - Construct a digital filter.

### Filter response types
* `Lowpass(Wn)` - Low pass filter with normalized cutoff frequency `Wn`.
* `Highpass(Wn)` - High pass filter with normalized cutoff frequency `Wn`.
* `Bandpass(Wn1, Wn2)` - Band pass filter with normalized pass band `(Wn1, Wn2)`.
* `Bandstop(Wn1, Wn2)` - Band stop filter with normalized stop band `(Wn1, Wn2)`.

### Filter design methods
* `Butterworth(N::Integer)` - `N` pole Butterworth filter.
