var documenterSearchIndex = {"docs": [

{
    "location": "contents/#",
    "page": "Contents",
    "title": "Contents",
    "category": "page",
    "text": ""
},

{
    "location": "contents/#Welcome-to-DSP.jl\'s-documentation!-1",
    "page": "Contents",
    "title": "Welcome to DSP.jl\'s documentation!",
    "category": "section",
    "text": "Contents:DSP.jl provides a number of common DSP routines in Julia. So far, the following submodules are implemented:Pages = [\"periodograms.md\",\n    \"estimation.md\",\n    \"windows.md\",\n    \"filters.md\",\n    \"util.md\",\n    \"convolutions.md\",\n    \"lpc.md\",\n    \"index.md\",\n]"
},

{
    "location": "periodograms/#",
    "page": "Periodograms - periodogram estimation",
    "title": "Periodograms - periodogram estimation",
    "category": "page",
    "text": ""
},

{
    "location": "periodograms/#DSP.Periodograms.arraysplit",
    "page": "Periodograms - periodogram estimation",
    "title": "DSP.Periodograms.arraysplit",
    "category": "function",
    "text": "arraysplit(s, n, m)\n\nSplit an array into arrays of length n with overlapping regions of length m. Iterating or indexing the returned AbstractVector always yields the same Vector with different contents.\n\n\n\n\n\n"
},

{
    "location": "periodograms/#DSP.Periodograms.periodogram-Union{Tuple{AbstractArray{T,1}}, Tuple{T}} where T<:Number",
    "page": "Periodograms - periodogram estimation",
    "title": "DSP.Periodograms.periodogram",
    "category": "method",
    "text": "periodogram(s; onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, window=nothing)\n\nComputes periodogram of a signal by FFT and returns a Periodogram object.\n\nFor real signals, the two-sided periodogram is symmetric and this function returns a one-sided (real only) periodogram by default. A two-sided periodogram can be obtained by setting onesided=false.\n\nnfft specifies the number of points to use for the Fourier transform. If length(s) < nfft, then the input is padded with zeros. By default, nfft is the closest size for which the Fourier transform can be computed with maximal efficiency.\n\nfs is the sample rate of the original signal, and window is an optional window function or vector to be applied to the original signal before computing the Fourier transform. The computed periodogram is normalized so that the area under the periodogram is equal to the uncentered variance (or average power) of the original signal.\n\n\n\n\n\n"
},

{
    "location": "periodograms/#DSP.Periodograms.welch_pgram",
    "page": "Periodograms - periodogram estimation",
    "title": "DSP.Periodograms.welch_pgram",
    "category": "function",
    "text": "welch_pgram(s, n=div(length(s), 8), noverlap=div(n, 2); onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, window=nothing)\n\nComputes the Welch periodogram of a signal s based on segments with n samples with overlap of noverlap samples, and returns a Periodogram object. For a Bartlett periodogram, set noverlap=0. See periodogram for description of optional keyword arguments.\n\n\n\n\n\n"
},

{
    "location": "periodograms/#DSP.Periodograms.mt_pgram",
    "page": "Periodograms - periodogram estimation",
    "title": "DSP.Periodograms.mt_pgram",
    "category": "function",
    "text": "mt_pgram(s; onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, nw=4, ntapers=iceil(2nw)-1, window=dpss(length(s), nw, ntapers))\n\nComputes the multitaper periodogram of a signal s.\n\nIf window is not specified, the signal is tapered with ntapers discrete prolate spheroidal sequences with time-bandwidth product nw. Each sequence is equally weighted; adaptive multitaper is not (yet) supported.\n\nIf window is specified, each column is applied as a taper. The sum of periodograms is normalized by the total sum of squares of window.\n\nSee also: dpss\n\n\n\n\n\n"
},

{
    "location": "periodograms/#DSP.Periodograms.spectrogram",
    "page": "Periodograms - periodogram estimation",
    "title": "DSP.Periodograms.spectrogram",
    "category": "function",
    "text": "spectrogram(s, n=div(length(s), 8), noverlap=div(n, 2); onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, window=nothing)\n\nComputes the spectrogram of a signal s based on segments with n samples with overlap of noverlap samples, and returns a Spectrogram object. See periodogram for description of optional keyword arguments.\n\n\n\n\n\n"
},

{
    "location": "periodograms/#DSP.Periodograms.stft",
    "page": "Periodograms - periodogram estimation",
    "title": "DSP.Periodograms.stft",
    "category": "function",
    "text": "stft(s, n=div(length(s), 8), noverlap=div(n, 2); onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, window=nothing)\n\nComputes the STFT of a signal s based on segments with n samples with overlap of noverlap samples, and returns a matrix containing the STFT coefficients. See periodogram for description of optional keyword arguments.\n\n\n\n\n\n"
},

{
    "location": "periodograms/#DSP.Periodograms.periodogram-Union{Tuple{AbstractArray{T,2}}, Tuple{T}} where T<:Real",
    "page": "Periodograms - periodogram estimation",
    "title": "DSP.Periodograms.periodogram",
    "category": "method",
    "text": "periodogram(s::AbstractMatrix; nfft=nextfastfft(size(s)), fs=1, radialsum=false, radialavg=false)\n\nComputes periodogram of a 2-d signal by FFT and returns a Periodogram2 object.\n\nReturns a 2-d periodogram by default. A radially summed or averaged periodogram is returned as a Periodogram object if radialsum or  radialavg is true, respectively.\n\nnfft specifies the number of points to use for the Fourier transform. If size(s) < nfft, then the input is padded with zeros. By default, nfft is the closest size for which the Fourier transform can be computed with maximal efficiency. fs is the sample rate of the original signal in both directions.\n\nFor radialsum=true the value of power[k] is proportional to frac1Nsum_kleq kk+1 Xk^2. For radialavg=true it is proportional to frac1N kleq kk+1 sum_kleq kk+1 Xk^2. The computation of |k\'| takes into account non-square signals by scaling the coordinates of the wavevector accordingly.\n\n\n\n\n\n"
},

{
    "location": "periodograms/#DSP.Periodograms.freq",
    "page": "Periodograms - periodogram estimation",
    "title": "DSP.Periodograms.freq",
    "category": "function",
    "text": "freq(p)\n\nReturns the frequency bin centers for a given Periodogram or Spectrogram object.\n\nReturns a tuple of frequency bin centers for a given Periodogram2 object.\n\nSee also: fftfreq, rfftfreq\n\n\n\n\n\n"
},

{
    "location": "periodograms/#DSP.Periodograms.power",
    "page": "Periodograms - periodogram estimation",
    "title": "DSP.Periodograms.power",
    "category": "function",
    "text": "power(p)\n\nFor a Periodogram, returns the computed power at each frequency as a Vector.\n\nFor a Spectrogram, returns the computed power at each frequency and time bin as a Matrix. Dimensions are frequency × time.\n\n\n\n\n\n"
},

{
    "location": "periodograms/#Base.Libc.time",
    "page": "Periodograms - periodogram estimation",
    "title": "Base.Libc.time",
    "category": "function",
    "text": "time(p)\n\nReturns the time bin centers for a given Spectrogram object.\n\n\n\n\n\n"
},

{
    "location": "periodograms/#Periodograms-periodogram-estimation-1",
    "page": "Periodograms - periodogram estimation",
    "title": "Periodograms - periodogram estimation",
    "category": "section",
    "text": "arraysplit\nperiodogram(s::AbstractVector{T}) where T <: Number\nwelch_pgram\nmt_pgram\nspectrogram\nstft\nperiodogram(s::AbstractMatrix{T}) where T <: Real\nfreq\npower\ntime"
},

{
    "location": "estimation/#",
    "page": "Estimation - parametric estimation functions",
    "title": "Estimation - parametric estimation functions",
    "category": "page",
    "text": ""
},

{
    "location": "estimation/#DSP.Estimation.esprit",
    "page": "Estimation - parametric estimation functions",
    "title": "DSP.Estimation.esprit",
    "category": "function",
    "text": "esprit(x::AbstractArray, M::Integer, p::Integer, Fs::Real=1.0)\n\nESPRIT [Roy1986] algorithm for frequency estimation. Estimation of Signal Parameters via Rotational Invariance Techniques\n\nGiven length N signal \"x\" that is the sum of p sinusoids of unknown frequencies, estimate and return an array of the p frequencies.\n\nArguments\n\nx::AbstractArray: complex length N signal array\nM::Integer: size of correlation matrix, must be <= N.     The signal subspace is computed from the SVD of an M x (N-M+1) signal matrix     formed from N-M+1 length-M shifts of the signal x in its columns.     For best performance for 1 sinusoid, use M = (N+1)/3 (according to van der Veen and Leus).     For faster execution (due to smaller SVD), use small M or small N-M\np::Integer: number of sinusoids to estimate.\nFs::Float64: sampling frequency, in Hz.\n\nReturns\n\nlength p real array of frequencies in units of Hz.\n\n[Roy1986]: R Roy, A Paulraj and T Kailath, ESPRIT - A subspace approach to estimation of parameters of cisoids in noise, IEEE Trans. Acoustics, Speech, Signal Process., 34, 1340-1342 (1986). url.\n\n\n\n\n\n"
},

{
    "location": "estimation/#Estimation-parametric-estimation-functions-1",
    "page": "Estimation - parametric estimation functions",
    "title": "Estimation - parametric estimation functions",
    "category": "section",
    "text": "esprit"
},

{
    "location": "windows/#",
    "page": "Windows - window functions",
    "title": "Windows - window functions",
    "category": "page",
    "text": ""
},

{
    "location": "windows/#DSP.Windows.rect",
    "page": "Windows - window functions",
    "title": "DSP.Windows.rect",
    "category": "function",
    "text": "  ┌──────────────────────────────────────────────────────────────────────┐ \n1 │ ▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▘│ \n  │                                                                      │ \n  │                                                                      │ \n  │                                                                      │ \n  │                                                                      │ \n  │                                                                      │ \n  │                                                                      │ \n  │                                                                      │ \n  │                                                                      │ \n  │                                                                      │ \n  │                                                                      │ \n  │                                                                      │ \n  │                                                                      │ \n  │                                                                      │ \n0 │                                                                      │ \n  └──────────────────────────────────────────────────────────────────────┘ \n  0                                                                     70\n\n\nrect(n::Integer; padding::Integer=0, zerophase::Bool=false)\nrect(dims; padding=0, zerophase=false)\n\nRectangular window of length n, padded with padding zeros. This window is 1 within the window, and 0 outside of it.\n\nProviding a dims Tuple rather than a single n constructs a 2D window. padding and zerophase can then be given either as a single value for both horizontal and vertical or a 2-tuple to specify them separately.\n\nIf zerophase is false (the default) the window is centered around index (n+1)/2, which is commonly used for FIR filter design. These are often used for FIR filter design, and usually odd-length. Note that for even-length windows this will cause the window center to be located between samples.\n\nIf zerophase is true the window is centered around index 1 (with the negative half wrapped to the end of the vector). Additionally this creates a \"periodic\" window, which means that if there is no padding then the left and right endpoints of the window wrap around to the same sample, so the window length is the same as an n+1-length non-zerophase window. Alternatively you can think of the continuous zerophase window being of width n and the non-zerophase window being of length n-1. zerophase windows are often used in FFT processing, and are usually even-length.\n\n\n\n\n\n"
},

{
    "location": "windows/#DSP.Windows.hanning",
    "page": "Windows - window functions",
    "title": "DSP.Windows.hanning",
    "category": "function",
    "text": "  ┌──────────────────────────────────────────────────────────────────────┐ \n1 │                             ▗▄▞▀▀▀▀▀▀▀▄▄                             │ \n  │                           ▄▞▘           ▀▄▖                          │ \n  │                         ▄▀                ▝▚▖                        │ \n  │                       ▗▞                    ▝▄                       │ \n  │                      ▞▘                      ▝▚▖                     │ \n  │                    ▗▀                          ▝▚                    │ \n  │                   ▞▘                             ▀▖                  │ \n  │                 ▗▞                                ▝▄                 │ \n  │                ▄▘                                   ▚▖               │ \n  │              ▗▞                                      ▝▄              │ \n  │             ▄▘                                         ▚▖            │ \n  │           ▗▀                                            ▝▚           │ \n  │         ▗▞▘                                               ▀▄         │ \n  │       ▄▀▘                                                   ▀▚▖      │ \n0 │ ▄▄▄▄▀▀                                                        ▝▀▚▄▄▄▖│ \n  └──────────────────────────────────────────────────────────────────────┘ \n  0                                                                     70\n\n\nhanning(n::Integer; padding::Integer=0, zerophase::Bool=false)\nhanning(dims; padding=0, zerophase=false)\n\nHanning window of length n with padding zeros. The Hanning (or Hann) window is a raised-cosine window that reaches zero at the endpoints.\n\nThe window is defined by sampling the continuous function:\n\n       1 + cos(2πx)\nw(x) = ──────────── = cos²(πx)\n            2\n\nin the range [-0.5, 0.5]\n\nThe hanning window satisfies the Constant Overlap-Add (COLA) property with an hop of 0.5, which means that adding together a sequence of delayed windows with 50% overlap will result in a constant signal. This is useful when synthesizing a signal from a number of overlapping frames (each with a roughly rectangular window), to eliminate windowing amplitude modulation.\n\nNote that the hanning window is the cosine window squared.\n\nProviding a dims Tuple rather than a single n constructs a 2D window. padding and zerophase can then be given either as a single value for both horizontal and vertical or a 2-tuple to specify them separately.\n\nIf zerophase is false (the default) the window is centered around index (n+1)/2, which is commonly used for FIR filter design. These are often used for FIR filter design, and usually odd-length. Note that for even-length windows this will cause the window center to be located between samples.\n\nIf zerophase is true the window is centered around index 1 (with the negative half wrapped to the end of the vector). Additionally this creates a \"periodic\" window, which means that if there is no padding then the left and right endpoints of the window wrap around to the same sample, so the window length is the same as an n+1-length non-zerophase window. Alternatively you can think of the continuous zerophase window being of width n and the non-zerophase window being of length n-1. zerophase windows are often used in FFT processing, and are usually even-length.\n\n\n\n\n\n"
},

{
    "location": "windows/#DSP.Windows.hamming",
    "page": "Windows - window functions",
    "title": "DSP.Windows.hamming",
    "category": "function",
    "text": "  ┌──────────────────────────────────────────────────────────────────────┐ \n1 │                             ▗▄▀▀▀▀▀▀▀▀▚▄                             │ \n  │                           ▄▀▘           ▀▚▖                          │ \n  │                         ▞▀                ▝▀▖                        │ \n  │                       ▄▀                    ▝▚                       │ \n  │                     ▗▀                        ▀▚                     │ \n  │                   ▗▞▘                           ▀▄                   │ \n  │                  ▄▘                               ▚▖                 │ \n  │                ▗▞                                  ▝▄                │ \n  │               ▞▘                                     ▀▖              │ \n  │             ▄▀                                        ▝▚▖            │ \n  │           ▗▞                                            ▝▄           │ \n  │         ▄▞▘                                               ▀▄▖        │ \n  │      ▗▄▀                                                    ▝▚▄      │ \n  │ ▄▄▄▞▀▘                                                         ▀▀▄▄▄▖│ \n0 │                                                                      │ \n  └──────────────────────────────────────────────────────────────────────┘ \n  0                                                                     70\n\n\nhamming(n::Integer; padding::Integer=0, zerophase::Bool=false)\nhamming(dims; padding=0, zerophase=false)\n\nHamming window of length n with padding zeros. The Hamming window does not reach zero at the endpoints and so has a shallower frequency roll-off when compared to the Hanning window, but is designed to cancel the first side-lobe.\n\nThe window is defined by sampling the continuous function:\n\nw(x) = 0.54 + 0.46*cos(2pi*x)\n\nin the range [-0.5, 0.5]\n\nProviding a dims Tuple rather than a single n constructs a 2D window. padding and zerophase can then be given either as a single value for both horizontal and vertical or a 2-tuple to specify them separately.\n\nIf zerophase is false (the default) the window is centered around index (n+1)/2, which is commonly used for FIR filter design. These are often used for FIR filter design, and usually odd-length. Note that for even-length windows this will cause the window center to be located between samples.\n\nIf zerophase is true the window is centered around index 1 (with the negative half wrapped to the end of the vector). Additionally this creates a \"periodic\" window, which means that if there is no padding then the left and right endpoints of the window wrap around to the same sample, so the window length is the same as an n+1-length non-zerophase window. Alternatively you can think of the continuous zerophase window being of width n and the non-zerophase window being of length n-1. zerophase windows are often used in FFT processing, and are usually even-length.\n\n\n\n\n\n"
},

{
    "location": "windows/#DSP.Windows.tukey",
    "page": "Windows - window functions",
    "title": "DSP.Windows.tukey",
    "category": "function",
    "text": "  ┌──────────────────────────────────────────────────────────────────────┐ \n1 │            ▗▞▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▄            │ \n  │           ▞▘                                             ▌           │ \n  │          ▗▘                                              ▝▌          │ \n  │         ▗▌                                                ▝▖         │ \n  │         ▞                                                  ▚         │ \n  │        ▗▘                                                  ▝▌        │ \n  │        ▞                                                    ▐        │ \n  │       ▐▘                                                     ▌       │ \n  │       ▌                                                      ▐▖      │ \n  │      ▐                                                        ▚      │ \n  │      ▌                                                        ▝▖     │ \n  │     ▞                                                          ▚     │ \n  │    ▞▘                                                           ▌    │ \n  │   ▗▘                                                            ▝▚   │ \n0 │ ▄▄▘                                                               ▚▄▖│ \n  └──────────────────────────────────────────────────────────────────────┘ \n  0                                                                     70\n\n\ntukey(n::Integer, α::Real; padding::Integer=0, zerophase::Bool=false)\ntukey(dims, α; padding=0, zerophase=false)\n\nTukey window of length n with padding zeros. The Tukey window has a flat top and reaches zero at the endpoints, with a sinusoidal transition area parameterized by α. For α == 0, the window is equivalent to a rectangular window. For α == 1, the window is a Hann window.\n\nThe window is defined by sampling the continuous function:\n\n       ⎛              ⎛    ⎛    1 - α⎞⎞\n       ⎜      1 + cos ⎜2πα ⎜x + ─────⎟⎟             1 - α\n       ⎜              ⎝    ⎝      2  ⎠⎠         x ≤ ─────\n       ⎜      ─────────────────────────               2\n       ⎜                  2\n       ⎜\nw(x) = ⎜      1                                 -α/2 < x ≤ α/2\n       ⎜\n       ⎜              ⎛    ⎛    1 - α⎞⎞\n       ⎜      1 + cos ⎜2πα ⎜x - ─────⎟⎟             1 - α\n       ⎜              ⎝    ⎝      2  ⎠⎠         x > ─────\n       ⎜      ─────────────────────────               2\n       ⎝                  2\n\nin the range [-0.5, 0.5]\n\nProviding a dims Tuple rather than a single n constructs a 2D window. α, padding and zerophase can then be given either as a single value for both horizontal and vertical or a 2-tuple to specify them separately.\n\nIf zerophase is false (the default) the window is centered around index (n+1)/2, which is commonly used for FIR filter design. These are often used for FIR filter design, and usually odd-length. Note that for even-length windows this will cause the window center to be located between samples.\n\nIf zerophase is true the window is centered around index 1 (with the negative half wrapped to the end of the vector). Additionally this creates a \"periodic\" window, which means that if there is no padding then the left and right endpoints of the window wrap around to the same sample, so the window length is the same as an n+1-length non-zerophase window. Alternatively you can think of the continuous zerophase window being of width n and the non-zerophase window being of length n-1. zerophase windows are often used in FFT processing, and are usually even-length.\n\n\n\n\n\n"
},

{
    "location": "windows/#DSP.Windows.cosine",
    "page": "Windows - window functions",
    "title": "DSP.Windows.cosine",
    "category": "function",
    "text": "  ┌──────────────────────────────────────────────────────────────────────┐ \n1 │                           ▗▄▞▀▀▀▀▀▀▀▀▀▀▀▄▄                           │ \n  │                        ▄▞▀▘               ▀▀▄▖                       │ \n  │                     ▗▞▀                      ▝▚▄                     │ \n  │                   ▄▀▘                           ▀▚▖                  │ \n  │                 ▄▀                                ▝▚▖                │ \n  │               ▄▀                                    ▝▚▖              │ \n  │             ▗▞                                        ▝▄             │ \n  │            ▞▘                                           ▀▖           │ \n  │          ▄▀                                              ▝▚▖         │ \n  │        ▗▞                                                  ▝▄        │ \n  │       ▄▘                                                     ▚▖      │ \n  │     ▗▞                                                        ▝▄     │ \n  │    ▄▘                                                           ▚▖   │ \n  │  ▗▞                                                              ▝▄  │ \n0 │ ▄▘                                                                 ▚▖│ \n  └──────────────────────────────────────────────────────────────────────┘ \n  0                                                                     70\n\n\ncosine(n::Integer; padding::Integer=0, zerophase::Bool=false)\ncosine(dims; padding=0, zerophase=false)\n\nCosine window of length n with padding zeros. The cosine window is the first lobe of a cosine function (with the zero crossings at +/- π as endpoints). Also called the sine window.\n\nThe window is defined by sampling the continuous function:\n\nw(x) = cos(πx)\n\nin the range [-0.5, 0.5]\n\nNote that the cosine window is the square root of the hanning window, so it is sometimes used when you are applying the window twice, such as the analysis and synthesis steps of an STFT.\n\nProviding a dims Tuple rather than a single n constructs a 2D window. padding and zerophase can then be given either as a single value for both horizontal and vertical or a 2-tuple to specify them separately.\n\nIf zerophase is false (the default) the window is centered around index (n+1)/2, which is commonly used for FIR filter design. These are often used for FIR filter design, and usually odd-length. Note that for even-length windows this will cause the window center to be located between samples.\n\nIf zerophase is true the window is centered around index 1 (with the negative half wrapped to the end of the vector). Additionally this creates a \"periodic\" window, which means that if there is no padding then the left and right endpoints of the window wrap around to the same sample, so the window length is the same as an n+1-length non-zerophase window. Alternatively you can think of the continuous zerophase window being of width n and the non-zerophase window being of length n-1. zerophase windows are often used in FFT processing, and are usually even-length.\n\n\n\n\n\n"
},

{
    "location": "windows/#DSP.Windows.lanczos",
    "page": "Windows - window functions",
    "title": "DSP.Windows.lanczos",
    "category": "function",
    "text": "  ┌──────────────────────────────────────────────────────────────────────┐ \n1 │                            ▗▄▞▀▀▀▀▀▀▀▀▀▄▄                            │ \n  │                         ▗▞▀▘             ▀▀▄                         │ \n  │                       ▄▀▘                   ▀▚                       │ \n  │                     ▄▀                        ▀▚▖                    │ \n  │                   ▄▀                            ▝▚▖                  │ \n  │                 ▗▀                                ▝▚                 │ \n  │               ▗▞▘                                   ▀▄               │ \n  │              ▞▘                                       ▀▖             │ \n  │            ▄▀                                          ▝▚▖           │ \n  │          ▗▞                                              ▝▄          │ \n  │         ▞▘                                                 ▀▖        │ \n  │       ▄▀                                                    ▝▚▖      │ \n  │     ▗▞                                                        ▝▄     │ \n  │   ▗▞▘                                                           ▀▄   │ \n0 │ ▄▞▘                                                               ▀▄▖│ \n  └──────────────────────────────────────────────────────────────────────┘ \n  0                                                                     70\n\n\nlanczos(n::Integer; padding::Integer=0, zerophase::Bool=false)\nlanczos(dims; padding=0, zerophase=false)\n\nLanczos window of length n with padding zeros. The Lanczos window is the main lobe of a sinc function.\n\nThe window is defined by sampling the continuous function:\n\n                  sin(2πx)\nw(x) = sinc(2x) = ────────\n                     2πx\n\nin the range [-0.5, 0.5]\n\nProviding a dims Tuple rather than a single n constructs a 2D window. padding and zerophase can then be given either as a single value for both horizontal and vertical or a 2-tuple to specify them separately.\n\nIf zerophase is false (the default) the window is centered around index (n+1)/2, which is commonly used for FIR filter design. These are often used for FIR filter design, and usually odd-length. Note that for even-length windows this will cause the window center to be located between samples.\n\nIf zerophase is true the window is centered around index 1 (with the negative half wrapped to the end of the vector). Additionally this creates a \"periodic\" window, which means that if there is no padding then the left and right endpoints of the window wrap around to the same sample, so the window length is the same as an n+1-length non-zerophase window. Alternatively you can think of the continuous zerophase window being of width n and the non-zerophase window being of length n-1. zerophase windows are often used in FFT processing, and are usually even-length.\n\n\n\n\n\n"
},

{
    "location": "windows/#DSP.Windows.triang",
    "page": "Windows - window functions",
    "title": "DSP.Windows.triang",
    "category": "function",
    "text": "  ┌──────────────────────────────────────────────────────────────────────┐ \n1 │                                ▄▞▀▀▀▄▖                               │ \n  │                            ▗▄▞▀      ▝▀▄▖                            │ \n  │                         ▗▄▀▘            ▝▀▄▖                         │ \n  │                      ▗▄▀▘                  ▝▀▄▖                      │ \n  │                   ▗▄▀▘                        ▝▀▄▖                   │ \n  │                ▗▄▀▘                              ▝▀▄▄                │ \n  │             ▗▄▀▘                                     ▀▚▄             │ \n  │          ▗▄▀▘                                           ▀▙▄          │ \n  │       ▄▄▀▘                                                 ▀▚▄       │ \n  │    ▄▞▀                                                        ▀▚▄    │ \n  │ ▄▞▀                                                              ▀▚▄ │ \n  │▀                                                                    ▀│ \n  │                                                                      │ \n  │                                                                      │ \n0 │                                                                      │ \n  └──────────────────────────────────────────────────────────────────────┘ \n  1                                                                      7\n\n\ntriang(n::Integer; padding::Integer=0, zerophase::Bool=false)\ntriang(dims; padding=0, zerophase=false)\n\nTriangular window of length n with padding zeros. The Triangular window does not reach zero at the endpoints. For odd n the triang window is the center n points of an n+2-point bartlett window (i.e. the samples just outside the window would be zero). For even n the window slope is the same as the n-1 window but delayed by a half sample so the zero points would be 1/2 sample past the ends of the window.\n\nThe window is defined by sampling the continuous function:\n\n        ⎛    2(n-1)\n        ⎜1 - ────── abs(x)     n is even\n        ⎜       n\nw(x) =  ⎜\n        ⎜    2(n-1)\n        ⎜1 - ────── abs(x)     n is odd\n        ⎝     n+1\n\nin the range [-0.5, 0.5].\n\nProviding a dims Tuple rather than a single n constructs a 2D window. padding and zerophase can then be given either as a single value for both horizontal and vertical or a 2-tuple to specify them separately.\n\nIf zerophase is false (the default) the window is centered around index (n+1)/2, which is commonly used for FIR filter design. These are often used for FIR filter design, and usually odd-length. Note that for even-length windows this will cause the window center to be located between samples.\n\nIf zerophase is true the window is centered around index 1 (with the negative half wrapped to the end of the vector). Additionally this creates a \"periodic\" window, which means that if there is no padding then the left and right endpoints of the window wrap around to the same sample, so the window length is the same as an n+1-length non-zerophase window. Alternatively you can think of the continuous zerophase window being of width n and the non-zerophase window being of length n-1. zerophase windows are often used in FFT processing, and are usually even-length.\n\nWhen zerophase is true substitute n+1 for n in the above window expressions.\n\n\n\n\n\n"
},

{
    "location": "windows/#DSP.Windows.bartlett",
    "page": "Windows - window functions",
    "title": "DSP.Windows.bartlett",
    "category": "function",
    "text": "  ┌──────────────────────────────────────────────────────────────────────┐ \n1 │                                 ▄▀▀▚▖                                │ \n  │                              ▗▄▀    ▝▚▄                              │ \n  │                            ▗▞▘         ▀▄                            │ \n  │                          ▄▞▘             ▀▄▖                         │ \n  │                        ▄▀                  ▝▚▖                       │ \n  │                     ▗▄▀                      ▝▄▄                     │ \n  │                   ▗▞▘                           ▀▄                   │ \n  │                 ▄▞▘                               ▀▄▖                │ \n  │               ▄▀                                    ▝▚▖              │ \n  │            ▗▞▀                                        ▝▀▄            │ \n  │          ▗▞▘                                             ▀▄          │ \n  │        ▄▀▘                                                 ▀▚▖       │ \n  │      ▄▀                                                      ▝▚▖     │ \n  │   ▗▞▀                                                          ▝▀▄   │ \n0 │ ▄▞▘                                                               ▀▄▖│ \n  └──────────────────────────────────────────────────────────────────────┘ \n  0                                                                     70\n\n\nbartlett(n::Integer; padding::Integer=0, zerophase::Bool=false)\nbartlett(dims; padding=0, zerophase=false)\n\nBartlett window of length n. The Bartlett window is a triangular window that reaches 0 at the endpoints. This is equivalent to convolving two rectangular windows of length (n-1)/2 and adding the zero endpoints. See triang for a window that does not reach zero at the endpoints.\n\nThe window is defined by sampling the continuous function:\n\n1 - abs(2x)\n\nin the range [-0.5, 0.5]\n\nProviding a dims Tuple rather than a single n constructs a 2D window. padding and zerophase can then be given either as a single value for both horizontal and vertical or a 2-tuple to specify them separately.\n\nIf zerophase is false (the default) the window is centered around index (n+1)/2, which is commonly used for FIR filter design. These are often used for FIR filter design, and usually odd-length. Note that for even-length windows this will cause the window center to be located between samples.\n\nIf zerophase is true the window is centered around index 1 (with the negative half wrapped to the end of the vector). Additionally this creates a \"periodic\" window, which means that if there is no padding then the left and right endpoints of the window wrap around to the same sample, so the window length is the same as an n+1-length non-zerophase window. Alternatively you can think of the continuous zerophase window being of width n and the non-zerophase window being of length n-1. zerophase windows are often used in FFT processing, and are usually even-length.\n\n\n\n\n\n"
},

{
    "location": "windows/#DSP.Windows.gaussian",
    "page": "Windows - window functions",
    "title": "DSP.Windows.gaussian",
    "category": "function",
    "text": "  ┌──────────────────────────────────────────────────────────────────────┐ \n1 │                              ▄▞▀▀▀▀▀▀▀▄▖                             │ \n  │                            ▄▀          ▝▚▖                           │ \n  │                          ▄▀              ▝▚▖                         │ \n  │                        ▗▞                  ▝▄                        │ \n  │                       ▄▘                     ▚                       │ \n  │                     ▗▞                        ▚▄                     │ \n  │                    ▄▘                           ▚▖                   │ \n  │                  ▗▀                              ▝▚                  │ \n  │                 ▞▘                                 ▀▖                │ \n  │               ▄▀                                    ▝▚▖              │ \n  │             ▄▀                                        ▝▚▖            │ \n  │           ▄▀                                            ▝▚▖          │ \n  │        ▄▞▀                                                ▝▀▄▖       │ \n  │   ▗▄▄▀▀                                                      ▝▀▚▄▄   │ \n0 │ ▀▀▘                                                               ▀▀▘│ \n  └──────────────────────────────────────────────────────────────────────┘ \n  0                                                                     70\n\n\ngaussian(n::Integer, σ::Real; padding::Integer=0, zerophase::Bool=false)\ngaussian(dims, σ; padding=0, zerophase=false)\n\nGives an n-sample gaussian window defined by sampling the function:\n\n        ⎛        2⎞\n        ⎜-1   ⎛x⎞ ⎟\n        ⎜── ⋅ ⎜─⎟ ⎟\n        ⎝ 2   ⎝σ⎠ ⎠\nw(x) = e\n\nin the range [-0.5,0.5]. This means that for σ=0.5 the endpoints of the window will correspond to 1 standard deviation away from the center.\n\nProviding a dims Tuple rather than a single n constructs a 2D window. σ, padding and zerophase can then be given either as a single value for both horizontal and vertical or a 2-tuple to specify them separately.\n\nIf zerophase is false (the default) the window is centered around index (n+1)/2, which is commonly used for FIR filter design. These are often used for FIR filter design, and usually odd-length. Note that for even-length windows this will cause the window center to be located between samples.\n\nIf zerophase is true the window is centered around index 1 (with the negative half wrapped to the end of the vector). Additionally this creates a \"periodic\" window, which means that if there is no padding then the left and right endpoints of the window wrap around to the same sample, so the window length is the same as an n+1-length non-zerophase window. Alternatively you can think of the continuous zerophase window being of width n and the non-zerophase window being of length n-1. zerophase windows are often used in FFT processing, and are usually even-length.\n\n\n\n\n\n"
},

{
    "location": "windows/#DSP.Windows.bartlett_hann",
    "page": "Windows - window functions",
    "title": "DSP.Windows.bartlett_hann",
    "category": "function",
    "text": "  ┌──────────────────────────────────────────────────────────────────────┐ \n1 │                              ▗▄▞▀▀▀▀▀▄▄                              │ \n  │                            ▄▀▘         ▀▚▖                           │ \n  │                          ▄▀              ▝▚▖                         │ \n  │                        ▄▀                  ▝▚▖                       │ \n  │                      ▗▀                      ▝▖                      │ \n  │                    ▗▞▘                        ▝▀▄                    │ \n  │                   ▄▘                             ▚▖                  │ \n  │                 ▗▞                                ▝▄                 │ \n  │                ▞▘                                   ▀▖               │ \n  │              ▄▀                                      ▝▚▖             │ \n  │            ▗▞                                          ▝▄            │ \n  │          ▗▞▘                                             ▀▄          │ \n  │        ▗▞▘                                                 ▀▄        │ \n  │      ▄▞▘                                                     ▀▄▖     │ \n0 │ ▄▄▄▀▀                                                          ▝▀▚▄▄▖│ \n  └──────────────────────────────────────────────────────────────────────┘ \n  0                                                                     70\n\n\nbartlett_hann(n::Integer; padding::Integer=0, zerophase::Bool=false)\nbartlett_hann(dims; padding=0, zerophase=false)\n\nBartlett-Hann window of length n with padding zeros. The Bartlett-Hann window is a weighted sum of the Bartlett and Hann windows.\n\nThe window is defined by sampling the continuous function:\n\nw(x) = 0.62 - 0.48*abs(x) + 0.38*cos(2π*x)\n\nin the range [-0.5, 0.5]\n\nProviding a dims Tuple rather than a single n constructs a 2D window. padding and zerophase can then be given either as a single value for both horizontal and vertical or a 2-tuple to specify them separately.\n\nIf zerophase is false (the default) the window is centered around index (n+1)/2, which is commonly used for FIR filter design. These are often used for FIR filter design, and usually odd-length. Note that for even-length windows this will cause the window center to be located between samples.\n\nIf zerophase is true the window is centered around index 1 (with the negative half wrapped to the end of the vector). Additionally this creates a \"periodic\" window, which means that if there is no padding then the left and right endpoints of the window wrap around to the same sample, so the window length is the same as an n+1-length non-zerophase window. Alternatively you can think of the continuous zerophase window being of width n and the non-zerophase window being of length n-1. zerophase windows are often used in FFT processing, and are usually even-length.\n\n\n\n\n\n"
},

{
    "location": "windows/#DSP.Windows.blackman",
    "page": "Windows - window functions",
    "title": "DSP.Windows.blackman",
    "category": "function",
    "text": "  ┌──────────────────────────────────────────────────────────────────────┐ \n1 │                               ▄▀▀▀▀▀▀▚▖                              │ \n  │                             ▄▀        ▝▚▖                            │ \n  │                           ▗▀            ▝▚                           │ \n  │                          ▞▘               ▀▖                         │ \n  │                        ▗▞                  ▝▄                        │ \n  │                       ▗▘                     ▚                       │ \n  │                      ▞▘                      ▝▚▖                     │ \n  │                     ▞                          ▝▖                    │ \n  │                   ▗▀                            ▝▚                   │ \n  │                  ▄▘                               ▚▖                 │ \n  │                ▗▞                                  ▝▄                │ \n  │               ▞▘                                     ▀▖              │ \n  │            ▗▄▀                                        ▝▚▄            │ \n  │         ▗▄▞▘                                             ▀▄▄         │ \n0 │ ▄▄▄▄▄▄▞▀▘                                                   ▀▀▄▄▄▄▄▄▖│ \n  └──────────────────────────────────────────────────────────────────────┘ \n  0                                                                     70\n\n\nblackman(n::Integer; padding::Integer=0, zerophase::Bool=false)\nblackman(dims; padding=0, zerophase=false)\n\nApproximates the \"Exact\" Blackman window. This is the generalized Blackman window with α = 0.16.\n\nThe window is defined by sampling the continuous function:\n\nw(x) = 0.42 + 0.5*cos(2π*x) + 0.08*cos(4π*x)\n\nin the range [-0.5, 0.5]\n\nProviding a dims Tuple rather than a single n constructs a 2D window. padding and zerophase can then be given either as a single value for both horizontal and vertical or a 2-tuple to specify them separately.\n\nIf zerophase is false (the default) the window is centered around index (n+1)/2, which is commonly used for FIR filter design. These are often used for FIR filter design, and usually odd-length. Note that for even-length windows this will cause the window center to be located between samples.\n\nIf zerophase is true the window is centered around index 1 (with the negative half wrapped to the end of the vector). Additionally this creates a \"periodic\" window, which means that if there is no padding then the left and right endpoints of the window wrap around to the same sample, so the window length is the same as an n+1-length non-zerophase window. Alternatively you can think of the continuous zerophase window being of width n and the non-zerophase window being of length n-1. zerophase windows are often used in FFT processing, and are usually even-length.\n\n\n\n\n\n"
},

{
    "location": "windows/#DSP.Windows.kaiser",
    "page": "Windows - window functions",
    "title": "DSP.Windows.kaiser",
    "category": "function",
    "text": "  ┌──────────────────────────────────────────────────────────────────────┐ \n1 │                               ▄▞▀▀▀▀▀▄▖                              │ \n  │                             ▄▀        ▝▚▖                            │ \n  │                           ▗▞            ▝▄                           │ \n  │                          ▗▘               ▚                          │ \n  │                         ▞▘                 ▀▖                        │ \n  │                        ▞                    ▝▖                       │ \n  │                      ▗▞                      ▐▖                      │ \n  │                     ▗▘                        ▝▚                     │ \n  │                    ▄▘                           ▚▖                   │ \n  │                   ▞                              ▝▖                  │ \n  │                 ▗▀                                ▝▚                 │ \n  │               ▗▞▘                                   ▀▄               │ \n  │             ▗▞▘                                       ▀▄             │ \n  │          ▗▄▞▘                                           ▀▄▄          │ \n0 │ ▄▄▄▄▄▄▄▞▀▘                                                 ▀▀▄▄▄▄▄▄▄▖│ \n  └──────────────────────────────────────────────────────────────────────┘ \n  0                                                                     70\n\n\nkaiser(n::Integer, α::Real; padding::Integer=0, zerophase::Bool=false)\nkaiser(dims, α; padding=0, zerophase=false)\n\nKaiser window of length n parameterized by α. The Kaiser window approximates the DPSS window (given by dpss), using a simplified definition relying on a Bessel function. Larger values for α give a wider main lobe but have lower sidelobes. Typically α is set around 3.\n\nThe window is defined by sampling the continuous function:\n\n         ⎛  ⎛   _________⎞⎞\n         ⎜  ⎜  ╱        2⎟⎟\nw(x) = I₀⎝πα⎝╲╱ 1 - (2x) ⎠⎠\n       ────────────────────\n               I₀(πα)\n\nin the range [-0.5, 0.5]\n\nWhere I₀(⋅) is the zeroth-order modified Bessel function of the first kind.\n\nProviding a dims Tuple rather than a single n constructs a 2D window. α, padding and zerophase can then be given either as a single value for both horizontal and vertical or a 2-tuple to specify them separately.\n\nIf zerophase is false (the default) the window is centered around index (n+1)/2, which is commonly used for FIR filter design. These are often used for FIR filter design, and usually odd-length. Note that for even-length windows this will cause the window center to be located between samples.\n\nIf zerophase is true the window is centered around index 1 (with the negative half wrapped to the end of the vector). Additionally this creates a \"periodic\" window, which means that if there is no padding then the left and right endpoints of the window wrap around to the same sample, so the window length is the same as an n+1-length non-zerophase window. Alternatively you can think of the continuous zerophase window being of width n and the non-zerophase window being of length n-1. zerophase windows are often used in FFT processing, and are usually even-length.\n\n\n\n\n\n"
},

{
    "location": "windows/#DSP.Windows.dpss",
    "page": "Windows - window functions",
    "title": "DSP.Windows.dpss",
    "category": "function",
    "text": "    ┌──────────────────────────────────────────────────────────────────────┐ \n0.2 │                              ▄▄▀▀▀▀▀▀▚▄▖                             │ \n    │                           ▗▞▀          ▝▀▄                           │ \n    │                         ▗▞▘               ▀▄                         │ \n    │                        ▄▘                   ▚▖                       │ \n    │                      ▗▀                      ▝▖                      │ \n    │                     ▞▘                        ▝▀▖                    │ \n    │                   ▄▀                            ▝▚▖                  │ \n    │                 ▗▞                                ▝▄                 │ \n    │                ▄▘                                   ▚▖               │ \n    │              ▗▞                                      ▝▄              │ \n    │            ▗▞▘                                         ▀▄            │ \n    │          ▗▞▘                                             ▀▄          │ \n    │        ▄▞▘                                                 ▀▄▖       │ \n    │     ▄▞▀                                                      ▝▀▄▖    │ \n  0 │ ▄▞▀▀                                                            ▝▀▀▄▖│ \n    └──────────────────────────────────────────────────────────────────────┘ \n    0                                                                     70\n\n\ndpss(n::Integer, nw::Real, ntapers::Integer=iceil(2*nw)-1;\n     padding::Integer=0, zerophase::Bool=false)\n\nThe first ntapers discrete prolate spheroid sequences (Slepian tapers) as an n × ntapers matrix. The signs of the tapers follow the convention that the first element of the skew-symmetric (odd) tapers is positive. The time-bandwidth product is given by nw.\n\nThe DPSS window maximizes the energy concentration in the main lobe.\n\nIf zerophase is false (the default) the window is centered around index (n+1)/2, which is commonly used for FIR filter design. These are often used for FIR filter design, and usually odd-length. Note that for even-length windows this will cause the window center to be located between samples.\n\nIf zerophase is true the window is centered around index 1 (with the negative half wrapped to the end of the vector). Additionally this creates a \"periodic\" window, which means that if there is no padding then the left and right endpoints of the window wrap around to the same sample, so the window length is the same as an n+1-length non-zerophase window. Alternatively you can think of the continuous zerophase window being of width n and the non-zerophase window being of length n-1. zerophase windows are often used in FFT processing, and are usually even-length.\n\n\n\n\n\n"
},

{
    "location": "windows/#DSP.Windows.dpsseig",
    "page": "Windows - window functions",
    "title": "DSP.Windows.dpsseig",
    "category": "function",
    "text": "dpsseig(A, nw)\n\nEigenvalues of the DPSS matrix, representing the ratios of the power within the main lobe to the total power (main and sidelobes). A is the output of dpss, and nw is the time-bandwidth product provided to dpss as input.\n\n\n\n\n\n"
},

{
    "location": "windows/#Windows-window-functions-1",
    "page": "Windows - window functions",
    "title": "Windows - window functions",
    "category": "section",
    "text": "rect\nhanning\nhamming\ntukey\ncosine\nlanczos\ntriang\nbartlett\ngaussian\nbartlett_hann\nblackman\nkaiser\ndpss\ndpsseig"
},

{
    "location": "filters/#",
    "page": "Filters - filter design and filtering",
    "title": "Filters - filter design and filtering",
    "category": "page",
    "text": ""
},

{
    "location": "filters/#Filters-filter-design-and-filtering-1",
    "page": "Filters - filter design and filtering",
    "title": "Filters - filter design and filtering",
    "category": "section",
    "text": "DSP.jl differentiates between filter coefficients and stateful filters. Filter coefficient objects specify the response of the filter in one of several standard forms. Stateful filter objects carry the state of the filter together with filter coefficients in an implementable form (PolynomialRatio, Biquad, or SecondOrderSections). When invoked on a filter coefficient object, filt does not preserve state."
},

{
    "location": "filters/#DSP.Filters.ZeroPoleGain",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.ZeroPoleGain",
    "category": "type",
    "text": "ZeroPoleGain(z, p, k)\n\nFilter representation in terms of zeros z, poles p, and gain k:\n\nH(x) = kfrac(x - verbz1) ldots (x - verbzend)(x - verbp1) ldots (x - verbpend)\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.PolynomialRatio",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.PolynomialRatio",
    "category": "type",
    "text": "PolynomialRatio(b, a)\n\nFilter representation in terms of the coefficients of the numerator b and denominator a of the transfer function:\n\nH(s) = fracverbb1 s^n-1 + ldots + verbbnverba1 s^n-1 + ldots + verban\n\nor equivalently:\n\nH(z) = fracverbb1 + ldots + verbbn z^-n+1verba1 + ldots + verban z^-n+1\n\nb and a may be specified as Polynomial objects or vectors ordered from highest power to lowest.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.Biquad",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.Biquad",
    "category": "type",
    "text": "Biquad(b0, b1, b2, a1, a2)\n\nFilter representation in terms of the transfer function of a single second-order section given by:\n\nH(s) = fracverbb0 s^2+verbb1 s+verbb2s^2+verba1 s + verba2\n\nor equivalently:\n\nH(z) = fracverbb0+verbb1 z^-1+verbb2 z^-21+verba1 z^-1 + verba2 z^-2\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.SecondOrderSections",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.SecondOrderSections",
    "category": "type",
    "text": "SecondOrderSections(biquads, gain)\n\nFilter representation in terms of a cascade of second-order sections and gain. biquads must be specified as a vector of Biquads.\n\n\n\n\n\n"
},

{
    "location": "filters/#coefficient-objects-1",
    "page": "Filters - filter design and filtering",
    "title": "Filter coefficient objects",
    "category": "section",
    "text": "DSP.jl supports common filter representations. Filter coefficients can be converted from one type to another using convert.ZeroPoleGain\nPolynomialRatio\nBiquad\nSecondOrderSections"
},

{
    "location": "filters/#DSP.Filters.DF2TFilter",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.DF2TFilter",
    "category": "type",
    "text": "DF2TFilter(coef[, si])\n\nConstruct a stateful direct form II transposed filter with coefficients coef. si is an optional array representing the initial filter state (defaults to zeros). If f is a PolynomialRatio, Biquad, or SecondOrderSections, filtering is implemented directly. If f is a ZeroPoleGain object, it is first converted to a SecondOrderSections object.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.FIRFilter",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.FIRFilter",
    "category": "type",
    "text": "FIRFilter(h[, ratio])\n\nConstruct a stateful FIRFilter object from the vector of filter taps h. ratio is an optional rational integer which specifies the input to output sample rate relationship (e.g. 147//160 for converting recorded audio from 48 KHz to 44.1 KHz).\n\n\n\n\n\nFIRFilter(h, rate[, Nϕ])\n\nReturns a polyphase FIRFilter object from the vector of filter taps h. rate is a floating point number that specifies the input to output sample-rate relationship fracfs_outfs_in. Nϕ is an optional parameter which specifies the number of phases created from h. Nϕ defaults to 32.\n\n\n\n\n\n"
},

{
    "location": "filters/#stateful-filter-objects-1",
    "page": "Filters - filter design and filtering",
    "title": "Stateful filter objects",
    "category": "section",
    "text": "DF2TFilterDSP.jl\'s FIRFilter type maintains state between calls to filt, allowing you to filter a signal of indefinite length in RAM-friendly chunks. FIRFilter contains nothing more that the state of the filter, and a FIRKernel. There are five different kinds of FIRKernel for single rate, up-sampling, down-sampling, rational resampling, and arbitrary sample-rate conversion. You need not specify the type of kernel. The FIRFilter constructor selects the correct kernel based on input parameters.FIRFilter"
},

{
    "location": "filters/#DSP.filt",
    "page": "Filters - filter design and filtering",
    "title": "DSP.filt",
    "category": "function",
    "text": "filt(f, x[, si])\n\nApply filter or filter coefficients f along the first dimension of array x. If f is a filter coefficient object, si is an optional array representing the initial filter state (defaults to zeros). If f is a PolynomialRatio, Biquad, or SecondOrderSections, filtering is implemented directly. If f is a ZeroPoleGain object, it is first converted to a SecondOrderSections object.  If f is a Vector, it is interpreted as an FIR filter, and a naïve or FFT-based algorithm is selected based on the data and filter length.\n\n\n\n\n\nfilt(b, a, x, [si])\n\nApply filter described by vectors a and b to vector x, with an optional initial filter state vector si (defaults to zeros).\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.filt!",
    "page": "Filters - filter design and filtering",
    "title": "DSP.filt!",
    "category": "function",
    "text": "filt!(out, f, x[, si])\n\nSame as filt() but writes the result into the out argument. Output array out may not be an alias of x, i.e. filtering may not be done in place.\n\n\n\n\n\nfilt!(out, b, a, x, [si])\n\nSame as filt but writes the result into the out argument, which may alias the input x to modify it in-place.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.filtfilt",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.filtfilt",
    "category": "function",
    "text": "filtfilt(coef, x)\n\nFilter x in the forward and reverse directions using filter coefficients coef. The initial state of the filter is computed so that its response to a step function is steady state. Before filtering, the data is extrapolated at both ends with an odd-symmetric extension of length 3*(max(length(b), length(a))-1).\n\nBecause filtfilt applies the given filter twice, the effective filter order is twice the order of coef. The resulting signal has zero phase distortion.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.fftfilt",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.fftfilt",
    "category": "function",
    "text": "fftfilt(h, x)\n\nApply FIR filter taps h along the first dimension of array x using an FFT-based overlap-save algorithm.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.fftfilt!",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.fftfilt!",
    "category": "function",
    "text": "fftfilt!(out, h, x)\n\nLike fftfilt but writes result into out array.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.tdfilt",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.tdfilt",
    "category": "function",
    "text": "tdfilt(h, x)\n\nApply filter or filter coefficients h along the first dimension of array x using a naïve time-domain algorithm\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.tdfilt!",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.tdfilt!",
    "category": "function",
    "text": "tdfilt!(out, h, x)\n\nLike tdfilt, but writes the result into array out. Output array out may not be an alias of x, i.e. filtering may not be done in place.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.resample",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.resample",
    "category": "function",
    "text": "resample(x, rate[, coef])\n\nResample x at rational or arbitrary rate. coef is an optional vector of FIR filter taps. If coef is not provided, the taps will be computed using a Kaiser window.\n\nInternally, resample uses a polyphase FIRFilter object, but performs additional operations to make resampling a signal easier. It compensates for for the FIRFilter\'s delay (ramp-up), and appends zeros to x. The result is that when the input and output signals are plotted on top of each other, they correlate very well, but one signal will have more samples that the other.\n\n\n\n\n\n"
},

{
    "location": "filters/#Filter-application-1",
    "page": "Filters - filter design and filtering",
    "title": "Filter application",
    "category": "section",
    "text": "filt\nfilt!\nfiltfilt\nfftfilt\nfftfilt!\ntdfilt\ntdfilt!\nresample"
},

{
    "location": "filters/#DSP.Filters.analogfilter",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.analogfilter",
    "category": "function",
    "text": "analogfilter(responsetype, designmethod)\n\nConstruct an analog filter. See below for possible response and filter types.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.digitalfilter",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.digitalfilter",
    "category": "function",
    "text": "digitalfilter(responsetype, designmethod)\n\nConstruct a digital filter. See below for possible response and filter types.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.iirnotch",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.iirnotch",
    "category": "function",
    "text": "iirnotch(Wn, bandwidth[; fs])\n\nSecond-order digital IIR notch filter at frequency Wn with bandwidth bandwidth. If fs is not specified, Wn is interpreted as a normalized frequency in half-cycles/sample.\n\n\n\n\n\n"
},

{
    "location": "filters/#Filter-design-1",
    "page": "Filters - filter design and filtering",
    "title": "Filter design",
    "category": "section",
    "text": "Most analog and digital filters are constructed by composing response types, which determine the frequency response of the filter, with design methods, which determine how the filter is constructed.analogfilter\ndigitalfilterFor some filters, the design method inherently implies a response type. Such filters are documented below.iirnotch"
},

{
    "location": "filters/#DSP.Filters.Lowpass",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.Lowpass",
    "category": "type",
    "text": "Lowpass(Wn[; fs])\n\nLow pass filter with cutoff frequency Wn. If fs is not specified, Wn is interpreted as a normalized frequency in half-cycles/sample.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.Highpass",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.Highpass",
    "category": "type",
    "text": "Highpass(Wn[; fs])\n\nHigh pass filter with cutoff frequency Wn. If fs is not specified, Wn is interpreted as a normalized frequency in half-cycles/sample.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.Bandpass",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.Bandpass",
    "category": "type",
    "text": "Bandpass(Wn1, Wn2[; fs])\n\nBand pass filter with normalized pass band (Wn1, Wn2). If fs is not specified, Wn1 and Wn2 are interpreted as normalized frequencies in half-cycles/sample.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.Bandstop",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.Bandstop",
    "category": "type",
    "text": "Bandstop(Wn1, Wn2[; fs])\n\nBand stop filter with normalized stop band (Wn1, Wn2). If fs is not specified, Wn1 and Wn2 are interpreted as normalized frequencies in half-cycles/sample.\n\n\n\n\n\n"
},

{
    "location": "filters/#response-types-1",
    "page": "Filters - filter design and filtering",
    "title": "Filter response types",
    "category": "section",
    "text": "Lowpass\nHighpass\nBandpass\nBandstop"
},

{
    "location": "filters/#design-methods-1",
    "page": "Filters - filter design and filtering",
    "title": "Filter design methods",
    "category": "section",
    "text": ""
},

{
    "location": "filters/#DSP.Filters.Butterworth",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.Butterworth",
    "category": "function",
    "text": "Butterworth(n)\n\nn pole Butterworth filter.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.Chebyshev1",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.Chebyshev1",
    "category": "function",
    "text": "Chebyshev1(n, ripple)\n\nn pole Chebyshev type I filter with ripple dB ripple in the passband.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.Chebyshev2",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.Chebyshev2",
    "category": "function",
    "text": "Chebyshev2(n, ripple)\n\nn pole Chebyshev type II filter with ripple dB ripple in the stopband.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.Elliptic",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.Elliptic",
    "category": "function",
    "text": "Elliptic(n, rp, rs)\n\nn pole elliptic (Cauer) filter with rp dB ripple in the passband and rs dB attentuation in the stopband.\n\n\n\n\n\n"
},

{
    "location": "filters/#IIR-filter-design-methods-1",
    "page": "Filters - filter design and filtering",
    "title": "IIR filter design methods",
    "category": "section",
    "text": "Butterworth\nChebyshev1\nChebyshev2\nElliptic"
},

{
    "location": "filters/#DSP.Filters.FIRWindow",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.FIRWindow",
    "category": "type",
    "text": "FIRWindow(window; scale=true)\n\nFIR filter design using window window, a vector whose length matches the number of taps in the resulting filter.\n\nIf scale is true (default), the designed FIR filter is scaled so that the following holds:\n\nFor Lowpass and Bandstop filters, the frequency response is unity at 0 (DC).\nFor Highpass filters, the frequency response is unity at the Nyquist frequency.\nFor Bandpass filters, the frequency response is unity in the center of the passband.\n\n\n\n\n\nFIRWindow(; transitionwidth, attenuation=60, scale=true)\n\nKaiser window FIR filter design. The required number of taps is calculated based on transitionwidth (in half-cycles/sample) and stopband attenuation (in dB). attenuation defaults to 60 dB.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.remez",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.remez",
    "category": "function",
    "text": "remez(numtaps::Integer, band_defs;\n      Hz::Real=1.0,\n      neg::Bool=false,\n      maxiter::Integer=25,\n      grid_density::Integer=16)\n\nCalculate the minimax optimal filter using the Remez exchange algorithm [McClellan1973a] [McClellan1973b].\n\nThis is the simplified API that accepts just 2 required arguments (numtaps, band_defs). For a scipy compatible version see the 3 arguments version (numtaps, bands, desired). \n\nCalculate the filter-coefficients for the finite impulse response (FIR) filter whose transfer function minimizes the maximum error between the desired gain and the realized gain in the specified frequency bands using the Remez exchange algorithm.\n\nArguments\n\nnumtaps::Integer: The desired number of taps in the filter.    The number of taps is the number of terms in the filter, or the filter    order plus one.\nbands_defs: A sequence of band definitions.   This sequence defines the bands. Each entry is a pair. The pair\'s   first item is a tuple of band edges (low, high). The pair\'s second item    defines the desired response and weight in that band. The weight is optional   and defaults to 1.0. Both the desired response and weight may be either scalars   or functions. If a function, the function should accept a real frequency and   return the real desired response or real weight. Examples:\nLPF with unity weights. [(0, 0.475) => 1, (0.5, 1.0) => 0]\nLPF with weight of 2 in the stop band. [(0, 0.475) => (1, 1), (0.5, 1.0) => (0, 2)]\nBPF with unity weights. [(0, 0.375) => 0, (0.4, 0.5) => 1, (0.525, 1.0) => 0]\nHilbert transformer. [(0.1, 0.95) => 1]; neg=true\nDifferentiator. [(0.01, 0.99) => (f -> f/2, f -> 1/f)]; neg=true\nHz::Real: The sampling frequency in Hz. Default is 1.\nneg::Bool: Whether the filter has negative symmetry or not. Default is false.   If false, the filter is even-symmetric. If true, the filter is odd-symmetric.   neg=true means that h[n]=-h[end+1-n]; neg=false means that h[n]=h[end+1-n].\nmaxiter::Integer: (optional)   Maximum number of iterations of the algorithm. Default is 25.\ngrid_density:Integer: (optional)   Grid density. The dense grid used in remez is of size   (numtaps + 1) * grid_density. Default is 16.\n\nReturns\n\nh::Array{Float64,1}: A rank-1 array containing the coefficients of the optimal   (in a minimax sense) filter.\n\n[McClellan1973a]: \n\nJ. H. McClellan and T. W. Parks, A unified approach to the design of optimum FIR linear phase digital filters, IEEE Trans. Circuit Theory, vol. CT-20, pp. 697-701, 1973.\n\n[McClellan1973b]: \n\nJ. H. McClellan, T. W. Parks and L. R. Rabiner, A Computer Program for Designing Optimum FIR Linear Phase Digital Filters, IEEE Trans. Audio Electroacoust., vol. AU-21, pp. 506-525, 1973.\n\nExamples\n\nConstruct a length 35 filter with a passband at 0.15-0.4 Hz  (desired response of 1), and stop bands at 0-0.1 Hz and 0.45-0.5 Hz (desired response of 0). Note: the behavior in the frequency ranges between  those bands - the transition bands - is unspecified.\n\njulia> bpass = remez(35, [(0, 0.1)=>0, (0.15, 0.4)=>1, (0.45, 0.5)=>0])\n\nYou can trade-off maximum error achieved for transition bandwidth.  The wider the transition bands, the lower the maximum error in the bands specified. Here is a bandpass filter with the same passband, but wider transition bands.\n\njulia> bpass2 = remez(35, [(0, 0.08)=>0, (0.15, 0.4)=>1, (0.47, 0.5)=>0])\n\nHere we compute the frequency responses and plot them in dB.\n\nusing PyPlot\nb = DSP.Filters.PolynomialRatio(bpass, [1.0])\nb2 = DSP.Filters.PolynomialRatio(bpass2, [1.0])\nf = range(0, stop=0.5, length=1000)\nplot(f, 20*log10.(abs.(freqz(b,f,1.0))))\nplot(f, 20*log10.(abs.(freqz(b2,f,1.0))))\ngrid()\n\n\n\n\n\nremez(numtaps::Integer, \n      bands::Vector, \n      desired::Vector; \n      weight::Vector=[], \n      Hz::Real=1.0, \n      filter_type::RemezFilterType=filter_type_bandpass,\n      maxiter::Integer=25, \n      grid_density::Integer=16)\n\nThis is the scipy compatible version that requires 3 arguments (numtaps, bands, desired).  For a simplified API, see the 2 argument version (numtaps, band_defs). The filters designed are equivalent, the inputs are just specified in a different way. Below the arguments and examples are described that differ from the simplified API version.\n\nArguments\n\nbands::Vector: A monotonic sequence containing the band edges in Hz.   All elements must be non-negative and less than half the sampling   frequency as given by Hz.\ndesired::Vector:A sequence half the size of bands containing the desired    gain in each of the specified bands.\nweight::Vector: (optional)   A relative weighting to give to each band region. The length of   weight has to be half the length of bands.\nfilter_type::RemezFilterType: Default is filter_type_bandpass.   The type of filter:\nfilter_type_bandpass : flat response in bands. This is the default.\nfilter_type_differentiator : frequency proportional response in bands.   Odd symetric as in filter_type_hilbert case, but with a linear sloping   desired response.\nfilter_type_hilbert : filter with odd symmetry, that is, type III             (for even order) or type IV (for odd order)             linear phase filters.\n\nExamples\n\nCompare the examples with the scipy API and the simplified API.\n\nSimple:\njulia> bpass = remez(35, [(0, 0.1)=>0, (0.15, 0.4)=>1, (0.45, 0.5)=>0])\nScipy:\njulia> bpass = remez(35, [0, 0.1, 0.15, 0.4, 0.45, 0.5], [0, 1, 0])\n\nSimple:\njulia> bpass2 = remez(35, [(0, 0.08)=>0, (0.15, 0.4)=>1, (0.47, 0.5)=>0])\nScipy:\njulia> bpass2 = remez(35, [0, 0.08, 0.15, 0.4, 0.47, 0.5], [0, 1, 0])\n\n\n\n\n\n"
},

{
    "location": "filters/#FIR-filter-design-methods-1",
    "page": "Filters - filter design and filtering",
    "title": "FIR filter design methods",
    "category": "section",
    "text": "FIRWindow\nremez"
},

{
    "location": "filters/#DSP.Filters.freqz",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.freqz",
    "category": "function",
    "text": "freqz(filter, w = range(0, stop=π, length=250))\n\nFrequency response of a digital filter at normalised frequency or frequencies w in radians/sample.\n\n\n\n\n\nfreqz(filter, hz, fs)\n\nFrequency response of a digital filter at frequency or frequencies hz with sampling rate fs.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.phasez",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.phasez",
    "category": "function",
    "text": "phasez(filter, w = range(0, stop=π, length=250))\n\nPhase response of a digital filter at normalised frequency or frequencies w in radians/sample.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.impz",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.impz",
    "category": "function",
    "text": "impz(filter, n=100)\n\nImpulse response of a digital filter with n points.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.stepz",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.stepz",
    "category": "function",
    "text": "stepz(filter, n=100)\n\nStep response of a digital filter with n points.\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.freqs",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.freqs",
    "category": "function",
    "text": "freqs(filter, w)\n\nFrequency response of an analog filter at normalised frequency or frequencies w in radians/sample.\n\n\n\n\n\nfreqs(filter, hz, fs)\n\nFrequency response of an analog filter at frequency or frequencies hz with sampling rate fs.\n\n\n\n\n\n"
},

{
    "location": "filters/#Filter-response-1",
    "page": "Filters - filter design and filtering",
    "title": "Filter response",
    "category": "section",
    "text": "freqz\nphasez\nimpz\nstepz\nfreqs"
},

{
    "location": "filters/#DSP.Filters.coefb",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.coefb",
    "category": "function",
    "text": "coefb(f)\n\nCoefficients of the numerator of a PolynomialRatio object, highest power first, i.e., the b passed to filt()\n\n\n\n\n\n"
},

{
    "location": "filters/#DSP.Filters.coefa",
    "page": "Filters - filter design and filtering",
    "title": "DSP.Filters.coefa",
    "category": "function",
    "text": "coefa(f)\n\nCoefficients of the denominator of a PolynomialRatio object, highest power first, i.e., the a passed to filt()\n\n\n\n\n\n"
},

{
    "location": "filters/#Miscellaneous-1",
    "page": "Filters - filter design and filtering",
    "title": "Miscellaneous",
    "category": "section",
    "text": "coefb\ncoefa"
},

{
    "location": "filters/#Examples-1",
    "page": "Filters - filter design and filtering",
    "title": "Examples",
    "category": "section",
    "text": "Construct a 4th order elliptic lowpass filter with normalized cutoff frequency 0.2, 0.5 dB of passband ripple, and 30 dB attentuation in the stopband and extract the coefficients of the numerator and denominator of the transfer function:responsetype = Lowpass(0.2)\ndesignmethod = Elliptic(4, 0.5, 30)\ntf = convert(PolynomialRatio, digitalfilter(responsetype, designmethod))\nnumerator_coefs = coefb(tf)\ndenominator_coefs = coefa(tf)Filter the data in x, sampled at 1000 Hz, with a 4th order Butterworth bandpass filter between 10 and 40 Hz:responsetype = Bandpass(10, 40; fs=1000)\ndesignmethod = Butterworth(4)\nfilt(digitalfilter(responsetype, designmethod), x)Filter the data in x, sampled at 50 Hz, with a 64 tap Hanning window FIR lowpass filter at 5 Hz:responsetype = Lowpass(5; fs=50)\ndesignmethod = FIRWindow(hanning(64))\nfilt(digitalfilter(responsetype, designmethod), x)"
},

{
    "location": "util/#",
    "page": "Util - utility functions",
    "title": "Util - utility functions",
    "category": "page",
    "text": ""
},

{
    "location": "util/#DSP.Unwrap.unwrap",
    "page": "Util - utility functions",
    "title": "DSP.Unwrap.unwrap",
    "category": "function",
    "text": "unwrap(m; kwargs...)\n\nAssumes m to be a sequence of values that has been wrapped to be inside the given range (centered around zero), and undoes the wrapping by identifying discontinuities. If a single dimension is passed to dims, then m is assumed to have wrapping discontinuities only along that dimension. If a range of dimensions, as in 1:ndims(m), is passed to dims, then m is assumed to have wrapping discontinuities across all ndims(m) dimensions.\n\nA common usage for unwrapping across a singleton dimension is for a phase measurement over time, such as when comparing successive frames of a short-time-fourier-transform, as each frame is wrapped to stay within (-pi, pi].\n\nA common usage for unwrapping across multiple dimensions is for a phase measurement of a scene, such as when retrieving the phase information of of an image, as each pixel is wrapped to stay within (-pi, pi].\n\nArguments\n\nm::AbstractArray{T, N}: Array to unwrap.\ndims=nothing: Dimensions along which to unwrap. If dims is an integer, then   unwrap is called on that dimension. If dims=1:ndims(m), then m is unwrapped   across all dimensions.\nrange=2pi: Range of wrapped array.\ncircular_dims=(false, ...):  When an element of this tuple is true, the   unwrapping process will consider the edges along the corresponding axis   of the array to be connected.\nrng=GLOBAL_RNG: Unwrapping of arrays with dimension > 1 uses a random   initialization. A user can pass their own RNG through this argument.\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Unwrap.unwrap!",
    "page": "Util - utility functions",
    "title": "DSP.Unwrap.unwrap!",
    "category": "function",
    "text": "unwrap!(m; kwargs...)\n\nIn-place version of unwrap.\n\n\n\n\n\nunwrap!(y, m; kwargs...)\n\nUnwrap m storing the result in y, see unwrap.\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Util.hilbert",
    "page": "Util - utility functions",
    "title": "DSP.Util.hilbert",
    "category": "function",
    "text": "hilbert(x)\n\nComputes the analytic representation of x, x_a = x + j hatx, where hatx is the Hilbert transform of x, along the first dimension of x.\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Util.fftfreq",
    "page": "Util - utility functions",
    "title": "DSP.Util.fftfreq",
    "category": "function",
    "text": "fftfreq(n, fs=1)\n\nReturn discrete fourier transform sample frequencies. The returned Frequencies object is an AbstractVector containing the frequency bin centers at every sample point. fs is the sample rate of the input signal.\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Util.rfftfreq",
    "page": "Util - utility functions",
    "title": "DSP.Util.rfftfreq",
    "category": "function",
    "text": "rfftfreq(n, fs=1)\n\nReturn discrete fourier transform sample frequencies for use with rfft. The returned Frequencies object is an AbstractVector containing the frequency bin centers at every sample point. fs is the sample rate of the input signal.\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Util.nextfastfft",
    "page": "Util - utility functions",
    "title": "DSP.Util.nextfastfft",
    "category": "function",
    "text": "nextfastfft(n)\n\nReturn the closest product of 2, 3, 5, and 7 greater than or equal to n. FFTW contains optimized kernels for these sizes and computes Fourier transforms of input that is a product of these sizes faster than for input of other sizes.\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Util.pow2db",
    "page": "Util - utility functions",
    "title": "DSP.Util.pow2db",
    "category": "function",
    "text": "pow2db(a)\n\nConvert a power ratio to dB (decibel), or 10log_10(a). The inverse of db2pow.\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Util.amp2db",
    "page": "Util - utility functions",
    "title": "DSP.Util.amp2db",
    "category": "function",
    "text": "amp2db(a)\n\nConvert an amplitude ratio to dB (decibel), or 20 log_10(a)=10log_10(a^2). The inverse of db2amp.\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Util.db2pow",
    "page": "Util - utility functions",
    "title": "DSP.Util.db2pow",
    "category": "function",
    "text": "db2pow(a)\n\nConvert dB to a power ratio. This function call also be called using a*dB, i.e. 3dB == db2pow(3). The inverse of pow2db.\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Util.db2amp",
    "page": "Util - utility functions",
    "title": "DSP.Util.db2amp",
    "category": "function",
    "text": "db2amp(a)\n\nConvert dB to an amplitude ratio. This function call also be called using a*dBa, i.e. 3dBa == db2amp(3). The inverse of amp2db.\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Util.rms",
    "page": "Util - utility functions",
    "title": "DSP.Util.rms",
    "category": "function",
    "text": "rms(s)\n\nReturn the root mean square of signal s.\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Util.rmsfft",
    "page": "Util - utility functions",
    "title": "DSP.Util.rmsfft",
    "category": "function",
    "text": "rmsfft(f)\n\nReturn the root mean square of signal s given the FFT transform f = fft(s). Equivalent to rms(ifft(f)).\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Util.meanfreq",
    "page": "Util - utility functions",
    "title": "DSP.Util.meanfreq",
    "category": "function",
    "text": "meanfreq(x, fs)\n\nCalculate the mean power frequency of x with a sampling frequency of fs, defined as:\n\nMPF = fracsum_i=1^F f_i X_i^2 sum_i=0^F X_i^2  Hz\n\nwhere F is the Nyquist frequency, and X is the power spectral density.\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Util.finddelay",
    "page": "Util - utility functions",
    "title": "DSP.Util.finddelay",
    "category": "function",
    "text": "finddelay(x, y)\n\nEstimate the delay of x with respect to y by locating the peak of their cross-correlation.\n\nThe output delay will be positive when x is delayed with respect y, negative if advanced, 0 otherwise.\n\nExample\n\njulia> finddelay([0, 0, 1, 2, 3], [1, 2, 3])\n2\n\njulia> finddelay([1, 2, 3], [0, 0, 1, 2, 3])\n-2\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Util.shiftsignal",
    "page": "Util - utility functions",
    "title": "DSP.Util.shiftsignal",
    "category": "function",
    "text": "shiftsignal(x, s)\n\nShift elements of signal x in time by a given amount s of samples and fill the spaces with zeros. For circular shifting, use circshift.\n\nExample\n\njulia> shiftsignal([1, 2, 3], 2)\n3-element Array{Int64,1}:\n 0\n 0\n 1\n\njulia> shiftsignal([1, 2, 3], -2)\n3-element Array{Int64,1}:\n 3\n 0\n 0\n\nSee also shiftsignal!.\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Util.shiftsignal!",
    "page": "Util - utility functions",
    "title": "DSP.Util.shiftsignal!",
    "category": "function",
    "text": "shiftsignal!(x, s)\n\nMutating version of shiftsignals(): shift x of s samples and fill the spaces with zeros in-place.\n\nSee also shiftsignal.\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Util.alignsignals",
    "page": "Util - utility functions",
    "title": "DSP.Util.alignsignals",
    "category": "function",
    "text": "alignsignals(x, y)\n\nUse finddelay() and shiftsignal() to time align x to y. Also return the delay of x with respect to y.\n\nExample\n\njulia> alignsignals([0, 0, 1, 2, 3], [1, 2, 3])\n([1, 2, 3, 0, 0], 2)\n\njulia> alignsignals([1, 2, 3], [0, 0, 1, 2, 3])\n([0, 0, 1], -2)\n\nSee also alignsignals!.\n\n\n\n\n\n"
},

{
    "location": "util/#DSP.Util.alignsignals!",
    "page": "Util - utility functions",
    "title": "DSP.Util.alignsignals!",
    "category": "function",
    "text": "alignsignals!(x, y)\n\nMutating version of alignsignals(): time align x to y in-place.\n\nSee also alignsignals.\n\n\n\n\n\n"
},

{
    "location": "util/#Util-utility-functions-1",
    "page": "Util - utility functions",
    "title": "Util - utility functions",
    "category": "section",
    "text": "DocTestSetup = quote\n    using DSP\nendunwrap\nunwrap!\nhilbert\nfftfreq\nrfftfreq\nnextfastfft\npow2db\namp2db\ndb2pow\ndb2amp\nrms\nrmsfft\nmeanfreq\nfinddelay\nshiftsignal\nshiftsignal!\nalignsignals\nalignsignals!DocTestSetup = nothing"
},

{
    "location": "convolutions/#",
    "page": "Convolutions - similarity methods",
    "title": "Convolutions - similarity methods",
    "category": "page",
    "text": ""
},

{
    "location": "convolutions/#DSP.conv",
    "page": "Convolutions - similarity methods",
    "title": "DSP.conv",
    "category": "function",
    "text": "conv(u,v)\n\nConvolution of two vectors. Uses FFT algorithm.\n\n\n\n\n\n"
},

{
    "location": "convolutions/#DSP.conv2",
    "page": "Convolutions - similarity methods",
    "title": "DSP.conv2",
    "category": "function",
    "text": "conv2(u,v,A)\n\n2-D convolution of the matrix A with the 2-D separable kernel generated by the vectors u and v. Uses 2-D FFT algorithm.\n\n\n\n\n\nconv2(B,A)\n\n2-D convolution of the matrix B with the matrix A. Uses 2-D FFT algorithm.\n\n\n\n\n\n"
},

{
    "location": "convolutions/#DSP.deconv",
    "page": "Convolutions - similarity methods",
    "title": "DSP.deconv",
    "category": "function",
    "text": "deconv(b,a) -> c\n\nConstruct vector c such that b = conv(a,c) + r. Equivalent to polynomial division.\n\n\n\n\n\n"
},

{
    "location": "convolutions/#DSP.xcorr",
    "page": "Convolutions - similarity methods",
    "title": "DSP.xcorr",
    "category": "function",
    "text": "xcorr(u,v)\n\nCompute the cross-correlation of two vectors.\n\n\n\n\n\n"
},

{
    "location": "convolutions/#Convolutions-similarity-methods-1",
    "page": "Convolutions - similarity methods",
    "title": "Convolutions - similarity methods",
    "category": "section",
    "text": "conv\nconv2\ndeconv\nxcorr"
},

{
    "location": "lpc/#",
    "page": "LPC - Linear Predictive Coding",
    "title": "LPC - Linear Predictive Coding",
    "category": "page",
    "text": ""
},

{
    "location": "lpc/#DSP.LPC.lpc",
    "page": "LPC - Linear Predictive Coding",
    "title": "DSP.LPC.lpc",
    "category": "function",
    "text": "lpc(x::AbstractVector, p::Int, [LPCBurg()])\n\nGiven input signal x and prediction order p, returns IIR coefficients a and average reconstruction error prediction_err. Note that this method does NOT return the leading 1 present in the true autocorrelative estimate; it omits it as it is implicit in every LPC estimate, and must be manually reintroduced if the returned vector should be treated as a polynomial.\n\nThe algorithm used is determined by the last optional parameter, and can be either LPCBurg or LPCLevinson.\n\n\n\n\n\n"
},

{
    "location": "lpc/#DSP.LPC.lpc-Tuple{AbstractArray{Number,1},Int64,LPCBurg}",
    "page": "LPC - Linear Predictive Coding",
    "title": "DSP.LPC.lpc",
    "category": "method",
    "text": "lpc(x::AbstractVector, p::Int, LPCBurg())\n\nLPC (Linear-Predictive-Code) estimation, using the Burg method. This function implements the mathematics published in [1].\n\n[1] - Enhanced Partial Tracking Using Linear Prediction (DAFX 2003 article, Lagrange et al) http://www.sylvain-marchand.info/Publications/dafx03.pdf\n\n\n\n\n\n"
},

{
    "location": "lpc/#DSP.LPC.lpc-Tuple{AbstractArray{Number,1},Int64,LPCLevinson}",
    "page": "LPC - Linear Predictive Coding",
    "title": "DSP.LPC.lpc",
    "category": "method",
    "text": "lpc(x::AbstractVector, p::Int, LPCLevinson())\n\nLPC (Linear-Predictive-Code) estimation, using the Levinson method. This function implements the mathematics described in [1].\n\n[1] - The Wiener (RMS) Error Criterion in Filter Design and Prediction (Studies in Applied Mathematics 1947 article, N. Levison)\n\n\n\n\n\n"
},

{
    "location": "lpc/#LPC-Linear-Predictive-Coding-1",
    "page": "LPC - Linear Predictive Coding",
    "title": "LPC - Linear Predictive Coding",
    "category": "section",
    "text": "lpc\nlpc(::AbstractVector{Number}, ::Int, ::LPCBurg)\nlpc(::AbstractVector{Number}, ::Int, ::LPCLevinson)"
},

{
    "location": "#",
    "page": "Index",
    "title": "Index",
    "category": "page",
    "text": ""
},

{
    "location": "#Index-1",
    "page": "Index",
    "title": "Index",
    "category": "section",
    "text": ""
},

]}
