# The periodogram module contains functions which compute non-parametric
# estimates of the periodogram P[s] of a signal s.  An overview of some 
# of the methods is available at:
# http://www.ee.lamar.edu/gleb/adsp/Lecture%2008%20-%20Nonparametric%20SE.pdf
module Periodogram

using DSP
using Base

export periodogram, welch_pgram, bartlett_pgram

# Compute the periodogram of a signal S, defined as 1/N*X[s(n)]^2, where X is the
# DTFT of the signal S.
function periodogram(s)
    s_fft = fft(s)
    1/length(s_fft)*real((conj(s_fft) .* s_fft))
end

# Compute an estimate of the power spectral density of a signal s via Welch's
# method.  The resulting periodogram has length N and is computed with an overlap
# region of length M.  The method is detailed in "The Use of Fast Fourier Transform
# for the Estimation of Power Spectra: A Method based on Time Averaging over Short,
# Modified Periodograms."  P. Welch, IEEE Transactions on Audio and Electroacoustics,
# vol AU-15, pp 70-73, 1967.
function welch_pgram(s, n, m)
    sig_split = arraysplit(s, n, m)
    1/length(sig_split)*sum([periodogram(x) for x in sig_split])
end

# Compute an estimate of the periodogram of a signal s via Bartlett's method.  
# The resulting periodogram has length N.  The method appears in "Smoothing 
# Periodograms from Time Series with Continuous Spectra", M.S. Bartlett, Nature 
# #4096, May 1, 1948.The estimate is equivalent to welch_pgram(s, n, 0), as 
# it is a special case of the Welch estimate of the periodogram.
function bartlett_pgram(s, n)
    welch_pgram(s, n, 0)
end

end # end module definition