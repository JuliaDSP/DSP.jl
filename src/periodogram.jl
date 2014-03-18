# The periodogram module contains functions which compute non-parametric
# estimates of the periodogram P[s] of a signal s.  An overview of some 
# of the methods is available at:
# http://www.ee.lamar.edu/gleb/adsp/Lecture%2008%20-%20Nonparametric%20SE.pdf
module Periodogram
export arraysplit, periodogram, welch_pgram, bartlett_pgram, spectrogram

# Split an array into subarrays of length N, with overlapping regions
# of length M.
function arraysplit(s, n::Integer, m::Integer)
    # n = m is a problem - the algorithm will not terminate.
    if !(0 <= m < n)
        error("m must be between zero and n.")
    end

    # the length of the non-overlapping array stride
    l = n - m
    
    # total number of strides is the total length of the signal divided
    # by the unique number of elements per stride.  extra elements at the
    # end of of the signal are dropped.
    k = int(length(s)/l - n/l + 1)
    [s[(a*l + 1):(a*l + n)] for a=0:(k-1)]
end

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

function spectrogram(x; nfft=int(length(x)/8), fs=1, noverlap=int(nfft/2))
  p=[periodogram(x) for x in arraysplit(x, nfft, noverlap)]
  p=hcat(p...)
  t=[1:size(p,2)]*(nfft-noverlap)/fs
  f=[1:size(p,1)]/size(p,1)*(fs/2)
  p, t, f
end

end # end module definition
