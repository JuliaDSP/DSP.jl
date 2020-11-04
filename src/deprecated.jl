# deprecations after 0.6
@deprecate freqz(filter::FilterCoefficients{:z}) freqresp(filter)
@deprecate freqz(filter::FilterCoefficients{:z}, w) freqresp(filter, w)
@deprecate freqs(filter::FilterCoefficients{:s}, w) freqresp(filter, w)
@deprecate freqz(filter::FilterCoefficients{:z}, hz, fs) freqresp(filter, hz, fs)
@deprecate freqs(filter::FilterCoefficients{:s}, hz, fs) freqresp(filter, hz, fs)
@deprecate phasez(filter::FilterCoefficients{:z}) phaseresp(filter)
@deprecate phasez(filter::FilterCoefficients{:z}, w) phaseresp(filter, w)
