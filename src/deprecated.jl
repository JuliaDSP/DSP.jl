# deprecations in 0.8
import .Util.nextfastfft
@deprecate nextfastfft(ns::Tuple) nextfastfft.(ns) false
@deprecate nextfastfft(ns...) nextfastfft.(ns) false

# deprecations in 0.7
@deprecate freqz(filter::FilterCoefficients{:z}) freqresp(filter, range(0, stop=π, length=250))
@deprecate freqz(filter::FilterCoefficients{:z}, w) freqresp(filter, w)
@deprecate freqs(filter::FilterCoefficients{:s}, w) freqresp(filter, w)
@deprecate freqz(filter::FilterCoefficients{:z}, hz, fs) freqresp(filter,  hz * ((2 * pi) / fs))
@deprecate freqs(filter::FilterCoefficients{:s}, hz, fs) freqresp(filter,  hz * ((2 * pi) / fs))
@deprecate phasez(filter::FilterCoefficients{:z}) phaseresp(filter, range(0, stop=π, length=250))
@deprecate phasez(filter::FilterCoefficients{:z}, w) phaseresp(filter, w)
@deprecate grpdelayz(filter::FilterCoefficients{:z}) grpdelay(filter, range(0, stop=π, length=250))
@deprecate grpdelayz(filter::FilterCoefficients{:z}, w) grpdelay(filter, w)
@deprecate impz(filter::FilterCoefficients{:z}) impresp(filter)
@deprecate impz(filter::FilterCoefficients{:z}, n) impresp(filter, n)
@deprecate stepz(filter::FilterCoefficients{:z}) stepresp(filter)
@deprecate stepz(filter::FilterCoefficients{:z}, n) stepresp(filter, n)
