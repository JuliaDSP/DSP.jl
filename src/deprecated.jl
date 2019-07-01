# deprecations after 0.5
Base.@deprecate_binding ZPKFilter ZeroPoleGain
Base.@deprecate_binding TFFilter PolynomialRatio
Base.@deprecate_binding BiquadFilter Biquad
Base.@deprecate_binding SOSFilter SecondOrderSections

Base.@deprecate Frequencies(nreal::Int, n::Int, multiplier::Float64) FFTW.Frequencies(nreal, n, multiplier)
Base.@deprecate fftfreq(n::Int, fs::Real=1) FFTW.fftfreq(n, fs)
Base.@deprecate rfftfreq(n::Int, fs::Real=1) FFTW.rfftfreq(n, fs)

@deprecate conv2 conv
