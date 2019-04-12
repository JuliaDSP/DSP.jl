# deprecations after 0.5
Base.@deprecate_binding ZPKFilter ZeroPoleGain
Base.@deprecate_binding TFFilter PolynomialRatio
Base.@deprecate_binding BiquadFilter Biquad
Base.@deprecate_binding SOSFilter SecondOrderSections

@deprecate conv2 conv
