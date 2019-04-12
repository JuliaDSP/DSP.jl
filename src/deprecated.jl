Base.@deprecate_binding ZPKFilter ZeroPoleGain
Base.@deprecate_binding TFFilter PolynomialRatio
Base.@deprecate_binding BiquadFilter Biquad
Base.@deprecate_binding SOSFilter SecondOrderSections

@deprecate firfilt(h, x) filt(h, x)
@deprecate conv2 conv
