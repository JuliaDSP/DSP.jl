const ZPKFilter = ZeroPoleGain
const TFFilter = PolynomialRatio
const BiquadFilter = Biquad
const SOSFilter = SecondOrderSections
export ZPKFilter, TFFilter, BiquadFilter, SOSFilter

@deprecate firfilt(h, x) filt(h, x)
@deprecate conv2 conv
