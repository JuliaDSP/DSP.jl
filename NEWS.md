# DSP v0.0.7 Release Notes

- Filter coefficient types have been renamed to distinguish them from implementations ([#96](https://github.com/JuliaDSP/DSP.jl/pull/96)):
  - The `Filter` abstract type is now `FilterCoefficients`
  - `ZPKFilter` is now `ZeroPoleGain`
  - `TFFilter` is now `PolynomialRatio`
  - `BiquadFilter` is now `Biquad`
  - `SOSFilter` is now `SecondOrderSections`
