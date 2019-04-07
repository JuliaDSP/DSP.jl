# Awaiting next release

- Addition of `lpc_burg` and `lpc_levinson` methods to calculate the Linear
  Predictive Coefficients describing a signal.

# DSP v0.0.8 Release Notes

- The `DF2TFilter` object provides a filter that preserves state
  between invocations ([#100](https://github.com/JuliaDSP/DSP.jl/pull/100)).

# DSP v0.0.7 Release Notes

- Filter coefficient types have been renamed to distinguish them from implementations ([#96](https://github.com/JuliaDSP/DSP.jl/pull/96)):
  - The `Filter` abstract type is now `FilterCoefficients`
  - `ZPKFilter` is now `ZeroPoleGain`
  - `TFFilter` is now `PolynomialRatio`
  - `BiquadFilter` is now `Biquad`
  - `SOSFilter` is now `SecondOrderSections`
