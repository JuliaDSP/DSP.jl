# `Filters` - filter design and filtering

DSP.jl differentiates between [filter coefficients](@ref coefficient-objects)
and [stateful filters](@ref stateful-filter-objects). Filter
coefficient objects specify the response of the filter in one of
several standard forms. Stateful filter objects carry the state of the
filter together with filter coefficients in an implementable form
(`PolynomialRatio`, `Biquad`, or `SecondOrderSections`).
When invoked on a filter coefficient object, `filt` does not preserve
state.

## [Filter coefficient objects](@id coefficient-objects)

DSP.jl supports common filter representations. Filter coefficients can
be converted from one type to another using `convert`.

```@docs
ZeroPoleGain
PolynomialRatio
Biquad
SecondOrderSections
```

These filter coefficient objects support the following arithmetic operations:
inversion (`inv`), multiplication (`*`) for series connection, and integral
power (`^`) for repeated multiplication with itself. For example:

```jldoctest; setup = :(using DSP)
julia> H = PolynomialRatio([1.0], [1.0, 0.3])
PolynomialRatio{:z, Float64}(LaurentPolynomial(1.0), LaurentPolynomial(0.3*z⁻¹ + 1.0))

julia> inv(H)
PolynomialRatio{:z, Float64}(LaurentPolynomial(0.3*z⁻¹ + 1.0), LaurentPolynomial(1.0))

julia> H * H
PolynomialRatio{:z, Float64}(LaurentPolynomial(1.0), LaurentPolynomial(0.09*z⁻² + 0.6*z⁻¹ + 1.0))

julia> H^2
PolynomialRatio{:z, Float64}(LaurentPolynomial(1.0), LaurentPolynomial(0.09*z⁻² + 0.6*z⁻¹ + 1.0))

julia> H^-2
PolynomialRatio{:z, Float64}(LaurentPolynomial(0.09*z⁻² + 0.6*z⁻¹ + 1.0), LaurentPolynomial(1.0))

```

## [Stateful filter objects](@id stateful-filter-objects)

```@docs
DF2TFilter
```

DSP.jl's `FIRFilter` type maintains state between calls to [`filt`](@ref), allowing
you to filter a signal of indefinite length in RAM-friendly chunks. `FIRFilter`
contains nothing more that the state of the filter, and a `FIRKernel`. There are
five different kinds of `FIRKernel` for single rate, up-sampling, down-sampling,
rational resampling, and arbitrary sample-rate conversion. You need not specify the
type of kernel. The `FIRFilter` constructor selects the correct kernel based on input
parameters.

```@docs
FIRFilter
```

## Filter application

```@docs
filt
filt!
filtfilt
fftfilt
fftfilt!
tdfilt
tdfilt!
resample
```

## Filter design

Most analog and digital filters are constructed by composing
[response types](@ref response-types), which determine the frequency
response of the filter, with [design methods](@ref design-methods),
which determine how the filter is constructed.

The response type is [`Lowpass`](@ref), [`Highpass`](@ref), [`Bandpass`](@ref)
or [`Bandstop`](@ref) and includes the edges of the bands.

The design method is [`Butterworth`](@ref), [`Chebyshev1`](@ref), [`Chebyshev2`](@ref),
[`Elliptic`](@ref), or [`FIRWindow`](@ref), and includes any
necessary parameters for the method that affect the shape of the response,
such as filter order, ripple, and attenuation.

[Filter order estimation methods](@ref order-est-methods)
are available in [`buttord`](@ref), [`cheb1ord`](@ref), [`cheb2ord`](@ref),
and [`ellipord`](@ref) if the corner frequencies for different IIR filter types are known.
[`remezord`](@ref) can be used for an initial FIR filter order estimate.

```@docs
analogfilter
digitalfilter
```

For some filters, the design method is more general or
inherently implies a response type;
these [direct design methods](@ref direct-design-methods)
include [`remez`](@ref) which designs equiripple FIR
filters of all types, and [`iirnotch`](@ref) which designs a
2nd order "biquad" IIR notch filter.

For a more general application of creating a digital filter from s-domain
representation of an analog filter, one can use [`bilinear`](@ref):

```@docs
bilinear
```

### [Filter response types](@id response-types)

```@docs
Lowpass
Highpass
Bandpass
ComplexBandpass
Bandstop
```

The interpretation of the frequencies `Wn`, `Wn1` and `Wn2` depends on wether an analog
or a digital filter is designed.
1. If an analog filter is designed using [`analogfilter`](@ref), the frequencies are
   interpreted as analog frequencies in radians/second.
1. If a digital filter is designed using [`digitalfilter`](@ref) and the sampling
   frequency `fs` is specified, the frequencies of the filter response type are
   normalized to `fs`. This requires that the sampling frequency and the filter response
   type use the same frequency unit (Hz, radians/second, ...). If `fs` is not specified,
   the frequencies of the filter response type are interpreted as normalized frequencies
   in half-cycles/sample.


### [Filter design methods](@id design-methods)

#### IIR filter design methods

```@docs
Butterworth
Chebyshev1
Chebyshev2
Elliptic
```

### [Filter order estimation methods](@id order-est-methods)

#### IIR filter order estimation methods

```@docs
buttord
cheb1ord
cheb2ord
ellipord
```

#### FIR filter order estimation methods

```@docs
remezord
```

#### FIR filter design methods

```@docs
FIRWindow
```

### [Direct filter design methods](@id direct-design-methods)

```@docs
remez
iirnotch
```

## Filter response

```@docs
freqresp
phaseresp
grpdelay
impresp
stepresp
```

## Miscellaneous

```@docs
coefb
coefa
```

## Examples

Construct a 4th order elliptic lowpass filter with normalized cutoff
frequency 0.2, 0.5 dB of passband ripple, and 30 dB attentuation in
the stopband and extract the coefficients of the numerator and
denominator of the transfer function:

```julia
responsetype = Lowpass(0.2)
designmethod = Elliptic(4, 0.5, 30)
tf = convert(PolynomialRatio, digitalfilter(responsetype, designmethod))
numerator_coefs = coefb(tf)
denominator_coefs = coefa(tf)
```

Filter the data in `x`, sampled at 1000 Hz, with a 4th order
Butterworth bandpass filter between 10 and 40 Hz:

```julia
responsetype = Bandpass(10, 40)
designmethod = Butterworth(4)
filt(digitalfilter(responsetype, designmethod; fs=1000), x)
```

Filter the data in `x`, sampled at 50 Hz, with a 64 tap Hanning
window FIR lowpass filter at 5 Hz:

```julia
responsetype = Lowpass(5)
designmethod = FIRWindow(hanning(64))
filt(digitalfilter(responsetype, designmethod; fs=50), x)
```

Estimate a Lowpass Elliptic filter order with a normalized
passband cutoff frequency of 0.2, a stopband cutoff frequency of 0.4,
3 dB of passband ripple, and 40 dB attenuation in the stopband:

```julia
(N, ωn) = ellipord(0.2, 0.4, 3, 40)
```
