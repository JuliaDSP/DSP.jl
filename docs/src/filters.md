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
resample
```
## Filter design

Most analog and digital filters are constructed by composing
[response types](@ref response-types), which determine the frequency
response of the filter, with [design methods](@ref design-methods),
which determine how the filter is constructed.

```@docs
analogfilter
digitalfilter
```

For some filters, the design method inherently implies a response type.
Such filters are documented below.

```@docs
iirnotch
```

### [Filter response types](@id response-types)

```@docs
Lowpass
Highpass
Bandpass
Bandstop
```

### [Filter design methods](@id design-methods)

#### IIR filter design methods

```@docs
Butterworth
Chebyshev1
Chebyshev2
Elliptic
```

#### FIR filter design methods

```@docs
FIRWindow
remez
```

## Filter response

```@docs
freqz
phasez
impz
stepz
freqs
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
responsetype = Bandpass(10, 40; fs=1000)
designmethod = Butterworth(4)
filt(digitalfilter(responsetype, designmethod), x)
```

Filter the data in `x`, sampled at 50 Hz, with a 64 tap Hanning
window FIR lowpass filter at 5 Hz:

```julia
responsetype = Lowpass(5; fs=50)
designmethod = FIRWindow(hanning(64))
filt(digitalfilter(responsetype, designmethod), x)
```
