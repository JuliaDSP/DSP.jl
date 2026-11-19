# `Convolutions` - similarity methods

```@docs
conv
conv!
deconv
xcorr
```

## Nonfinite inputs

Convolution algorithms can propagate `NaN` and `Inf` differently because
FFT-based methods reassociate operations. Select an algorithm explicitly
when the propagation of nonfinite values matters:

```jldoctest; setup = :(using DSP)
julia> x = [NaN, 1.0, 2.0, 3.0, 4.0, 5.0]; h = [1.0, 1.0, 1.0];

julia> println(conv(x, h; algorithm=:direct))
[NaN, NaN, NaN, 6.0, 9.0, 12.0, 9.0, 5.0]

julia> println(conv(x, h; algorithm=:fft))
[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN]
```
