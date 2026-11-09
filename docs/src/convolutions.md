# `Convolutions` - similarity methods

```@docs
conv
conv!
deconv
xcorr
```

## Nonfinite inputs

Convolution algorithms can propagate `NaN` and `Inf` differently. With
`algorithm=:direct`, a nonfinite value participates only in the convolution
sums that overlap its position. FFT-based algorithms can spread it further:
`:fft_simple` transforms the whole input, whereas `:fft_overlapsave` transforms
overlapping blocks, so contamination can extend across an affected block.
The affected positions can therefore depend on the algorithm and block sizes.

For example, explicitly selecting `:direct` keeps the finite tail in this
short signal, while `:fft_simple` produces only `NaN`s:

```jldoctest; setup = :(using DSP)
julia> x = [NaN, 1.0, 2.0, 3.0, 4.0, 5.0]; h = [1.0, 1.0, 1.0];

julia> direct = conv(x, h; algorithm=:direct);

julia> findall(isnan, direct) == [1, 2, 3]
true

julia> direct[4:end] == [6.0, 9.0, 12.0, 9.0, 5.0]
true

julia> all(isnan, conv(x, h; algorithm=:fft_simple))
true
```

An `Inf` input can likewise produce `NaN`s in FFT-based convolution rather
than the `Inf`s produced by direct summation. The default `algorithm=:auto`
uses performance heuristics based on types and sizes; it does not scan for
nonfinite values. If their propagation matters to an application, check the
inputs and choose an algorithm explicitly (see [`conv!`](@ref)).

These differences do not make `:direct` a guarantee of finite or more accurate
results. Even finite inputs can overflow to `Inf` or produce `NaN`, and
rounding and overflow can depend on arithmetic order. In particular, `NaN`
propagation should not be used as an algorithm-independent missing-data mask.
