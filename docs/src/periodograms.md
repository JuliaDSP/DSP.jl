# `Periodograms` - periodogram estimation

## Basic functions
Common procedures like computing the [short-time Fourier transform](@ref stft),
[`periodogram`](@ref)s, and [`spectrogram`](@ref)s are documented below.

```@docs
arraysplit
periodogram(s::AbstractVector{T}) where T <: Number
welch_pgram
spectrogram
stft
periodogram(s::AbstractMatrix{T}) where T <: Real
freq
power
time
coherence
DSP.Periodograms.Coherence
```

## Multitaper periodogram estimation

```@docs
mt_pgram
mt_pgram!
mt_spectrogram
mt_spectrogram!
mt_cross_power_spectra
mt_cross_power_spectra!
DSP.Periodograms.CrossPowerSpectra
mt_coherence
mt_coherence!
```

## Configuration objects

```@docs
MTConfig
MTSpectrogramConfig
MTCrossSpectraConfig
MTCoherenceConfig
```
