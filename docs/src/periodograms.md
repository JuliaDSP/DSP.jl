# `Periodograms` - periodogram estimation

## Basic functions
Common procedures like computing the [short-time Fourier transform](@ref stft),
[`periodogram`](@ref)s, and [`spectrogram`](@ref)s are documented below.

```@docs
arraysplit
DSP.Periodograms.ArraySplit
periodogram(s::AbstractVector{T}) where T <: Number
periodogram(s::AbstractMatrix{T}) where T <: Real
DSP.Periodograms.Periodogram
DSP.Periodograms.Periodogram2
welch_pgram
welch_pgram!
spectrogram
DSP.Periodograms.Spectrogram
stft
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
WelchConfig
MTConfig
MTSpectrogramConfig
MTCrossSpectraConfig
MTCoherenceConfig
```
