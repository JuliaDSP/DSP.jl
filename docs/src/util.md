# `Util` - utility functions

!!! note
    As of version 0.6.1 of DSP.jl, `fftfreq` and `rfftfreq` have been moved from
    DSP.jl to [AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl)
    version 0.5 and above. You can also access these functions through
    [FFTW.jl](https://github.com/JuliaMath/FFTW.jl) version 1.1 and above.

```@docs
unwrap
unwrap!
hilbert
nextfastfft
pow2db
amp2db
db2pow
db2amp
rms
rmsfft
meanfreq
finddelay
shiftin!
shiftsignal
shiftsignal!
alignsignals
alignsignals!
diric
```
