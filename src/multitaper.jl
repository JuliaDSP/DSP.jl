#####
##### Multitapered periodogram
#####

struct MTConfig{T,R1,F,P,T1,T2,W,R2}
    n_samples::Int
    fs::R1
    nfft::Int
    ntapers::Int
    freq::F
    plan::P
    fft_input_tmp::T1
    fft_output_tmp::T2
    window::W
    onesided::Bool
    r::R2 # inverse normalization; e.g. equal to `fs*N` for an unwindowed/untapered periodogram of a signal of length `N`
          # e.g. equal to `fs*ntapers*ones(ntapers)` when tapered by properly normalized window functions (i.e. norm-2 equal to 1 for each taper)
          # can be adjusted to weight the tapers, e.g. by the eigenvalues of the DPSS windows
    function MTConfig{T}(n_samples, nfft, ntapers, freq, fs, plan, fft_input_tmp,
                         fft_output_tmp, window, onesided, r) where {T}
        n_samples > 0 || throw(ArgumentError("`n_samples` must be positive"))
        nfft >= n_samples || throw(ArgumentError("Must have `nfft >= n_samples`"))
        ntapers > 0 || throw(ArgumentError("`ntapers` must be positive"))
        fs > 0 || throw(ArgumentError("`fs` must be positive"))
        if size(plan) != size(fft_input_tmp)
            throw(DimensionMismatch("""Must have `size(plan) == size(fft_input_tmp)`;
                got `size(plan)` = $(size(plan)) and `size(fft_input_tmp)` = $(size(fft_input_tmp))"""))
        end
        if size(fft_input_tmp) != (nfft,)
            throw(DimensionMismatch("""Must have `size(fft_input_tmp) == (nfft,)`;
                got `size(fft_input_tmp)` = $(size(fft_input_tmp)) and `nfft` = $(nfft)"""))
        end
        if size(fft_output_tmp) != (length(freq),)
            throw(DimensionMismatch("""Must have `size(fft_output_tmp) == (length(freq),)`;
                got `size(fft_output_tmp)` = $(size(fft_output_tmp)) and `length(freq)` = $(length(freq))"""))
        end
        if size(window) != (n_samples, ntapers)
            throw(DimensionMismatch("""Must have `size(window) == (n_samples, ntapers)`;
                got `size(window)` = $(size(window)) and `(ntapers, n_samples)` = $((n_samples, ntapers))"""))
        end
        if size(r) != (ntapers,)
            throw(DimensionMismatch("""Must have `size(r) == (ntapers,)`;
            got `size(r)` = $(size(r)) and `(ntapers,)` = $((ntapers,))"""))
        end
        if onesided && T <: Complex
            throw(ArgumentError("cannot compute one-sided FFT of a complex signal"))
        end
        return new{T,typeof(fs),typeof(freq),typeof(plan),typeof(fft_input_tmp),
                   typeof(fft_output_tmp),typeof(window),typeof(r)}(n_samples, fs, nfft,
                                                                    ntapers, freq, plan,
                                                                    fft_input_tmp,
                                                                    fft_output_tmp, window,
                                                                    onesided, r)
    end
end

# provides an `MTConfig` with `keep_only_large_evals` and `weight_by_evals` options used in some of the reference tests.
function dpss_config(::Type{T}, n_samples; nw=4, ntapers = 2nw-1, fs=1, keep_only_large_evals=false,
                     weight_by_evals=false, kwargs...) where {T}

    window = dpss(n_samples, nw, ntapers)

    if keep_only_large_evals
        evals = dpsseig(window, nw)
        keep_evals_mask = evals .> 0.9
        window = window[:, keep_evals_mask]
        evals = evals[keep_evals_mask]
        ntapers = size(window, 2)
    else
        evals = nothing
    end

    if weight_by_evals
        if evals === nothing
            evals = dpsseig(window, nw)
        end
        taper_weights = evals ./ sum(evals)
    else
        taper_weights = fill(1/ntapers, ntapers)
    end

    return MTConfig{T}(n_samples; window=window, nw=nw, ntapers=ntapers, taper_weights=taper_weights, fs=fs, kwargs...)
end


"""
    MTConfig{T}(n_samples; fs=1,
            nfft = nextpow(2, n_samples),
            window = nothing,
            nw = 4,
            ntapers = 2 * nw - 1,
            taper_weights = fill(1/ntapers, ntapers),
            onesided::Bool=T<:Real,
            fft_flags = FFTW.MEASURE)

Creates a config object which holds the configuration state
and temporary variables used in multitaper computations,
e.g. [`mt_pgram!`](@ref), [`mt_spectrogram`](@ref), [`MTSpectrogramConfig`](@ref),
[`MTCrossSpectraConfig`](@ref), and [`MTCoherenceConfig`](@ref).

An `MTConfig` can be re-used between computations as long
as none of the input arguments change.

* `n_samples`: the number of samples to be used as input when computing multitaper periodograms
  with this configuration. Used for pre-allocating temporary buffers.
* `fs`: the number of samples per second of the input signal
* `nfft`: length of input vector to the FFT; if `nfft > n_samples`, then the
  input signal will be zero-padded until it is of length `nfft`.
* `window`: window function to use for tapering. If left at the default of `nothing`,
  `window` will be set to `dpss(n_samples, nw, ntapers)`.
* `ntapers`: the number of tapers to use.
* `taper_weights = fill(1/ntapers, ntapers)`: how to weight the contribution of each taper.
  The default setting is to simply average them.
* `onesided`: whether or not to compute a "one-sided" FFT by using that real signal data
  yields conjugate-symmetry in Fourier space.
* `fft_flags`: flags to control how the FFT plan is generated.
"""
function MTConfig{T}(n_samples; fs=1, nfft=nextpow(2, n_samples), window=nothing, nw=4,
                     ntapers=2 * nw - 1, taper_weights = fill(1/ntapers, ntapers),
                    onesided::Bool=T <: Real,
                     fft_flags=FFTW.MEASURE) where {T}
    if onesided && T <: Complex
        throw(ArgumentError("cannot compute one-sided FFT of a complex signal"))
    end
    n_samples > 0 || throw(ArgumentError("`n_samples` must be positive"))
    nfft >= n_samples || throw(ArgumentError("Must have `nfft >= n_samples`"))
    freq = onesided ? rfftfreq(nfft, fs) : fftfreq(nfft, fs)
    fft_input_tmp = Vector{T}(undef, nfft)
    fft_output_tmp = Vector{fftouttype(T)}(undef, length(freq))
    plan = onesided ? plan_rfft(fft_input_tmp; flags=fft_flags) :
           plan_fft(fft_input_tmp; flags=fft_flags)
    if window === nothing
        r = fs ./ taper_weights
        window = dpss(n_samples, nw, ntapers)
    else
        r = fs .* vec(sum(abs2, window; dims=1)) ./ taper_weights
    end

    return MTConfig{T}(n_samples, nfft, ntapers, freq, fs, plan, fft_input_tmp,
                       fft_output_tmp, window, onesided, r)
end

function allocate_output(config::MTConfig{T}) where {T}
    return Vector{fftabs2type(T)}(undef, length(config.freq))
end


# Internal function used in [`mt_pgram!`](@ref) and [`mt_cross_power_spectra!`](@ref).
function mt_fft_tapered!(fft_output, signal, taper_index, config)
    # Create the input: tapered + zero-padded version of the signal
    fft_input = config.fft_input_tmp
    @inbounds for i in eachindex(signal)
        fft_input[i] = config.window[i, taper_index] * signal[i]
    end
    fft_input[(length(signal) + 1):(config.nfft)] .= 0

    # do the FFT
    mul!(fft_output, config.plan, fft_input)
    return nothing
end


"""
    mt_pgram(s; onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, nw=4, ntapers=iceil(2nw)-1, window=dpss(length(s), nw, ntapers))
    mt_pgram(signal::AbstractVector, config::MTConfig)

Computes the multitaper periodogram of a signal `s`.

If `window` is not specified, the signal is tapered with
`ntapers` discrete prolate spheroidal sequences with
time-bandwidth product `nw`. Each sequence is equally weighted;
adaptive multitaper is not (yet) supported.

If `window` is specified, each column is applied as a taper. The
sum of periodograms is normalized by the total sum of squares of
`window`.

Returns a `Periodogram`.

See also [`mt_pgram!`](@ref) and [`MTConfig`](@ref).
"""
mt_pgram

function mt_pgram(s::AbstractVector{T}; onesided::Bool=eltype(s)<:Real,
                  nfft::Int=nextfastfft(length(s)), fs::Real=1,
                  nw::Real=4, ntapers::Int=ceil(Int, 2nw)-1,
                  window::Union{AbstractMatrix,Nothing}=nothing) where T<:Number
    config = MTConfig{T}(length(s); fs = fs,
        nfft = nfft, window = window, nw = nw, ntapers = ntapers, onesided = onesided,
        fft_flags = FFTW.ESTIMATE)
    out = allocate_output(config)
    return mt_pgram!(out, s, config)
end

function mt_pgram!(output, s::AbstractVector{T}; onesided::Bool=eltype(s)<:Real,
    nfft::Int=nextfastfft(length(s)), fs::Real=1,
    nw::Real=4, ntapers::Int=ceil(Int, 2nw)-1,
    window::Union{AbstractMatrix,Nothing}=nothing) where T<:Number
    config = MTConfig{T}(length(s); fs = fs,
    nfft = nfft, window = window, nw = nw, ntapers = ntapers, onesided = onesided,
    fft_flags = FFTW.ESTIMATE)
    return mt_pgram!(output, s, config)
end

function mt_pgram(signal::AbstractVector, config::MTConfig)
    out = allocate_output(config)
    return mt_pgram!(out, signal, config)
end

"""
    mt_pgram!(output, s::AbstractVector{T}; onesided::Bool=eltype(s)<:Real,
        nfft::Int=nextfastfft(length(s)), fs::Real=1,
        nw::Real=4, ntapers::Int=ceil(Int, 2nw)-1,
        window::Union{AbstractMatrix,Nothing}=nothing) where T<:Number
    mt_pgram!(output::AbstractVector, signal::AbstractVector, config::MTConfig) -> Periodogram

Computes a multitapered periodogram, storing the output in `output`. Arguments:

* `signal::AbstractVector`: should be of length `config.n_samples`
* `output::AbstractVector`: should be of length `length(config.freq)`

Optionally pass an [`MTConfig`](@ref) object to preallocate temporary variables
and choose configuration settings; otherwise, keyword arguments may be passed
to choose those settings.

Returns a `Periodogram`.

See also [`mt_pgram`](@ref) and [`MTConfig`](@ref).
"""
mt_pgram!

function mt_pgram!(output::AbstractVector, signal::AbstractVector, config::MTConfig)
    if length(output) != length(config.freq)
        throw(DimensionMismatch("""Expected `output` to be of length `length(config.freq)`;
                got `length(output) = $(length(output)) and `length(config.freq)` = $(length(config.freq))"""))
    end
    if length(signal) != config.n_samples
        throw(DimensionMismatch("""Expected `signal` to be of length `config.n_samples`;
                got `length(signal) == $(length(signal)) and `config.n_samples` = $(config.n_samples)"""))
    end
    fft_output = config.fft_output_tmp

    output .= 0
    for taper_index in 1:(config.ntapers)
        mt_fft_tapered!(fft_output, signal, taper_index, config)
        fft2pow!(output, fft_output, config.nfft, config.r[taper_index], config.onesided)
    end
    return Periodogram(output, config.freq)
end

#####
##### Multitapered spectrogram
#####

struct MTSpectrogramConfig{T,C<:MTConfig{T}}
    n_samples::Int
    n_overlap_samples::Int
    time::Float64Range
    mt_config::C
end

"""
    MTSpectrogramConfig(n_samples, mt_config::MTConfig{T}, n_overlap_samples) where {T}
    MTSpectrogramConfig{T}(n_samples, samples_per_window, n_overlap_samples; fs=1, kwargs...) where {T}

Creates a `MTSpectrogramConfig` which holds configuration and temporary variables for [`mt_spectrogram`](@ref) and [`mt_spectrogram!`](@ref).
Any keyword arguments accepted by [`MTConfig`](@ref) may be passed here, or an `MTConfig` object itself.
"""
function MTSpectrogramConfig(n_samples::Int, mt_config::MTConfig{T}, n_overlap_samples::Int) where {T}
    samples_per_window = mt_config.n_samples
    if samples_per_window <= n_overlap_samples
        throw(ArgumentError("Need `samples_per_window > n_overlap_samples`; got `samples_per_window` = $(samples_per_window) and `n_overlap_samples` = $(n_overlap_samples)."))
    end

    f = samples_per_window / 2
    hop = samples_per_window - n_overlap_samples

    len = n_samples < samples_per_window ? 0 : div(n_samples - samples_per_window, samples_per_window - n_overlap_samples) + 1
    time = range(f, step = hop, length = len) / mt_config.fs
    return MTSpectrogramConfig{T,typeof(mt_config)}(n_samples, n_overlap_samples, time,
                                                    mt_config)
end

# add a method for if the user specifies the type, i.e. `MTSpectrogramConfig{T}`
function MTSpectrogramConfig{T}(n_samples::Int, mt_config::MTConfig{T}, n_overlap_samples::Int) where {T}
    return MTSpectrogramConfig(n_samples, mt_config, n_overlap_samples)
end

# Create the `MTConfig` if it's not passed
function MTSpectrogramConfig{T}(n_samples::Int, samples_per_window::Int,
                                n_overlap_samples::Int; fs=1, kwargs...) where {T}
    return MTSpectrogramConfig(n_samples, MTConfig{T}(samples_per_window; fs=fs, kwargs...), n_overlap_samples)
end

"""
    mt_spectrogram!(output, signal::AbstractVector{T}, n::Int=length(signal) >> 3,
        n_overlap::Int=n >> 1; fs=1, onesided::Bool=T <: Real, kwargs...) where {T}
    mt_spectrogram!(destination::AbstractMatrix, signal::AbstractVector, config::MTSpectrogramConfig)

Computes a multitaper spectrogram using the parameters specified in `config`. Arguments:

* `destination`: `length(config.mt_config.freq)` x `length(config.time)` matrix. This can be created by `DSP.allocate_output(config)`.
* `signal`: vector of length `config.n_samples`
* `config`: optionally, pass an [`MTSpectrogramConfig`](@ref) object to hold temporary variables and configuration settings. Otherwise, settings arguments may be passed directly.

Returns a `Spectrogram`.

See also [`mt_spectrogram`](@ref).
"""
mt_spectrogram!

function mt_spectrogram!(output, signal::AbstractVector{T}, n::Int=length(signal) >> 3,
    n_overlap::Int=n >> 1; fs=1, onesided::Bool=T <: Real,
    kwargs...) where {T}
    config = MTSpectrogramConfig{T}(length(signal), n, n_overlap; fs=fs, onesided=onesided,
        fft_flags=FFTW.ESTIMATE, kwargs...)
    return mt_spectrogram!(output, signal, config)
end

@views function mt_spectrogram!(destination::AbstractMatrix, signal::AbstractVector,
                                config::MTSpectrogramConfig)
    if size(destination) != (length(config.mt_config.freq), length(config.time))
        throw(DimensionMismatch("""Expected `destination` to be of size `(length(config.mt_config.freq), length(config.time))`;
                got `size(destination) == $(size(destination)) and
                `(length(config.mt_config.freq), length(config.time))` = $((length(config.mt_config.freq), length(config.time)))"""))
    end
    if length(signal) != config.n_samples
        throw(DimensionMismatch("""Expected `signal` to be of length `config.n_samples`;
                got `length(signal) = $(length(signal)) and `config.n_samples` = $(config.n_samples)"""))
    end
    samples_per_window = config.mt_config.n_samples
    subepochs = arraysplit(signal, samples_per_window, config.n_overlap_samples)
    @assert length(subepochs) == length(config.time)
    for (time_index, subepoch) in enumerate(subepochs)
        mt_pgram!(destination[:, time_index], subepoch, config.mt_config)
    end
    return Spectrogram(destination, config.mt_config.freq, config.time)
end

function allocate_output(config::MTSpectrogramConfig{T}) where {T}
    return Matrix{fftabs2type(T)}(undef, length(config.mt_config.freq), length(config.time))
end

function mt_spectrogram(signal::AbstractVector, config::MTSpectrogramConfig)
    output = allocate_output(config)
    return mt_spectrogram!(output, signal, config)
end

"""
    mt_spectrogram(signal::AbstractVector{T}, n::Int=length(s) >> 3,
                                  n_overlap::Int=n >> 1; fs=1,
                                  onesided::Bool=T <: Real, kwargs...) where {T}
    mt_spectrogram(signal::AbstractVector, config::MTSpectrogramConfig)

Compute a multitaper spectrogram, returning a `Spectrogram` object.
Optionally pass a [`MTSpectrogramConfig`](@ref) object; otherwise, any additional keyword
arguments accepted by [`MTConfig`](@ref) may be passed to configure the tapering.

Returns a `Spectrogram`.

See also [`mt_spectrogram!`](@ref).
"""
mt_spectrogram

function mt_spectrogram(signal::AbstractVector{T}, n::Int=length(signal) >> 3,
                        n_overlap::Int=n >> 1; fs=1, onesided::Bool=T <: Real,
                        kwargs...) where {T}
    config = MTSpectrogramConfig{T}(length(signal), n, n_overlap; fs=fs, onesided=onesided,
                                    fft_flags=FFTW.ESTIMATE, kwargs...)
    X = allocate_output(config)
    return mt_spectrogram!(X, signal, config)
end




"""
    mt_spectrogram(signal::AbstractVector, mt_config::MTConfig,
                   n_overlap::Int=mt_config.n_samples >> 1) -> Spectrogram

Compute a multitaper spectrogram using an an [`MTConfig`](@ref) object to choose
the window size. Note that all the workspace variables are contained in the [`MTConfig`](@ref) object, and this object can be re-used between spectrogram computations, so that only the output needs to be allocated.

Returns a `Spectrogram`.

See also [`mt_spectrogram!`](@ref).

## Example

```julia
signal1 = sin.(10_000) .+ randn(10_000)
mt_config = MTConfig{Float64}(1000) # 1000 samples per window
spec1 = mt_spectrogram(signal1, mt_config, 500)
signal2 = cos.(2_000) .+ randn(2_000)
spec2 = mt_spectrogram(signal2, mt_config, 250) # same `mt_config`, different signal, different overlap
```
"""
function mt_spectrogram(signal::AbstractVector, mt_config::MTConfig{T},
                        n_overlap::Int=mt_config.n_samples >> 1) where {T}
    config = MTSpectrogramConfig{T}(length(signal), mt_config, n_overlap)
    X = allocate_output(config)
    return mt_spectrogram!(X, signal, config)
end

#####
##### Multitapered cross spectral matrix
#####

# the following code was initially ported from MNE-python before being partially rewritten.
# MNE-python's 3-clause BSD license can be found here: <https://github.com/mne-tools/mne-python/blob/c6b22e82b64da0ed0f52a5f28bd4711e47090f9e/LICENSE.txt>.


"""
    CrossPowerSpectra{T,F,A<:AbstractArray{T, 3}}

Access the power (an `n_channels` x `n_channels` x `length(freq)` array) via the function [`power`](@ref), and the frequencies by the function [`freq`](@ref).

See also [`mt_cross_power_spectra`](@ref) and [`mt_cross_power_spectra!`](@ref).
"""
struct CrossPowerSpectra{T,F,A<:AbstractArray{T, 3}}
    power::A
    freq::F
end

freq(c::CrossPowerSpectra) = c.freq
power(c::CrossPowerSpectra) = c.power

function check_onesided_real(mt_config::MTConfig{T}) where T
    if !((T <: Real) && mt_config.onesided)
        throw(ArgumentError("Only real data is supported (with the default choice of `onesided=true`) for this operation."))
    end
    return nothing
end

struct MTCrossSpectraConfig{T,T1,T2,T3,T4,F,T5,T6,C<:MTConfig{T}}
    n_channels::Int
    normalization_weights::T1
    x_mt::T2
    demean::Bool
    mean_per_channel::T3
    demeaned_signal::T4
    freq::F
    freq_range::T5
    freq_inds::T6
    ensure_aligned::Bool
    mt_config::C
    function MTCrossSpectraConfig{T,T1,T2,T3,T4,F,T5,T6,C}(n_channels, normalization_weights,
            x_mt, demean, mean_per_channel, demeaned_signal, freq, freq_range,
            freq_inds, ensure_aligned, mt_config) where {T,T1,T2,T3,T4,F,T5,T6,C}
        check_onesided_real(mt_config) # this restriction is artifical; the code needs to be generalized
        return new{T,T1,T2,T3,T4,F,T5,T6,C}(n_channels, normalization_weights, x_mt,
            demean, mean_per_channel, demeaned_signal, freq, freq_range,
            freq_inds, ensure_aligned, mt_config)
    end
end

"""
    MTCrossSpectraConfig{T}(n_channels, n_samples; fs=1, demean=false, freq_range=nothing,
                            ensure_aligned = T == Float32 || T == Complex{Float32}, kwargs...) where {T}
    MTCrossSpectraConfig(n_channels, mt_config::MTConfig{T}; demean=false, freq_range=nothing,
                         ensure_aligned = T == Float32 || T == Complex{Float32})

Creates a configuration object used for [`mt_cross_power_spectra`](@ref) and [`mt_cross_power_spectra!`](@ref).

* `n_channels`: the number of channels of the input.
* `n_samples`: the number of samples for each channel of the input.
* `demean`: if `true`, the channelwise mean will be subtracted from the input signals before the cross spectral powers are computed.
* `freq_range`: if `nothing`, all frequencies are retained. Otherwise, only frequencies between `first(freq_range)` and `last(freq_range)` are retained.
* `ensure_aligned = T == Float32 || T == Complex{Float32}`: perform an extra copy to ensure that the FFT output is memory-aligned
* Additionally, either pass an [`MTConfig`](@ref) object, or keyword arguments such as `fs` accepted by [`MTConfig`](@ref).

Returns a `CrossPowerSpectra` object.
"""
function MTCrossSpectraConfig{T}(n_channels, n_samples; fs=1, demean=false,
                                 freq_range=nothing,
                                 ensure_aligned = T == Float32 || T == Complex{Float32},
                                 kwargs...) where {T}
    mt_config = MTConfig{T}(n_samples; fs=fs, kwargs...)
    return MTCrossSpectraConfig{T}(n_channels, mt_config; demean=demean, freq_range=freq_range, ensure_aligned=ensure_aligned)
end

# extra method to ensure it's ok to pass the redundant type parameter {T}
MTCrossSpectraConfig{T}(n_channels, mt_config::MTConfig{T}; demean=false,
        freq_range=nothing, ensure_aligned = T == Float32 || T == Complex{Float32}) where {T} = MTCrossSpectraConfig(n_channels, mt_config; demean=demean,
        freq_range=freq_range, ensure_aligned = ensure_aligned)

function MTCrossSpectraConfig(n_channels, mt_config::MTConfig{T}; demean=false,
        freq_range=nothing, ensure_aligned = T == Float32 || T == Complex{Float32}) where {T}

    n_samples = mt_config.n_samples
    if demean
        mean_per_channel = Vector{T}(undef, n_channels)
        demeaned_signal = Matrix{T}(undef, n_channels, n_samples)
    else
        mean_per_channel = nothing
        demeaned_signal = nothing
    end

    normalization_weights = 2 ./ mt_config.r

    x_mt = Array{fftouttype(T),3}(undef, length(mt_config.freq), mt_config.ntapers,
                                  n_channels)
    if freq_range !== nothing
        freq_mask = first(freq_range) .< mt_config.freq .< last(freq_range)
        freq_inds = findall(freq_mask)
        freq = mt_config.freq[freq_mask]
    else
        freq_inds = eachindex(mt_config.freq)
        freq = mt_config.freq
    end
    return MTCrossSpectraConfig{T,typeof(normalization_weights),typeof(x_mt),
                                typeof(mean_per_channel),typeof(demeaned_signal),
                                typeof(freq),typeof(freq_range),typeof(freq_inds),
                                typeof(mt_config)}(n_channels, normalization_weights, x_mt, demean,
                                                   mean_per_channel, demeaned_signal, freq,
                                                   freq_range, freq_inds, ensure_aligned, mt_config)
end

function allocate_output(config::MTCrossSpectraConfig{T}) where {T}
    return Array{fftouttype(T),3}(undef, config.n_channels, config.n_channels,
                                  length(config.freq))
end

"""
    mt_cross_power_spectra!(output, signal::AbstractMatrix; fs=1, kwargs...)
    mt_cross_power_spectra!(output, signal::AbstractMatrix, config::MTCrossSpectraConfig)

Computes multitapered cross power spectra between channels of a signal. Arguments:

* `output`: `n_channels` x `n_channels` x `length(config.freq)`. Can be created by `DSP.allocate_output(config)`.
* `signal`: `n_channels` x `n_samples`
* `config`: `MTCrossSpectraConfig{T}`: optionally pass a [`MTCrossSpectraConfig`](@ref) to
  preallocate temporary and choose configuration settings.
  Otherwise, one may pass any keyword arguments accepted by this object.

Produces a `CrossPowerSpectra` object holding the `n_channels` x `n_channels` x `n_frequencies`
output array and the corresponding frequencies (accessed by [`freq`](@ref)).

See also [`mt_cross_power_spectra`](@ref) and [`MTCrossSpectraConfig`](@ref).
"""
mt_cross_power_spectra!

function mt_cross_power_spectra!(output, signal::AbstractMatrix{T}; fs=1, kwargs...) where {T}
    n_channels, n_samples = size(signal)
    config = MTCrossSpectraConfig{T}(n_channels, n_samples; fs=fs, fft_flags=FFTW.ESTIMATE,
                                     kwargs...)
    return mt_cross_power_spectra!(output, signal, config)
end

@views function mt_cross_power_spectra!(output, signal::AbstractMatrix,
                                   config::MTCrossSpectraConfig)
    if size(signal) != (config.n_channels, config.mt_config.n_samples)
        throw(DimensionMismatch("Size of `signal` does not match `(config.n_channels, config.mt_config.n_samples)`;
        got `size(signal)`=$(size(signal)) but `(config.n_channels, config.mt_config.n_samples)`=$((config.n_channels, config.mt_config.n_samples))"))
    end
    if size(output) != (config.n_channels, config.n_channels, length(config.freq_inds))
        throw(DimensionMismatch("Size of `output` does not match `(config.n_channels, config.n_channels, length(config.freq_inds))`;
        got `size(output)`=$(size(output)) but `(config.n_channels, config.n_channels, length(config.freq_inds))`=$((config.n_channels, config.n_channels, length(config.freq_inds)))"))
    end

    if config.demean
        mean!(config.mean_per_channel, signal)
        config.demeaned_signal .= signal .- config.mean_per_channel
        signal = config.demeaned_signal
    end

    x_mt = config.x_mt

    if config.ensure_aligned
        mt_fft_tapered_multichannel_ensure_aligned!(x_mt, signal, config)
    else
        mt_fft_tapered_multichannel!(x_mt, signal, config)
    end
    x_mt[1, :, :] ./= sqrt(2)
    if iseven(config.mt_config.nfft)
        x_mt[end, :, :] ./= sqrt(2)
    end
    cs_inner!(output, config.normalization_weights, x_mt, config)
    return CrossPowerSpectra(output, config.freq)
end

function mt_fft_tapered_multichannel_ensure_aligned!(x_mt, signal, config)
    fft_output = config.mt_config.fft_output_tmp
    for k in 1:(config.n_channels), taper in 1:config.mt_config.ntapers
        # we do this in two steps so that we are sure `fft_output` has the memory alignment FFTW expects (without needing the `FFTW.UNALIGNED` flag)
        mt_fft_tapered!(fft_output, signal[k, :], taper, config.mt_config)
        x_mt[:, taper, k] .= fft_output
    end
end

@views function mt_fft_tapered_multichannel!(x_mt, signal, config)
    for k in 1:(config.n_channels), taper in 1:config.mt_config.ntapers
        mt_fft_tapered!(x_mt[:, taper, k], signal[k, :], taper, config.mt_config)
    end
end

function cs_inner!(output, normalization_weights, x_mt, config)
    freq_inds = config.freq_inds
    n_channels = config.n_channels
    @boundscheck checkbounds(output, 1:n_channels, 1:n_channels, 1:length(freq_inds))
    @boundscheck checkbounds(normalization_weights, 1:length(normalization_weights))
    @boundscheck checkbounds(x_mt, freq_inds, 1:length(normalization_weights), 1:n_channels)
    output .= zero(eltype(output))
    # Up to the `normalization_weights` scaling, we have
    # J_k^l(f) = x_mt[k, f, l]
    # Ŝ^lm(f) = output[l, m, f]
    # using the notation from https://en.wikipedia.org/wiki/Multitaper#The_method
    # so the formula `Ŝ^lm(f) = \sum_k conj(J_k^l(f)) * (J_k^m(f))`` becomes the following loop:
    @inbounds for (fi, f) in enumerate(freq_inds),
                  m in 1:n_channels,
                  l in 1:n_channels,
                  k in eachindex(normalization_weights)

        output[l, m, fi] += normalization_weights[k] * x_mt[f, k, l] * conj(x_mt[f, k, m])
    end
    return nothing
end

"""
    mt_cross_power_spectra(signal::AbstractMatrix{T}; fs=1, kwargs...) where {T}
    mt_cross_power_spectra(signal::AbstractMatrix, config::MTCrossSpectraConfig)

Computes multitapered cross power spectra between channels of a signal. Arguments:

* `signal`: `n_channels` x `n_samples`
* Optionally pass an [`MTCrossSpectraConfig`](@ref) object to preallocate temporary variables
and choose configuration settings. Otherwise, any keyword arguments accepted by [`MTCrossSpectraConfig`](@ref) may be passed here.

Produces a `CrossPowerSpectra` object holding the `n_channels` x `n_channels` x `n_frequencies`
output array (accessed by [`power`](@ref)) and the corresponding frequencies (accessed by [`freq`](@ref)).

See also [`mt_cross_power_spectra!`](@ref) and [`MTCrossSpectraConfig`](@ref).
"""
mt_cross_power_spectra

function mt_cross_power_spectra(signal::AbstractMatrix{T}; fs=1, kwargs...) where {T}
    n_channels, n_samples = size(signal)
    config = MTCrossSpectraConfig{T}(n_channels, n_samples; fs=fs, fft_flags=FFTW.ESTIMATE,
                                     kwargs...)
    return mt_cross_power_spectra(signal, config)
end

function mt_cross_power_spectra(signal::AbstractMatrix, config::MTCrossSpectraConfig)
    output = allocate_output(config)
    return mt_cross_power_spectra!(output, signal, config)
end

#####
##### Multitapered coherences
#####

struct MTCoherenceConfig{T,T1,C<:MTCrossSpectraConfig{T}}
    cs_matrix::T1
    cs_config::C
end

"""
    MTCoherenceConfig{T}(n_channels, n_samples; fs=1, demean=false, freq_range=nothing, kwargs...) where T
    MTCoherenceConfig(cs_config::MTCrossSpectraConfig{T}) where {T}
    MTCoherenceConfig(n_channels, mt_config::MTConfig{T}; demean=false, freq_range=nothing,
        ensure_aligned = T == Float32 || T == Complex{Float32}) where {T}

Creates a configuration object for coherences from a [`MTCrossSpectraConfig`](@ref). Provides a helper method
with the same arugments as `MTCrossSpectraConfig` to construct the `MTCrossSpectraConfig` object.

See also [`mt_coherence`](@ref) and [`mt_coherence!`](@ref).
"""
function MTCoherenceConfig(cs_config::MTCrossSpectraConfig{T}) where {T}
    cs_matrix = allocate_output(cs_config)
    return MTCoherenceConfig{T,typeof(cs_matrix),typeof(cs_config)}(cs_matrix, cs_config)
end

# add a method to cover the case in which the user specifies the `{T}` here
MTCoherenceConfig{T}(cs_config::MTCrossSpectraConfig{T}) where {T} = MTCoherenceConfig(cs_config)

function MTCoherenceConfig{T}(n_channels, n_samples; fs=1, demean=false,
                              freq_range=nothing, kwargs...) where {T}
    cs_config = MTCrossSpectraConfig{T}(n_channels, n_samples; fs=fs, demean=demean,
                                        freq_range=freq_range, kwargs...)
    cs_matrix = allocate_output(cs_config)
    return MTCoherenceConfig{T,typeof(cs_matrix),typeof(cs_config)}(cs_matrix, cs_config)
end

# ensure it's OK to pass the extra {T} type parameter
function MTCoherenceConfig{T}(n_channels, mt_config::MTConfig{T}; demean=false,
    freq_range=nothing, ensure_aligned = T == Float32 || T == Complex{Float32}) where T
    return MTCoherenceConfig(n_channels, mt_config; demean=demean,
    freq_range=freq_range, ensure_aligned = ensure_aligned)
end

function MTCoherenceConfig(n_channels, mt_config::MTConfig{T}; demean=false,
    freq_range=nothing, ensure_aligned = T == Float32 || T == Complex{Float32}) where {T}
    cs_config = MTCrossSpectraConfig(n_channels, mt_config; demean=demean,
              freq_range=freq_range, ensure_aligned=ensure_aligned)
    cs_matrix = allocate_output(cs_config)
    return MTCoherenceConfig{T,typeof(cs_matrix),typeof(cs_config)}(cs_matrix, cs_config)
end


function allocate_output(config::MTCoherenceConfig{T}) where {T}
    return Array{real(T)}(undef, config.cs_config.n_channels, config.cs_config.n_channels, length(config.cs_config.freq))
end

"""
    coherence_from_cs!(output::AbstractArray, cs_matrix)

Compute the pairwise channel coherences from a cross spectral matrix, storing the result in `output`.
"""
function coherence_from_cs!(output::AbstractArray{T}, cs_matrix) where T
    n_channels = size(cs_matrix, 2)
    n_freqs = size(cs_matrix, 3)
    @boundscheck checkbounds(output, 1:n_channels, 1:n_channels, 1:n_freqs)
    @boundscheck checkbounds(cs_matrix, 1:n_channels, 1:n_channels, 1:n_freqs)
    output .= zero(T)
    @inbounds for ch2 in 1:n_channels
        for ch1 in (ch2 + 1):n_channels # lower triangular matrix
            for f in 1:n_freqs
                output[ch1, ch2, f] += abs(cs_matrix[ch1, ch2, f]) /
                                         sqrt(real(cs_matrix[ch1, ch1, f] *
                                               cs_matrix[ch2, ch2, f]))
            end
        end
    end
    output .+= PermutedDimsArray(output, (2, 1, 3)) # symmetrize
    # diagonal elements should be `1` for any frequency
    @inbounds for i in 1:n_channels
        output[i, i, :] .= one(T)
    end
    return nothing
end

"""
    Coherence

Holds an `n_channels` x `n_channels` x `length(freq)` array consisting of the pairwise coherences between each channel for each frequency which is accessed by [`coherence`](@ref), as well
as the frequencies accessed by [`freq`](@ref).

See also [`mt_coherence`](@ref) and [`mt_coherence!`](@ref).
"""
struct Coherence{T,F,A<:AbstractArray{T, 3}}
    coherence::A
    freq::F
end

freq(c::Coherence) = c.freq

"""
    coherence(c::Coherence)

Given an `Coherence` object, returns an `n_channels` x `n_channels` x `length(freq(c))`
array consisting of the pairwise coherences between each channel for each frequency.
"""
coherence(c::Coherence) = c.coherence

"""
    mt_coherence!(output, signal::AbstractMatrix; fs=1, freq_range=nothing, demean=false, kwargs...)
    mt_coherence!(output, signal::AbstractMatrix, config::MTCoherenceConfig)

Computes the pairwise coherences between channels.

* `output`: `n_channels` x `n_channels` matrix
* `signal`: `n_samples` x `n_channels` matrix
* `config`: optional configuration object that pre-allocates temporary variables and choose settings.

Returns a `Coherence` object.

See also [`mt_coherence`](@ref) and [`MTCoherenceConfig`](@ref).
"""
mt_coherence!

function mt_coherence!(output, signal::AbstractMatrix,
                       config::MTCoherenceConfig)
    if size(signal) != (config.cs_config.n_channels, config.cs_config.mt_config.n_samples)
        throw(DimensionMismatch("Size of `signal` does not match `(config.cs_config.n_channels, config.cs_config.mt_config.n_samples)`;
            got `size(signal)`=$(size(signal)) but `(config.cs_config.n_channels, config.cs_config.mt_config.n_samples)`=$((config.cs_config.n_channels, config.cs_config.mt_config.n_samples))"))
    end
    if size(output) != (config.cs_config.n_channels, config.cs_config.n_channels, length(config.cs_config.freq))
        throw(DimensionMismatch("Size of `output` does not match `(config.cs_config.n_channels, config.cs_config.n_channels, length(config.cs_config.freq))`;
        got `size(output)`=$(size(output)) but `(config.cs_config.n_channels, config.cs_config.n_channels, length(config.cs_config.freq))`=$((config.cs_config.n_channels, config.cs_config.n_channels, length(config.cs_config.freq)))"))
    end
    cs = mt_cross_power_spectra!(config.cs_matrix, signal, config.cs_config)
    coherence_from_cs!(output, power(cs))

    return Coherence(output, config.cs_config.freq)
end

function mt_coherence!(output, signal::AbstractMatrix{T}; fs=1, freq_range=nothing, demean=false,
    kwargs...) where {T}
    n_channels, n_samples = size(signal)
    config = MTCoherenceConfig{T}(n_channels, n_samples; fs=fs, demean=demean,
                freq_range=freq_range,
                fft_flags=FFTW.ESTIMATE, kwargs...)
    return mt_coherence!(output, signal, config)
end


"""
    mt_coherence(signal::AbstractMatrix{T}; fs=1, freq_range = nothing, demean=false, kwargs...) where T
    mt_coherence(signal::AbstractMatrix, config::MTCoherenceConfig)

Arguments:

* `signal`: `n_channels` x `n_samples` matrix
* Optionally pass an `MTCoherenceConfig` to pre-allocate temporary variables and choose configuration settings, otherwise, see [`MTCrossSpectraConfig`](@ref) for the meaning of the keyword arguments.

Returns a `Coherence` object.

See also [`mt_coherence`](@ref) and [`MTCoherenceConfig`](@ref).
"""
mt_coherence

function mt_coherence(signal::AbstractMatrix{T}; fs=1, freq_range=nothing, demean=false,
                      kwargs...) where {T}
    n_channels, n_samples = size(signal)
    config = MTCoherenceConfig{T}(n_channels, n_samples; fs=fs, demean=demean,
                                  freq_range=freq_range,
                                  fft_flags=FFTW.ESTIMATE, kwargs...)
    return mt_coherence(signal, config)
end


function mt_coherence(signal::AbstractMatrix, config::MTCoherenceConfig)
    cohs = allocate_output(config)
    return mt_coherence!(cohs, signal, config)
end
