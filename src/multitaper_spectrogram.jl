#####
##### Multitapered periodogram
#####

struct MTConfig{T,R1,F,P,T1,T2,W,R2}
    n_samples::Int
    fs::R1
    n_for_fft::Int
    n_freq_bins::Int
    n_tapers::Int
    freq::F
    plan::P
    fft_input_tmp::T1
    fft_output_tmp::T2
    window::W
    onesided::Bool
    r::R2
    function MTConfig{T}(n_samples, n_for_fft, n_freq_bins, n_tapers, freq, fs,
                         plan, fft_input_tmp, fft_output_tmp, window, onesided, r) where {T}
        n_samples > 0 || throw(ArgumentError("`n_samples` must be positive"))
        n_for_fft >= n_samples || throw(ArgumentError("Must have `n_for_fft >= n_samples`"))
        n_tapers > 0 || throw(ArgumentError("`n_tapers` must be positive"))
        fs > 0 || throw(ArgumentError("`fs` must be positive"))
        if size(plan) != size(fft_input_tmp)
            throw(DimensionMismatch("""Must have `size(plan) == size(fft_input_tmp)`;
                got `size(plan)` = $(size(plan)) and `size(fft_input_tmp)` = $(size(fft_input_tmp))"""))
        end
        if size(fft_input_tmp) != (n_for_fft,)
            throw(DimensionMismatch("""Must have `size(fft_input_tmp) == (n_for_fft,)`;
                got `size(fft_input_tmp)` = $(size(fft_input_tmp)) and `n_for_fft` = $(n_for_fft)"""))
        end
        if size(fft_output_tmp) != (n_freq_bins, n_tapers)
            throw(DimensionMismatch("""Must have `size(fft_output_tmp) == (n_freq_bins, n_tapers)`;
                got `size(fft_output_tmp)` = $(size(fft_output_tmp)) and `(n_freq_bins, n_tapers)` = $((n_freq_bins, n_tapers))"""))
        end
        if size(window) != (n_samples, n_tapers)
            throw(DimensionMismatch("""Must have `size(window) == (n_samples, n_tapers)`;
                got `size(window)` = $(size(window)) and `(n_tapers, n_samples)` = $((n_samples, n_tapers))"""))
        end
        if onesided && T <: Complex
            throw(ArgumentError("cannot compute one-sided FFT of a complex signal"))
        end
        return new{T,typeof(fs),typeof(freq),typeof(plan),typeof(fft_input_tmp),
                   typeof(fft_output_tmp),typeof(window),typeof(r)}(n_samples, fs,
                                                                                                                                          n_for_fft,
                                                                                                                                          n_freq_bins,
                                                                                                                                          n_tapers,
                                                                                                                                          freq,
                                                                                                                                          plan,
                                                                                                                                          fft_input_tmp,
                                                                                                                                          fft_output_tmp,
                                                                                                                                          window,
                                                                                                                                          onesided,
                                                                                                                                          r)
    end
end

"""
    MTConfig{T}(;n_samples, fs
            n_for_fft = nextpow(2, n_samples),
            window = nothing,
            time_bandwidth_product = 4,
            n_tapers = 2 * time_bandwidth_product - 1,
            onesided::Bool=T<:Real)

Creates a config object which holds the configuration state
and temporary variables used in multitaper computations,
e.g. [`mt_pgram!`](@ref), [`mt_spectrogram`](@ref), and [`mt_cross_spectral!`](@ref).

An `MTConfig` can be re-used between computations as long
as none of the input arguments change.

* `n_samples`: the number of samples to be used as input when computing multitaper periodograms
  with this configuration. Used for pre-allocating temporary buffers.
* `fs`: the number of samples per second of the input signal
* `n_for_fft`: length of input vector to the FFT; if `n_for_fft > n_samples`, then the
  input signal will be zero-padded until it is of length `n_for_fft`.
* `window`: window function to use for tapering. If left at the default of `nothing`,
  `window` will be set to `dpss(n_samples, time_bandwidth_product, n_tapers)`.
* `n_tapers`: the number of tapers to use.
* `onesided`: whether or not to compute a "one-sided" FFT by using that real signal data
  yields conjugate-symmetry in Fourier space.
"""
function MTConfig{T}(n_samples, fs; n_for_fft=nextpow(2, n_samples),
                     window=nothing, time_bandwidth_product=4,
                     n_tapers=2 * time_bandwidth_product - 1,
                     onesided::Bool=T <: Real) where {T}
    if onesided && T <: Complex
        throw(ArgumentError("cannot compute one-sided FFT of a complex signal"))
    end
    n_samples > 0 || throw(ArgumentError("`n_samples` must be positive"))
    n_for_fft >= n_samples || throw(ArgumentError("Must have `n_for_fft >= n_samples`"))
    n_tapers > 0 || throw(ArgumentError("`n_tapers` must be positive"))
    freq = onesided ? rfftfreq(n_for_fft, fs) : fftfreq(n_for_fft, fs)
    n_freq_bins = length(freq)
    fft_input_tmp = Vector{T}(undef, n_for_fft)
    fft_output_tmp = Matrix{fftouttype(T)}(undef, n_freq_bins, n_tapers)
    plan = onesided ? plan_rfft(fft_input_tmp) : plan_fft(fft_input_tmp)
    if window === nothing
        r = fs * n_tapers
        window = dpss(n_samples, time_bandwidth_product, n_tapers)
    else
        r = fs * sum(abs2, window)
    end

    return MTConfig{T}(n_samples, n_for_fft, n_freq_bins, n_tapers, freq, fs, plan,
                       fft_input_tmp, fft_output_tmp, window, onesided, r)
end

allocate_output(config::MTConfig{T}) where T = Vector{fftabs2type(T)}(undef, config.n_freq_bins)

"""
    tapered_spectra!(output, signal::AbstractVector, config) where {T}

Internal function used in [`mt_pgram!`](@ref) and [`mt_spectrogram!`](@ref).

* `output`: `length(frequencies)` x `n_tapers`
* `signal`: `n_samples`
* `config.fft_input_tmp`: `n_for_fft`
* `config.window`: `n_samples` x `n_tapers`
* `config.plan`: 1D plan for `n_for_fft`

"""
@views function tapered_spectra!(output, signal::AbstractVector, config::MTConfig)
    if length(signal) != config.n_samples
        throw(DimensionMismatch("Expected `length(signal) == config.n_samples`; got `length(signal)` = $(length(signal)) and `config.n_samples` = $(config.n_samples)"))
    end
    if size(output) != (config.n_freq_bins, config.n_tapers)
        throw(DimensionMismatch("Expected `size(output) == (config.n_freq_bins, config.n_tapers)`; got `size(output)` = $(size(output)) and `(config.n_freq_bins, config.n_tapers)` = $((config.n_freq_bins, config.n_tapers))"))

    end
    input = config.fft_input_tmp
    @inbounds for j in 1:(config.n_tapers)
        for i in 1:length(signal)
            input[i] = config.window[i, j] * signal[i]
        end
        # zero-pad
        input[(length(signal) + 1):(config.n_for_fft)] .= 0

        mul!(output[:, j], config.plan, input)
    end
    return nothing
end

"""
    mt_pgram!(output::AbstractVector, signal::AbstractVector, config::MTConfig) -> Periodogram

Computes a multitapered periodogram with parameters specifed by `config`,
storing the output in `output`.

* `signal::AbstractVector`: should be of length `config.n_samples`
* `output::AbstractVector`: should be of length `config.n_freq_bins`
"""
function mt_pgram!(output::AbstractVector, signal::AbstractVector, config::MTConfig)
    if length(output) != config.n_freq_bins
        throw(DimensionMismatch("""Expected `output` to be of length `config.n_freq_bins`;
                got `length(output) = $(length(output)) and `config.n_freq_bins` = $(config.n_freq_bins)"""))
    end
    if length(signal) != config.n_samples
        throw(DimensionMismatch("""Expected `signal` to be of length `config.n_samples`;
                got `length(signal) == $(length(signal)) and `config.n_samples` = $(config.n_samples)"""))
    end
    tmp = config.fft_output_tmp
    tapered_spectra!(tmp, signal, config)

    output .= 0
    for j in 1:(config.n_tapers)
        fft2pow!(output, @view(tmp[:, j]), config.n_for_fft, config.r, config.onesided)
    end
    return Periodogram(output, config.freq)
end

#####
##### Multitapered spectrogram
#####

struct MTSpectrogramConfig{T,C<:MTConfig{T}}
    n_samples::Int
    n_overlap_samples::Int
    n_time_points::Int
    time::FloatRange{Float64}
    mt_config::C
end

"""
    MTSpectrogramConfig{T}(; n_samples, n_overlap_samples, mt_config::MTConfig{T}) where {T}
    MTSpectrogramConfig{T}(; n_samples, n_overlap_samples, samples_per_window, fs, kwargs...) where {T}

Creates a `MTSpectrogramConfig` which holds configuration and temporary variables for [`mt_spectrogram!`](@ref).
Any keyword arguments accepted by [`MTConfig`](@ref) may be passed here, or an `MTConfig` object itself.
"""
function MTSpectrogramConfig{T}(n_samples, n_overlap_samples;
                                          mt_config::MTConfig{T}) where {T}
    samples_per_window = mt_config.n_samples
    fs = mt_config.fs
    # Dimensions of the spectrogram power matrix
    n_time_points = div(n_samples - samples_per_window,
                        samples_per_window - n_overlap_samples) + 1

    window_duration = samples_per_window / fs
    overlap_duration = n_overlap_samples / fs

    hop = window_duration - overlap_duration
    time = range(window_duration / 2; step=hop, length=n_time_points)
    return MTSpectrogramConfig{T,typeof(mt_config)}(n_samples, n_overlap_samples,
                                                              n_time_points, time,
                                                              mt_config)
end

function MTSpectrogramConfig{T}(n_samples, n_overlap_samples, samples_per_window;
                                          fs, kwargs...) where {T}
    mt_config = MTConfig{T}(samples_per_window, fs; kwargs...)
    return MTSpectrogramConfig{T}(n_samples, n_overlap_samples;
                                            mt_config=mt_config)
end

"""
    mt_spectrogram!(destination::AbstractMatrix,
                                signal::AbstractVector,
                                config::MTSpectrogramConfig) -> Spectrogram

Computes a multitaper spectrogram using the parameters specified in `config`.

* `destination`: `config.mt_config.n_freq_bins` x `config.n_time_points` matrix. This can be created by `allocate_output(config)`.
* `signal`: vector of length `config.n_samples`
* `config`: an [`MTSpectrogramConfig`](@ref) object to hold temporary variables and configuration settings.
"""
@views function mt_spectrogram!(destination::AbstractMatrix,
                                          signal::AbstractVector,
                                          config::MTSpectrogramConfig)
    if size(destination) != (config.mt_config.n_freq_bins, config.n_time_points)
        throw(DimensionMismatch("""Expected `destination` to be of size `(config.mt_config.n_freq_bins, config.n_time_points)`;
                got `size(destination) == $(size(destination)) and
                `(config.mt_config.n_freq_bins, config.n_time_points)` = $((config.mt_config.n_freq_bins, config.n_time_points))"""))
    end
    if length(signal) != config.n_samples
        throw(DimensionMismatch("""Expected `signal` to be of length `config.n_samples`;
                got `length(signal) = $(length(signal)) and `config.n_samples` = $(config.n_samples)"""))
    end
    samples_per_window = config.mt_config.n_samples
    subepochs = arraysplit(signal, samples_per_window, config.n_overlap_samples)
    for (time_index, subepoch) in enumerate(subepochs)
        mt_pgram!(destination[:, time_index], subepoch, config.mt_config)
    end
    return Spectrogram(destination, config.mt_config.freq, config.time)
end

allocate_output(config::MTSpectrogramConfig{T}) where {T} = Matrix{fftabs2type(T)}(undef, config.mt_config.n_freq_bins, config.n_time_points)

"""
    mt_spectrogram(signal::AbstractVector{T}, n::Int=length(s) >> 3,
                                  n_overlap::Int=n >> 1; fs::Int=1,
                                  onesided::Bool=T <: Real, kwargs...) where {T} -> Spectrogram

Compute a multitaper spectrogram. Any additional keyword arguments accepted by [`MTConfig`](@ref)
may be passed to configure the tapering. See also [`mt_spectrogram!`](@ref).
"""
function mt_spectrogram(signal::AbstractVector{T}, n::Int=length(signal) >> 3,
                                  n_overlap::Int=n >> 1; fs::Int=1,
                                  onesided::Bool=T <: Real, kwargs...) where {T}
    config = MTSpectrogramConfig{T}(length(signal), n_overlap, n; fs=fs,
                                              onesided=onesided, kwargs...)
    X = allocate_output(config)
    return mt_spectrogram!(X, signal, config)
end

#####
##### Multitapered cross spectral matrix
#####

struct CrossSpectral{T, F}
    values::T
    freq::F
end

freq(c::CrossSpectral) = c.freq

struct MTCrossSpectraConfig{T,T1,T2,T3,T4,F,T5,T6,C<:MTConfig{T}}
    n_channels::Int
    weighted_evals::T1
    x_mt::T2
    demean::Bool
    channelwise_mean::T3
    demeaned_signal::T4
    freq::F
    freq_range::T5
    freq_inds::T6
    mt_config::C
end

"""
    MTCrossSpectraConfig{T}(n_samples, n_channels, fs; demean=true, low_bias=true, freq_range=nothing, time_bandwidth_product=4, kwargs...) where T

Creates a configuration object used for [`mt_cross_spectral!`](@ref) as well as [`mt_coherence!`](@ref).


* `n_samples`: the number of samples to be used as input when computing multitaper periodograms
    with this configuration. Used for pre-allocating temporary buffers.
* `fs`: the number of samples per second of the input signal
* `demean`: if `true`, the channelwise mean will be subtracted from the input signals before the cross spectral values are computed.
* `low_bias`: if `true`, keeps only tapers with eigenvalues > `0.9`.
* `freq_range`: if `nothing`, all frequencies are retained. Otherwise, only frequencies between `first(freq_range)` and `last(freq_range)` are retained.

Any keywords accepted by [`MTConfig`](@ref) may be passed here.
"""
function MTCrossSpectraConfig{T}(n_samples, n_channels, fs; demean=true, low_bias=true, freq_range=nothing, time_bandwidth_product=4, n_tapers = 2 * time_bandwidth_product - 1, kwargs...) where T
    if demean
        channelwise_mean = Matrix{T}(undef, 1, n_channels)
        demeaned_signal = Matrix{T}(undef, n_samples, n_channels)
    else
        channelwise_mean = nothing
        demeaned_signal = nothing
    end
    
    window = dpss(n_samples, time_bandwidth_product, n_tapers)
    evals = dpsseig(window, time_bandwidth_product)

    if low_bias
        keep_evals_mask = evals .> 0.9
        window, evals = window[:, keep_evals_mask], evals[keep_evals_mask]
    end
    n_tapers = size(window, 2)

    denom = fs * sum(evals) / 2
    weighted_evals = evals ./ denom


    mt_config = MTConfig{T}(n_samples, fs; window=window, n_tapers=n_tapers, time_bandwidth_product=time_bandwidth_product, kwargs...)

    x_mt = Array{fftouttype(T), 3}(undef, mt_config.n_freq_bins, mt_config.n_tapers, n_channels)

    if freq_range !== nothing
        freq_mask = first(freq_range) .< mt_config.freq .< last(freq_range)
        freq_inds = findall(freq_mask)
        freq = mt_config.freq[freq_mask]
    else
        freq_inds = eachindex(mt_config.freq)
        freq = mt_config.freq
    end
    return MTCrossSpectraConfig{T, typeof(weighted_evals), typeof(x_mt), typeof(channelwise_mean), typeof(demeaned_signal), typeof(freq), typeof(freq_range), typeof(freq_inds), typeof(mt_config)}(n_channels, weighted_evals, x_mt, demean, channelwise_mean, demeaned_signal, freq, freq_range, freq_inds,  mt_config)
end

allocate_output(config::MTCrossSpectraConfig{T}) where {T} = Array{fftouttype(T), 3}(undef, config.n_channels, config.n_channels, length(config.freq))

"""
    mt_cross_spectral!(csd_return, signal::AbstractMatrix{T}, config::MTCrossSpectraConfig{T}) where {T} -> CrossSpectral

Computes a multitapered cross spectral matrix.

* `output`: `n_channels` x `n_channels` x `length(config.freq)`. Can be created by `allocate_output(config)`.
* `signal`: `n_samples` x `n_channels`
* `config`: `MTCrossSpectraConfig{T}`

"""
@views function mt_cross_spectral!(output, signal::AbstractMatrix{T}, config::MTCrossSpectraConfig{T}) where {T}
    if config.demean
        mean!(config.channelwise_mean, signal)
        config.demeaned_signal .= signal .- config.channelwise_mean
        signal = config.demeaned_signal
    end

    x_mt = config.x_mt

    for k in axes(signal, 2)
        tapered_spectra!(x_mt[:,:,k], signal[:, k], config.mt_config)
    end
    x_mt[:, 1, :] ./= sqrt(2)
    if iseven(size(signal, 1))
        x_mt[:, end, :] ./= sqrt(2)
    end
    csd_inner!(output, config.weighted_evals, x_mt, config)
    return CrossSpectral(output, config.freq)
end

function csd_inner!(output, weighted_evals, x_mt, config)
    freq_inds = config.freq_inds
    n_channels = config.n_channels
    @boundscheck checkbounds(output, 1:n_channels, 1:n_channels, 1:length(freq_inds))
    @boundscheck checkbounds(weighted_evals, 1:length(weighted_evals))
    @boundscheck checkbounds(x_mt, freq_inds, 1:length(weighted_evals), 1:n_channels)

    # Up to the `weighted_evals` scaling, we have 
    # J_k^l(f) = x_mt[k, f, l]
    # Ŝ^lm(f) = output[l, m, f, epoch]
    # using the notation from https://en.wikipedia.org/wiki/Multitaper#The_method
    # so the formula `Ŝ^lm(f) = \sum_k conj(J_k^l(f)) * (J_k^m(f))`` becomes the following loop:
    @inbounds for (fi, f) in enumerate(freq_inds),
                  m in 1:n_channels,
                  l in 1:n_channels,
                  k in 1:length(weighted_evals)
        output[l, m, fi] += weighted_evals[k] *
                                      x_mt[f, k, l] *
                                      conj(x_mt[f, k, m])
    end
    return nothing
end

"""
    mt_cross_spectral(signal::AbstractMatrix{T}; fs, kwargs...) where {T} -> CrossSpectral

Compute a multitapered cross-spectral matrix.

* `signal`: `n_samples` x `n_channels`

Any keyword arguments accepted by [`MTCrossSpectraConfig`](@ref) may be passed here.
See also [`mt_cross_spectral!`](@ref).
"""
function mt_cross_spectral(signal::AbstractMatrix{T}; fs, kwargs...) where {T}
    n_samples, n_channels = size(signal)
    config = MTCrossSpectraConfig{T}(n_samples, n_channels, fs; kwargs...)
    output = allocate_output(config)
    return mt_cross_spectral!(output, signal, config)
end

#####
##### Multitapered coherences
#####

struct MTCoherenceConfig{T, T1,C<:MTCrossSpectraConfig{T}}
    cs_matrix::T1
    csd_config::C
end

"""
    MTCoherenceConfig(csd_config::MTCrossSpectraConfig{T}) where {T}
    MTCoherenceConfig{T}(n_samples, n_channels, fs; demean=true, low_bias=true, freq_range=nothing, kwargs...) where T

Creates a configuration object for coherences from a [`MTCrossSpectraConfig`](@ref). Provides a helper method
with the same arugments as `MTCrossSpectraConfig` to construct the `MTCrossSpectraConfig` object.
"""
function MTCoherenceConfig(csd_config::MTCrossSpectraConfig{T}) where {T}
    cs_matrix = allocate_output(csd_config)
    return MTCoherenceConfig{T, typeof(cs_matrix), typeof(csd_config)}(cs_matrix, csd_config)
end

function MTCoherenceConfig{T}(n_samples, n_channels, fs; demean=true, low_bias=true, freq_range=nothing, kwargs...) where T
    csd_config = MTCrossSpectraConfig{T}(n_samples, n_channels, fs; demean=demean, low_bias=low_bias, freq_range=freq_range, kwargs...)
    cs_matrix = allocate_output(csd_config)
    MTCoherenceConfig{T, typeof(cs_matrix), typeof(csd_config)}(cs_matrix, csd_config)
end

allocate_output(config::MTCoherenceConfig{T}) where T = Matrix{real(T)}(undef, config.csd_config.n_channels, config.csd_config.n_channels)

"""
    coherence_from_cs!(output, cs::CrossSpectral)

Compute the pairwise channel coherences from a cross spectral matrix, storing the result in `output`
and returning a `Symmetric` view of `output`.
"""
function coherence_from_cs!(output, cs::CrossSpectral)
    cs_matrix = cs.values
    n_channels = size(cs_matrix, 2)
    n_freqs = size(cs_matrix, 3)
    @boundscheck checkbounds(output, 1:n_channels, 1:n_channels)
    @boundscheck checkbounds(cs_matrix, 1:n_channels, 1:n_channels, 1:n_freqs)
    output .= 0
    @inbounds for ch2 in 1:n_channels
        for ch1 in (ch2 + 1):n_channels # lower triangular matrix
            for f in 1:n_freqs # average over frequency
                output[ch1, ch2] += real(abs(cs_matrix[ch1, ch2, f]) /
                                            (sqrt(cs_matrix[ch1, ch1, f] *
                                                cs_matrix[ch2, ch2, f]) * n_freqs))
            end
        end
    end
    return Symmetric(output, :L)
end

"""
    mt_coherence!(output, signal::AbstractMatrix{T}, config::MTCoherenceConfig{T}) where T

Computes the pairwise coherences between channels.

* `output`: `n_channels` x `n_channels` matrix 
* `signal`: `n_samples` x `n_channels` matrix
* `config`: configuration object that holds temporary variables.

Returns a `Symmetric` view of `output`.
"""
function mt_coherence!(output, signal::AbstractMatrix{T}, config::MTCoherenceConfig{T}) where T
    if size(signal) != (config.csd_config.mt_config.n_samples, config.csd_config.n_channels)
        throw(ArgumentError("Size of signal does not match `(config.csd_config.mt_config.n_samples, config.csd_config.n_channels)`;
            got `size(signal)`=$(size(signal)) but `(config.csd_config.mt_config.n_samples, config.csd_config.n_channels)`=$((config.csd_config.mt_config.n_samples, config.csd_config.n_channels))"))
    end
    if size(output) != (config.csd_config.n_channels, config.csd_config.n_channels)
        throw(ArgumentError("Size of output does not match `(config.csd_config.n_channels, config.csd_config.n_channels)`;
        got `size(output)`=$(size(output)) but `(config.csd_config.n_channels, config.csd_config.n_channels)`=$((config.csd_config.n_channels, config.csd_config.n_channels))"))
    end
    cs = mt_cross_spectral!(config.cs_matrix, signal, config.csd_config)
    return coherence_from_cs!(output, cs)
end

"""
    mt_coherence(signal::AbstractMatrix{T}; fs, freq_range = nothing, demean=true, low_bias=true, kwargs...) where T

Input: `signal`: `n_samples` x `n_channels` matrix

Output: `n_channels` x `n_channels` matrix of pairwise coherences between channels.

See [`MTCrossSpectraConfig`](@ref) for the meaning of the keyword arugments.
"""
function mt_coherence(signal::AbstractMatrix{T}; fs, freq_range = nothing, demean=true, low_bias=true, kwargs...) where T
    n_samples, n_channels = size(signal)
    config = MTCoherenceConfig{T}(n_samples, n_channels, fs; demean=demean, low_bias=low_bias, freq_range=freq_range, kwargs...)
    cohs = allocate_output(config)
    return mt_coherence!(cohs, signal, config)
end
