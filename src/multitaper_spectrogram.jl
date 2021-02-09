struct MultitaperedSpectrogramConfig{T, P, F}
    n_samples::Int
    sample_rate::Int
    samples_per_window::Int
    n_overlap_samples::Int
    n_for_fft::Int
    n_time_points::Int
    n_freq_bins::Int
    time_bandwidth_product::Int
    n_tapers::Int
    dpss_window::Matrix{T}
    input_tmp::Vector{T}
    tmp::Vector{Complex{T}}
    plan::P
    freq::F
    time::FloatRange{Float64}
    function MultitaperedSpectrogramConfig{T}(; n_samples, sample_rate, samples_per_window,
                               n_overlap_samples,
                               n_for_fft = nextpow(2, samples_per_window)
                               time_bandwidth_product=4,
                               ) where T

        # Dimensions of the spectrogram power matrix
        n_time_points = div(n_samples - samples_per_window,
                            samples_per_window - n_overlap_samples) + 1

        freq = rfftfreq(n_for_fft, sample_rate)
        n_freq_bins = length(freq)

        n_tapers = 2 * time_bandwidth_product - 1
        
        input_tmp = zeros(T, n_for_fft)
        tmp = zeros(Complex{T}, n_freq_bins)
        plan = plan_rfft(input_tmp)
        dpss_window = dpss(samples_per_window, time_bandwidth_product, n_tapers)

        window_duration = samples_per_window / sample_rate
        overlap_duration = n_overlap_samples / sample_rate

        hop = window_duration - overlap_duration
        time = range(window_duration/2, step=hop, length = n_time_points)

        return new{T, typeof(plan), typeof(freq)}(n_samples, sample_rate, samples_per_window, n_overlap_samples,
                   n_for_fft, n_time_points, n_freq_bins, time_bandwidth_product, n_tapers,
                   dpss_window, input_tmp, tmp, plan, freq, time)
    end
end

function multitapered_spectrogram!(destination::AbstractMatrix{T}, signal::AbstractVector{T},
                                   config::MultitaperedSpectrogramConfig) where {T <: Real}
    expected_destination_size = (config.n_time_points, config.n_freq_bins)
    if size(destination) != expected_destination_size
        throw(ArgumentError("size(destination) == $(size(destination)) != $expected_destination_size == (config.n_time_points, config.n_freq_bins)"))
    end
    expected_signal_length = config.n_samples
    if length(signal) != expected_signal_length
        throw(ArgumentError("length(signal) == $(length(signal)) != $expected_signal_length == config.n_samples"))
    end

    subepochs = arraysplit(signal, config.samples_per_window, config.n_overlap_samples)
    for (time_index, subepoch) in enumerate(subepochs)
        power_mt_pgram!(view(destination, time_index, :), subepoch, config)
    end

    return Spectrogram(destination, config.freq, config.time)
end

function multitapered_spectrogram(signal::AbstractVector{T},
                                  n::Int=length(s)>>3, n_overlap::Int=n>>1; 
                                  fs::Int=1) where {T <: Real}

    config = MultitaperedSpectrogramConfig{T}(;
        n_samples=length(signal),
        sample_rate=fs,
        samples_per_window = n,
        n_overlap_samples = n_overlap)

    X = Array{T,2}(undef, config.n_time_points, config.n_freq_bins)
    return multitapered_spectrogram!(X, signal, config)
end

# Simplified, in-place version of
# ```
# power(mt_pgram(signal; fs=config.sample_rate, nfft=config.n_for_fft))
# ```
# Benchmarks suggest that this uses a third of the memory, a tenth of the allocations,
# and executes in a tenth of the time. It also allows us to avoid allocating an output
# array and copying it to the destination.
function power_mt_pgram!(output, signal, config::MultitaperedSpectrogramConfig{T}) where T
    @assert length(signal) == size(config.dpss_window, 1)
    fill!(output, zero(T))
    input = config.input_tmp
    plan = config.plan
    temp = config.tmp
    r = config.sample_rate * config.n_tapers
    m = iseven(config.n_for_fft) ? 1 : 2
    @inbounds for j in 1:(config.n_tapers)
        for i in 1:length(signal)
            input[i] = config.dpss_window[i, j] * signal[i]
        end
        mul!(temp, plan, input)
        output[1] += abs2(temp[1]) / r
        for i in 2:(config.n_freq_bins - 1)
            output[i] += 2 * abs2(temp[i]) / r
        end
        output[end] += m * abs2(temp[end]) / r
    end
    return output
end
