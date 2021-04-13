const epsilon = 10^-3

@testset "Configuration objects" begin
    @test_throws ArgumentError MTConfig{Float64}(100; nfft = 99)
    @test_throws ArgumentError MTConfig{Float64}(100; fs=-1)
    @test_throws ArgumentError MTConfig{Float64}(100; ntapers = -1)
    @test_throws DimensionMismatch MTConfig{Float64}(100; window = rand(1, 200))
    @test_throws ArgumentError MTConfig{Complex{Float64}}(100; onesided=true)

    @testset "`MTSpectrogramConfig` with `n_samples`=$(n_samples), `n_samples_per_window`=$(n_samples_per_window), `n_overlap`=$(n_overlap) " for n_samples in (20:20:100), n_samples_per_window in (20:20:100), n_overlap in (0:20:(n_samples_per_window-1))
        config = MTSpectrogramConfig{Float64}(n_samples, n_samples_per_window, n_overlap)
        @test length(config.time) == length(arraysplit(1:n_samples, n_samples_per_window, n_overlap))
    end

    @testset "`MTConfig` inner constructor errors" begin
        n_samples = 100
        nfft = 200
        ntapers = 5
        fs = 1
        window = rand(ntapers, n_samples)
        for T in (Float32, Float64, Complex{Float32}, Complex{Float64})
            fft_input_tmp = Vector{T}(undef, nfft)
            onesided = T <: Real
            fft_flags = FFTW.ESTIMATE
            freq = onesided ? rfftfreq(nfft, fs) : fftfreq(nfft, fs)
            fft_output_tmp = Matrix{fftouttype(T)}(undef, length(freq), ntapers)
            r = 1.0
            plan = onesided ? plan_rfft(fft_input_tmp; flags=fft_flags) :
                plan_fft(fft_input_tmp; flags=fft_flags)
            let n_samples = 201
                @test_throws ArgumentError MTConfig{T}(n_samples, nfft, ntapers, freq, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r)
            end
            let n_samples = -1
                @test_throws ArgumentError MTConfig{T}(n_samples, nfft, ntapers, freq, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r)
            end
            let n_tapers = -1
                @test_throws DimensionMismatch MTConfig{T}(n_samples, nfft, ntapers, freq, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r)
            end
            let fs = -1
                @test_throws ArgumentError MTConfig{T}(n_samples, nfft, ntapers, freq, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r)
            end
            let fft_input_tmp = Vector{T}(undef, 2*nfft)
                @test_throws DimensionMismatch MTConfig{T}(n_samples, nfft, ntapers, freq, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r)
            end
            let nfft = 2*nfft
                fft_input_tmp = Vector{T}(undef, 2*nfft)
                @test_throws DimensionMismatch MTConfig{T}(n_samples, nfft, ntapers, freq, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r)
            end
            let  window = rand(2*ntapers, n_samples)
                @test_throws DimensionMismatch MTConfig{T}(n_samples, nfft, ntapers, freq, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r)
            end
        end
    end

end


@testset "Coherence (synthetic data)" begin
    fs = 1000.0
    n_samples = 1024
    t = (0:1023) ./ fs
    sin_1 = sin.(2 * π * 12.0 * t)  # 12 Hz sinusoid signal
    sin_2 = sin.(2 * π * 12.0 * t .+ π)
    noise = rand(1024) * 2 .- 1

    same_signal = Matrix{Float64}(undef, 2, n_samples)
    same_signal[1, :] = sin_1
    same_signal[2, :] = sin_1
    coh = mt_coherence(same_signal; fs=fs, freq_range = (10, 15))
    same_signal_coherence = coh[2,1]

    # test in-place gets the same result
    config = MTCoherenceConfig{eltype(same_signal)}(size(same_signal)...; fs=fs, freq_range = (10, 15))
    out = allocate_output(config)
    coh2 = mt_coherence!(out, same_signal, config)
    @test coh ≈ coh2

    # check that manually demeaning with `demean=false` gives the same result
    same_signal_demeaned = same_signal .- mean(same_signal; dims = 2)
    coh3 = mt_coherence(same_signal_demeaned; demean = false)
    @test coh3 ≈ coh2

    @test_throws DimensionMismatch mt_coherence!(zeros(size(out,1), 2*size(out, 2)), same_signal, config)
    @test_throws DimensionMismatch mt_coherence!(out, vcat(same_signal, same_signal), config)


    @test abs(same_signal_coherence - 1) < epsilon

    phase_shift = Matrix{Float64}(undef, 2, n_samples)
    phase_shift[1, :] = sin_1
    phase_shift[2, :] = sin_2
    coh = mt_coherence(phase_shift; fs=fs, freq_range = (10, 15))
    phase_shift_coherence = coh[2,1]
    @test abs(phase_shift_coherence - 1) < epsilon
    @test coh[1,1] ≈ 1
    @test coh[2,2] ≈ 1

    # test in-place gets the same result
    config = MTCoherenceConfig{eltype(phase_shift)}(size(phase_shift)...; fs=fs, freq_range = (10, 15))
    out = allocate_output(config)
    coh2 = mt_coherence!(out, phase_shift, config)
    @test coh ≈ coh2

    different_signal = Matrix{Float64}(undef, 2, n_samples)
    different_signal[1, :] = sin_1
    different_signal[2, :] = noise
    coh = mt_coherence(different_signal; fs=fs, freq_range = (10, 15))
    different_signal_coherence = coh[2,1]
    # .8 is arbitrary, but represents a high coherence, so if a sine wave is
    # that coherent with noise, there is a problem
    @test different_signal_coherence < 0.8

    # test in-place gets the same result
    config = MTCoherenceConfig{eltype(different_signal)}(size(different_signal)...; fs=fs, freq_range = (10, 15))
    out = allocate_output(config)
    coh2 = mt_coherence!(out, different_signal, config)
    @test coh ≈ coh2

    less_noisy = Matrix{Float64}(undef, 2, n_samples)
    less_noisy[1, :] = sin_1
    less_noisy[2, :] .= sin_1 .+ noise
    coh = mt_coherence(less_noisy; fs=fs, freq_range = (10, 15))
    less_noisy_coherence = coh[2, 1]

    # test in-place gets the same result
    config = MTCoherenceConfig{eltype(less_noisy)}(size(less_noisy)...; fs=fs, freq_range = (10, 15))
    out = allocate_output(config)
    coh2 = mt_coherence!(out, less_noisy, config)
    @test coh ≈ coh2

    more_noisy = Matrix{Float64}(undef, 2, n_samples)
    more_noisy[1, :] = sin_1
    more_noisy[2, :] .= sin_1 .+ 3 * noise
    coh = mt_coherence(more_noisy; fs=fs, freq_range = (10, 15))
    more_noisy_coherence = coh[2, 1]
    @test less_noisy_coherence < same_signal_coherence
    @test more_noisy_coherence < less_noisy_coherence
    @test different_signal_coherence < more_noisy_coherence

    # test in-place gets the same result
    config = MTCoherenceConfig{eltype(more_noisy)}(size(more_noisy)...; fs=fs, freq_range = (10, 15))
    out = allocate_output(config)
    coh2 = mt_coherence!(out, more_noisy, config)
    @test coh ≈ coh2

    several_signals = Matrix{Float64}(undef, 3, n_samples)
    several_signals[1, :] = sin_1
    several_signals[2, :] = sin_2
    several_signals[3, :] = noise
    several_signals_coherences = mt_coherence(several_signals; fs=fs, freq_range = (10, 15))
    @test length(several_signals_coherences) == 9
    @test several_signals_coherences[2, 1] ≈ phase_shift_coherence
    @test several_signals_coherences[3, 1] ≈ different_signal_coherence

    # test in-place gets the same result
    config = MTCoherenceConfig{Float64}(size(several_signals)...; fs=fs, freq_range = (10, 15))
    out = allocate_output(config)
    @test eltype(out) == Float64
    coh2 = mt_coherence!(out, several_signals, config)
    @test coh2 ≈ several_signals_coherences

    # Float32 output:
    config = MTCoherenceConfig{Float32}(size(several_signals)...; fs=fs, freq_range = (10, 15))
    out = allocate_output(config)
    @test eltype(out) == Float32
    coh3 = mt_coherence!(out, several_signals, config)
    @test coh3 ≈ several_signals_coherences

    # Float32 input and output:
    coh4 = mt_coherence!(out, Float32.(several_signals), config)
    @test coh4 ≈ several_signals_coherences
end

@testset "`mt_coherence` reference test" begin
    fs = 1000.0
    n_samples = 1024
    t = (0:1023) ./ fs
    sin_1 = sin.(2 * π * 12.0 * t)  # 12 Hz sinusoid signal
    sin_2 = sin.(2 * π * 12.0 * t .+ π)
    noise = vec(readdlm(joinpath(@__DIR__, "data", "noise.txt"))) # generated by `rand(1024) * 2 .- 1`
    more_noisy = Array{Float64}(undef, 1, 2, n_samples)
    more_noisy[1, 1, :] = sin_1
    more_noisy[1, 2, :] .= sin_1 .+ 3 * noise

    ## to generate the reference result:
    # using PyMNE
    # mne_coherence_matrix, _ = PyMNE.connectivity.spectral_connectivity(more_noisy, method="coh", sfreq=fs,mode="multitaper",fmin=10,fmax=15,verbose=false)
    # coh = dropdims(mean(mne_coherence_matrix; dims=3); dims=3)[2, 1]
    coh = 0.982356762670818
    result = mt_coherence(dropdims(more_noisy;dims=1); fs=fs, freq_range = (10,15))
    @test result[2, 1] ≈ coh
end

@testset "`mt_cross_spectral`" begin
    fs = 1000.0
    n_samples = 1024
    t = (0:1023) ./ fs
    data = Array{Float64, 3}(undef, 1, 2, n_samples)
    data[1,1,:] = sin.(2 * π * 12.0 * t) # 12 Hz sinusoid signal
    data[1,2,:] = sin.(2 * π * 12.0 * t .+ π)
    ## Script to generate reference data:
    # using PyMNE, DelimitedFiles
    # result = PyMNE.time_frequency.csd_array_multitaper(data, fs, n_fft = nextpow(2, n_samples), low_bias=true, adaptive=false)
    # csd_array_multitaper_frequencies = result.frequencies
    # csd_array_multitaper_values = reduce((x,y) -> cat(x,y; dims=3), [ r.get_data() for r in result])
    # writedlm("data/csd_array_multitaper_frequencies.txt", csd_array_multitaper_frequencies)
    # writedlm("data/csd_array_multitaper_values_re.txt", real(csd_array_multitaper_values))
    # writedlm("data/csd_array_multitaper_values_im.txt", imag(csd_array_multitaper_values))
    csd_array_multitaper_frequencies = readdlm(joinpath(@__DIR__, "data", "csd_array_multitaper_frequencies.txt"))
    csd_array_multitaper_values_re = reshape(readdlm(joinpath(@__DIR__, "data", "csd_array_multitaper_values_re.txt")), (2,2,512))
    csd_array_multitaper_values_im = reshape(readdlm(joinpath(@__DIR__, "data", "csd_array_multitaper_values_im.txt")), (2,2,512))
    csd_array_multitaper_values = csd_array_multitaper_values_re + im*csd_array_multitaper_values_im

    signal = dropdims(data; dims=1)
    @test signal isa Matrix{Float64}
    result = mt_cross_spectral(signal; fs=fs)
    @test freq(result)[2:end] ≈ csd_array_multitaper_frequencies
    @test result.values[:,:,2:end] ≈ csd_array_multitaper_values

    # Test in-place. Full precision:
    config = MTCrossSpectraConfig{Float64}(size(signal)...; fs=fs)
    out = allocate_output(config)
    @test eltype(out) == Complex{Float64}
    result2 = mt_cross_spectral!(out, signal, config)
    @test freq(result) ≈ freq(result2)
    @test result.values ≈ result2.values

    # Float32 output:
    config = MTCrossSpectraConfig{Float32}(size(signal)...; fs=fs)
    out = allocate_output(config)
    @test eltype(out) == Complex{Float32}
    result2 = mt_cross_spectral!(out, signal, config)
    @test freq(result2) ≈ freq(result)
    @test result2.values ≈ result.values

    # Float32 input and output:
    result3 = mt_cross_spectral!(out, Float32.(signal), config)
    @test freq(result3) ≈ freq(result)
    @test result3.values ≈ result.values

    @test_throws DimensionMismatch mt_cross_spectral!(similar(out, size(out, 1) + 1, size(out, 2), size(out, 3)), signal, config)
    @test_throws DimensionMismatch mt_cross_spectral!(out, vcat(signal, signal), config)

end
