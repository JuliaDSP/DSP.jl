const epsilon = 10^-3

@testset "Configuration objects" begin
    @test_throws ArgumentError MTConfig{Float64}(100; nfft = 99)
    @test_throws ArgumentError MTConfig{Float64}(100; fs=-1)
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
        window = rand(n_samples, ntapers)
        for T in (Float32, Float64, Complex{Float32}, Complex{Float64})
            fft_input_tmp = Vector{T}(undef, nfft)
            onesided = T <: Real
            fft_flags = FFTW.ESTIMATE
            freqs = onesided ? rfftfreq(nfft, fs) : fftfreq(nfft, fs)
            fft_output_tmp = Vector{fftouttype(T)}(undef, length(freqs))
            r = fs*ntapers*ones(ntapers)
            plan = onesided ? plan_rfft(fft_input_tmp; flags=fft_flags) :
                plan_fft(fft_input_tmp; flags=fft_flags)
            # Test that the current configuration is valid so we know if it errors later
            # it's because we changed it, not that it was always broken
            @test MTConfig{T}(n_samples, nfft, ntapers, freqs, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r) isa MTConfig
            # now do a series of changes (in let-blocks to introduce new local bindings)
            # and check that they are each invalid
            let n_samples = 201
                @test_throws ArgumentError MTConfig{T}(n_samples, nfft, ntapers, freqs, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r)
            end
            let n_samples = -1
                @test_throws ArgumentError MTConfig{T}(n_samples, nfft, ntapers, freqs, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r)
            end
            let ntapers = -1
                @test_throws ArgumentError MTConfig{T}(n_samples, nfft, ntapers, freqs, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r)
            end
            let fs = -1
                @test_throws ArgumentError MTConfig{T}(n_samples, nfft, ntapers, freqs, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r)
            end
            let fft_input_tmp = Vector{T}(undef, 2*nfft)
                @test_throws DimensionMismatch MTConfig{T}(n_samples, nfft, ntapers, freqs, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r)
            end
            let nfft = 2*nfft, fft_input_tmp = Vector{T}(undef, 2*nfft)
                @test_throws DimensionMismatch MTConfig{T}(n_samples, nfft, ntapers, freqs, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r)
            end
            let window = rand(2*ntapers, n_samples)
                @test_throws DimensionMismatch MTConfig{T}(n_samples, nfft, ntapers, freqs, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r)
            end
            let r = 1.0
                @test_throws DimensionMismatch MTConfig{T}(n_samples, nfft, ntapers, freqs, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r)
            end
            if T <: Complex
                let onesided=true, freqs=rfftfreq(nfft, fs), fft_output_tmp=Vector{fftouttype(T)}(undef, length(freqs)), plan=plan_fft(fft_input_tmp; flags=fft_flags)
                    @test_throws ArgumentError MTConfig{T}(n_samples, nfft, ntapers, freqs, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r)
                end
            end
            # let's check that our original config is still valid, meaning we haven't
            # accidentally broken the config by not properly scoping our changes
            # (and therefore invalidating the tests that follow).
            @test MTConfig{T}(n_samples, nfft, ntapers, freqs, fs, plan, fft_input_tmp, fft_output_tmp, window, onesided, r) isa MTConfig
        end
    end

end

@testset "`dpss_config`" begin
    dpss_config = DSP.Periodograms.dpss_config(Float64, 1024; keep_only_large_evals=false, weight_by_evals=false)
    mt_config = MTConfig{Float64}(1024)

    @test dpss_config.window ≈ mt_config.window
    @test dpss_config.ntapers == mt_config.ntapers
    @test dpss_config.r ≈ mt_config.r

    dpss_config = DSP.Periodograms.dpss_config(Float64, 1024; keep_only_large_evals=false, weight_by_evals=true)
    @test dpss_config.window ≈ mt_config.window
    @test dpss_config.ntapers == mt_config.ntapers
    # should be different since we're weighting by eigenvalues now
    @test sum(abs, dpss_config.r - mt_config.r) > 0.1
end

avg_coh(x) = dropdims(mean(coherence(x); dims=3); dims=3)

@testset "Coherence (synthetic data)" begin
    fs = 1000.0
    n_samples = 1024
    t = (0:1023) ./ fs
    freq_range = (10, 15)
    sin_1 = sinpi.(2 * 12.0 * t)  # 12 Hz sinusoid signal
    sin_2 = sinpi.(2 * 12.0 * t .+ 1)
    noise = rand(1024) * 2 .- 1

    T = Float64
    same_signal = Matrix{T}(undef, 2, n_samples)
    same_signal[1, :] = sin_1
    same_signal[2, :] = sin_1
    coh = avg_coh(mt_coherence(same_signal; demean=true, fs, freq_range))
    same_signal_coherence = coh[2, 1]

    # test in-place gets the same result
    config = MTCoherenceConfig{T}(size(same_signal)...; demean=true, fs, freq_range)
    out = allocate_output(config)
    coh2 = avg_coh(mt_coherence!(out, same_signal, config))
    @test coh ≈ coh2

    # check that the output object has the right frequency information and properties
    coh_object = mt_coherence!(out, same_signal, config)
    @test freq(coh_object) == coh_object.freq == config.cs_config.freq
    @test all(x -> 10 <= x <= 15, freq(coh_object))
    # Check symmetries
    for i in 1:2, j in 1:2, f in eachindex(freq(coh_object))
        @test coherence(coh_object)[i, j, f] == coherence(coh_object)[j, i, f]
        if i == j
            @test coherence(coh_object)[i,j, f] == 1
        end
    end

    # check that manually demeaning with `demean=false` gives the same result
    same_signal_demeaned = same_signal .- mean(same_signal; dims = 2)
    coh3 = avg_coh(mt_coherence(same_signal_demeaned; demean = false))
    @test coh3 ≈ coh2


    @test_throws DimensionMismatch mt_coherence!(zeros(size(out,1), 2*size(out, 2), size(out, 3)), same_signal, config)
    @test_throws DimensionMismatch mt_coherence!(out, vcat(same_signal, same_signal), config)


    @test abs(same_signal_coherence - 1) < epsilon

    phase_shift = Matrix{T}(undef, 2, n_samples)
    phase_shift[1, :] = sin_1
    phase_shift[2, :] = sin_2
    coh = avg_coh(mt_coherence(phase_shift; fs, freq_range))
    phase_shift_coherence = coh[2,1]
    @test abs(phase_shift_coherence - 1) < epsilon
    @test coh[1,1] ≈ 1
    @test coh[2,2] ≈ 1

    # test in-place gets the same result
    config = MTCoherenceConfig{T}(size(phase_shift)...; fs, freq_range)
    out = allocate_output(config)
    coh2 = avg_coh(mt_coherence!(out, phase_shift, config))
    @test coh ≈ coh2

    # test out-of-place with config the same result
    config = MTCoherenceConfig{T}(size(phase_shift)...; fs, freq_range)
    coh3 = avg_coh(mt_coherence(phase_shift, config))
    @test coh ≈ coh3

    # test in-place without config the same result
    coh4 = avg_coh(mt_coherence!(out, phase_shift; fs, freq_range))
    @test coh ≈ coh4

    # Construct the config via an `MTConfig`
    mt_config = MTConfig{T}(size(phase_shift, 2); fs)
    config = MTCoherenceConfig(size(phase_shift, 1), mt_config; freq_range)
    coh5 = avg_coh(mt_coherence(phase_shift, config))
    @test coh ≈ coh5
    # test that including type parameter produces the same result
    config = @test_nowarn MTCoherenceConfig{T}(size(phase_shift, 1), mt_config; freq_range)
    cohT = avg_coh(mt_coherence(phase_shift, config))
    @test coh5 == cohT

    # Construct the config via an `MTCrossSpectraConfig`
    cs_config = MTCrossSpectraConfig(size(phase_shift, 1), mt_config; freq_range)
    # also test with type parameter
    config, configT = MTCoherenceConfig(cs_config), MTCoherenceConfig{T}(cs_config)
    coh5, cohT = avg_coh.(mt_coherence.((phase_shift,), (config, configT)))
    @test coh ≈ coh5 == cohT

    different_signal = Matrix{T}(undef, 2, n_samples)
    different_signal[1, :] = sin_1
    different_signal[2, :] = noise
    coh = avg_coh(mt_coherence(different_signal; fs, freq_range))
    different_signal_coherence = coh[2,1]
    # .8 is arbitrary, but represents a high coherence, so if a sine wave is
    # that coherent with noise, there is a problem
    @test different_signal_coherence < 0.8

    # test in-place gets the same result
    config = MTCoherenceConfig{T}(size(different_signal)...; fs, freq_range)
    out = allocate_output(config)
    coh2 = avg_coh(mt_coherence!(out, different_signal, config))
    @test coh ≈ coh2

    less_noisy = Matrix{T}(undef, 2, n_samples)
    less_noisy[1, :] = sin_1
    less_noisy[2, :] .= sin_1 .+ noise
    coh = avg_coh(mt_coherence(less_noisy; fs, freq_range))
    less_noisy_coherence = coh[2, 1]

    # test in-place gets the same result
    config = MTCoherenceConfig{T}(size(less_noisy)...; fs, freq_range)
    out = allocate_output(config)
    coh2 = avg_coh(mt_coherence!(out, less_noisy, config))
    @test coh ≈ coh2

    more_noisy = Matrix{T}(undef, 2, n_samples)
    more_noisy[1, :] = sin_1
    more_noisy[2, :] .= sin_1 .+ 3 * noise
    coh = avg_coh(mt_coherence(more_noisy; fs, freq_range))
    more_noisy_coherence = coh[2, 1]
    @test less_noisy_coherence < same_signal_coherence
    @test more_noisy_coherence < less_noisy_coherence
    @test different_signal_coherence < more_noisy_coherence

    # test in-place gets the same result
    config = MTCoherenceConfig{T}(size(more_noisy)...; fs, freq_range)
    out = allocate_output(config)
    coh2 = avg_coh(mt_coherence!(out, more_noisy, config))
    @test coh ≈ coh2

    several_signals = Matrix{T}(undef, 3, n_samples)
    several_signals[1, :] = sin_1
    several_signals[2, :] = sin_2
    several_signals[3, :] = noise
    several_signals_coherences = avg_coh(mt_coherence(several_signals; fs, freq_range))
    @test length(several_signals_coherences) == 9
    @test several_signals_coherences[2, 1] ≈ phase_shift_coherence
    @test several_signals_coherences[3, 1] ≈ different_signal_coherence

    # test in-place gets the same result
    config = MTCoherenceConfig{T}(size(several_signals)...; fs, freq_range)
    out = allocate_output(config)
    @test eltype(out) == Float64
    coh2 = avg_coh(mt_coherence!(out, several_signals, config))
    @test coh2 ≈ several_signals_coherences

    T = Float32
    # Float32 output:
    config = MTCoherenceConfig{T}(size(several_signals)...; fs, freq_range)
    out = allocate_output(config)
    @test eltype(out) == Float32
    coh3 = avg_coh(mt_coherence!(out, several_signals, config))
    @test coh3 ≈ several_signals_coherences

    # Float32 input and output:
    coh4 = avg_coh(mt_coherence!(out, Float32.(several_signals), config))
    @test coh4 ≈ several_signals_coherences
end

@testset "`mt_coherence` reference test" begin
    fs = 1000.0
    n_samples = 1024
    t = (0:1023) ./ fs
    sin_1 = sinpi.(2 * 12.0 * t)  # 12 Hz sinusoid signal
    sin_2 = sinpi.(2 * 12.0 * t .+ 1)
    noise = vec(readdlm(joinpath(@__DIR__, "data", "noise.txt"))) # generated by `rand(1024) * 2 .- 1`
    more_noisy = Array{Float64,3}(undef, 1, 2, n_samples)
    more_noisy[1, 1, :] = sin_1
    more_noisy[1, 2, :] .= sin_1 .+ 3 * noise

    ## to generate the reference result:
    # using PyMNE
    # mne_coherence_matrix, _ = PyMNE.connectivity.spectral_connectivity(more_noisy, method="coh",sfreq=fs,mode="multitaper",fmin=10,fmax=15,verbose=false)
    # coh = dropdims(mean(mne_coherence_matrix; dims=3); dims=3)[2, 1]
    coh = 0.982356762670818

    mt_config = DSP.Periodograms.dpss_config(Float64, n_samples; fs, keep_only_large_evals=true, weight_by_evals=true)
    config = MTCoherenceConfig(2, mt_config; freq_range = (10,15), demean=true)
    result = avg_coh(mt_coherence(dropdims(more_noisy;dims=1), config))
    @test result[2, 1] ≈ coh
end

@testset "`mt_cross_power_spectra`" begin
    fs = 1000.0
    n_samples = 1024
    t = (0:1023) ./ fs
    data = Array{Float64, 3}(undef, 1, 2, n_samples)
    data[1, 1, :] = sinpi.(2 * 12.0 * t) # 12 Hz sinusoid signal
    data[1, 2, :] = sinpi.(2 * 12.0 * t .+ 1)
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
    mt_config = DSP.Periodograms.dpss_config(Float64, n_samples; fs, keep_only_large_evals=true, weight_by_evals=true)
    config = MTCrossSpectraConfig(2, mt_config; demean=true)
    result = mt_cross_power_spectra(signal, config)
    @test signal isa Matrix{Float64}
    @test freq(result)[2:end] ≈ csd_array_multitaper_frequencies
    @test power(result)[:,:,2:end] ≈ csd_array_multitaper_values

    # Test in-place. Full precision:
    config = MTCrossSpectraConfig(size(signal, 1), mt_config; demean=true)
    out = allocate_output(config)
    @test eltype(out) == Complex{Float64}
    result2 = mt_cross_power_spectra!(out, signal, config)
    @test freq(result) ≈ freq(result2)
    @test power(result) ≈ power(result2)

    # Float32 output:
    mt_config32 = DSP.Periodograms.dpss_config(Float32, n_samples; fs, keep_only_large_evals=true, weight_by_evals=true)
    config = MTCrossSpectraConfig(size(signal, 1), mt_config32; demean=true)
    out = allocate_output(config)
    @test eltype(out) == Complex{Float32}
    result2 = mt_cross_power_spectra!(out, signal, config)
    @test freq(result2) ≈ freq(result)
    @test power(result2) ≈ power(result)

    # Float32 input and output:
    result3 = mt_cross_power_spectra!(out, Float32.(signal), config)
    @test freq(result3) ≈ freq(result)
    @test power(result3) ≈ power(result)

    @test_throws DimensionMismatch mt_cross_power_spectra!(similar(out, size(out, 1) + 1, size(out, 2), size(out, 3)), signal, config)
    @test_throws DimensionMismatch mt_cross_power_spectra!(out, vcat(signal, signal), config)

end


@testset "`mt_cross_power_spectra` agrees with `mt_pgram`" begin
    fs = 1000.0
    n_samples = 1024
    t = (0:1023) ./ fs
    noise = vec(readdlm(joinpath(@__DIR__, "data", "noise.txt")))
    signal = sinpi.(2 * 12.0 * t) .+ 3 * noise
    cs = mt_cross_power_spectra(reshape(signal, 1, :); fs)
    p = mt_pgram(signal; fs)

    @test freq(cs) ≈ freq(p)
    @test dropdims(power(cs); dims=(1,2)) ≈ power(p)

    # out-of-place with config
    config = MTCrossSpectraConfig{Float64}(1, length(signal); fs)
    cs = mt_cross_power_spectra(reshape(signal, 1, :), config)
    @test freq(cs) ≈ freq(p)
    @test dropdims(power(cs); dims=(1,2)) ≈ power(p)

    # in-place without config
    out = allocate_output(config)
    cs = mt_cross_power_spectra!(out, reshape(signal, 1, :); fs)
    @test freq(cs) ≈ freq(p)
    @test dropdims(power(cs); dims=(1,2)) ≈ power(p)

    # rm once two-sided FFTs supported
    @test_throws ArgumentError mt_cross_power_spectra(reshape(complex.(signal), 1, :); fs)
end
