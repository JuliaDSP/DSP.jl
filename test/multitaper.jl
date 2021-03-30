const epsilon = 10^-3

@testset "`MTConfig`" begin
    @test_throws ArgumentError MTConfig{Float64}(100, 1; n_for_fft = 99)
    @test_throws ArgumentError MTConfig{Float64}(100, -1)
    @test_throws ArgumentError MTConfig{Float64}(100, 1; n_tapers = -1)
    @test_throws DimensionMismatch MTConfig{Float64}(100, 1; window = rand(1, 200))
    @test_throws ArgumentError MTConfig{Complex{Float64}}(100, 1; onesided=true)
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
    same_signal_coherence = mt_coherence(same_signal; fs=fs, freq_range = (10, 15))[2, 1]
    @test abs(same_signal_coherence - 1) < epsilon

    phase_shift = Matrix{Float64}(undef, 2, n_samples)
    phase_shift[1, :] = sin_1
    phase_shift[2, :] = sin_2
    phase_shift_coherence = mt_coherence(phase_shift; fs=fs, freq_range = (10, 15))[2, 1]
    @test abs(phase_shift_coherence - 1) < epsilon

    different_signal = Matrix{Float64}(undef, 2, n_samples)
    different_signal[1, :] = sin_1
    different_signal[2, :] = noise
    different_signal_coherence = mt_coherence(different_signal; fs=fs, freq_range = (10, 15))[2, 1]
    # .8 is arbitrary, but represents a high coherence, so if a sine wave is
    # that coherent with noise, there is a problem
    @test different_signal_coherence < 0.8

    less_noisy = Matrix{Float64}(undef, 2, n_samples)
    less_noisy[1, :] = sin_1
    less_noisy[2, :] .= sin_1 .+ noise
    less_noisy_coherence = mt_coherence(less_noisy; fs=fs, freq_range = (10, 15))[2, 1]

    more_noisy = Matrix{Float64}(undef, 2, n_samples)
    more_noisy[1, :] = sin_1
    more_noisy[2, :] .= sin_1 .+ 3 * noise
    more_noisy_coherence = mt_coherence(more_noisy; fs=fs, freq_range = (10, 15))[2, 1]
    @test less_noisy_coherence < same_signal_coherence
    @test more_noisy_coherence < less_noisy_coherence
    @test different_signal_coherence < more_noisy_coherence

    several_signals = Matrix{Float64}(undef, 3, n_samples)
    several_signals[1, :] = sin_1
    several_signals[2, :] = sin_2
    several_signals[3, :] = noise
    several_signals_coherences = mt_coherence(several_signals; fs=fs, freq_range = (10, 15))
    @test length(several_signals_coherences) == 9
    @test several_signals_coherences[2, 1] ≈ phase_shift_coherence
    @test several_signals_coherences[3, 1] ≈ different_signal_coherence
end
