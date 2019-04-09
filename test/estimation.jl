using DSP, Test

@testset "esprit" begin
    # create a sum of sinusoids in noise, and estimate their frequencies
    Fs = 10000.0           # sampling frequency in Hz
    duration = 1           # length of signal, in seconds
    n = Int(Fs * duration) # number of samples
    t = collect((1:n)/Fs)  # time vector
    frequencies = [1000.0 1250.0]
    amplitudes = [2.0 1.5]
    phases = [0.7 -1.0]
    x = exp.( 1im*2*π*t*frequencies .+ phases) * amplitudes'
    noise = randn(n, 2)*[1;1im]
    sigma = 0.1
    noise *= sigma
    x += noise
    M = 300
    p = 2                  # number of sinusoids to estimate
    frequencies_estimated = esprit(x, M, p, Fs)
    @test isapprox(frequencies', frequencies_estimated; atol = 1e-2)
end

