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
    frequencies_estimated = sort(esprit(x, M, p, Fs))
    @test isapprox(frequencies', frequencies_estimated; atol = 1e-2)
end

@testset "jacobsen" begin
    fs = 100
    t = range(0, 5; step=1/fs)
    function test_complex_jacobsen(fs, fc, f=0, t=t)
        sc = cispi.(2 * fc * t .+ f)
        f_est_complex = jacobsen(sc, fs)
        isapprox(f_est_complex, fc; atol=1e-5)
    end
    # test at two arbitrary frequencies
    @test test_complex_jacobsen(fs, -40.3, 1 / 1.4)
    @test test_complex_jacobsen(fs, 14.3, 1 / 3)
    # test near fs/2
    @test test_complex_jacobsen(fs, 49.90019)
    # test near -fs/2
    @test test_complex_jacobsen(fs, -49.90019)
    # test near +zero
    @test test_complex_jacobsen(fs, 0.04)
    # test near -zero
    @test test_complex_jacobsen(fs, -0.1)
    # tests for real signals: test only around fs/4, where the
    # expected error is small.
    fr = 28.3
    sr = cospi.(2 * fr * t .+ 1 / 4.2)
    f_est_real = jacobsen(sr, fs)
    @test isapprox(f_est_real, fr; atol=1e-5)
    fr = 23.45
    sr = sinpi.(2 * fr * t .+ 3 / 2.2)
    f_est_real = jacobsen(sr, fs)
    @test isapprox(f_est_real, fr; atol=1e-5)
end

@testset "quinn" begin
    ### real input
    fs = 100
    t = range(0, 5; step=1/fs)
    fr = 28.3
    sr = cospi.(2 * fr * t .+ 1 / 4.2)
    function test_quinn(f, s, args...)
        (f_est_real, maxiter) = quinn(s, args...)
        @test maxiter == false
        @test isapprox(f_est_real, f; atol=1e-3)
        return nothing
    end
    test_quinn(fr, sr, 50, fs)
    # use default initial guess
    test_quinn(fr, sr, fs) # initial guess given by Jacobsen
    # use default fs
    test_quinn(fr / fs, sr) # fs = 1.0, initial guess given by Jacobsen
    ### complex input
    fc = -40.3
    sc = cispi.(2 * fc * t .+ 1 / 1.4)
    test_quinn(fc, sc, -20, fs)
    # use default initial guess
    test_quinn(fc, sc, fs) # initial guess given by Jacobsen
    # use default fs
    test_quinn(fc / fs, sc) # fs = 1.0, initial guess by Jacobsen
end
