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
    # real input: test at two arbitrary frequencies
    fs = 100
    t = range(0, 5, step = 1/fs)
    fr = 28.3
    sr = cos.(2π*fr*t .+ π/4.2)
    f_est_real = jacobsen(sr, fs)
    @test isapprox(f_est_real, fr, atol = 1e-5)
    fr = 12.45
    sr = sin.(2π*fr*t .+ 3π/2.2)
    f_est_real = jacobsen(sr, fs)
    @test isapprox(f_est_real, fr, atol = 1e-5)
    # test at higher extreme of DFT
    fr = 49.9002
    sr = cos.(2π*fr*t)
    f_est_real = jacobsen(sr, fs)
    @test isapprox(f_est_real, fr, atol = 1e-5)
    # test at lower extreme of DFT
    @test jacobsen(zeros(10)) == 0.0
    # complex input: test at two arbitrary frequencies
    fc = -40.3
    sc = cis.(2π*fc*t .+ π/1.4)
    f_est_complex = jacobsen(sc, fs)
    @test isapprox(f_est_complex, fc, atol = 1e-5)
    fc = 14.3
    sc = cis.(2π*fc*t .+ π/3)
    f_est_complex = jacobsen(sc, fs)
    @test isapprox(f_est_complex, fc, atol = 1e-2)
    # test at higher extreme of DFT
    fr = 49.90019
    sr = cis.(2π*fr*t)
    f_est_real = jacobsen(sr, fs)
    @test isapprox(f_est_real, fr, atol = 1e-5)
    # test at lower extreme of DFT
    fr = -49.90019
    sr = cis.(2π*fr*t)
    f_est_real = jacobsen(sr, fs)
    @test isapprox(f_est_real, fr, atol = 1e-5)
end

@testset "quinn" begin
    # real input
    fs = 100
    t = range(0, 5, step = 1/fs)
    fr = 28.3
    sr = cos.(2π*fr*t .+ π/4.2)
    ### real input
    (f_est_real, maxiter) = quinn(sr, 50, fs)
    @test maxiter == false
    @test isapprox(f_est_real, fr, atol = 1e-3)
    # use default initial guess
    (f_est_real, maxiter) = quinn(sr, fs) # assumes f0 = 0.0
    @test maxiter == false
    @test isapprox(f_est_real, fr, atol = 1e-3)
    # use default fs
    (f_est_real, maxiter) = quinn(sr) # assumes fs = 1.0, f0 = 0.0
    @test maxiter == false
    @test isapprox(f_est_real, fr/fs, atol = 1e-3)
    ### complex input
    fc = -40.3
    sc = cis.(2π*fc*t .+ π/1.4)
    (f_est_real, maxiter) = quinn(sc, -20, fs)
    @test maxiter == false
    @test isapprox(f_est_real, fc, atol = 1e-3)
    # use default initial guess
    (f_est_real, maxiter) = quinn(sc, fs) # assumes f0 = 0.0
    @test maxiter == false
    @test isapprox(f_est_real, fc, atol = 1e-3)
    # use default fs
    (f_est_real, maxiter) = quinn(sc) # assumes fs = 1.0, f0 = 0.0
    @test maxiter == false
    @test isapprox(f_est_real, fc/fs, atol = 1e-3)
end
