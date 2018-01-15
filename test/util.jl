using DSP, Compat, Compat.Test

@testset "unwrap" begin
    @test unwrap([0.1, 0.2, 0.3, 0.4]) ≈ [0.1, 0.2, 0.3, 0.4]
    @test unwrap([0.1, 0.2 + 2pi, 0.3, 0.4]) ≈ [0.1, 0.2, 0.3, 0.4]
    @test unwrap([0.1, 0.2 - 2pi, 0.3, 0.4]) ≈ [0.1, 0.2, 0.3, 0.4]
    @test unwrap([0.1, 0.2 - 2pi, 0.3 - 2pi, 0.4]) ≈ [0.1, 0.2, 0.3, 0.4]
    @test unwrap([0.1 + 2pi, 0.2, 0.3, 0.4]) ≈ [0.1 + 2pi, 0.2 + 2pi, 0.3 + 2pi, 0.4 + 2pi]
    @test unwrap([0.1, 0.2 + 6pi, 0.3, 0.4]) ≈ [0.1, 0.2, 0.3, 0.4]

    test_v = [0.1, 0.2, 0.3 + 2pi, 0.4]
    res_v = unwrap(test_v)
    @test test_v ≈ [0.1, 0.2, 0.3 + 2pi, 0.4]
    unwrap!(test_v)
    @test test_v ≈ [0.1, 0.2, 0.3, 0.4]

    # test multi-dimensional unwrapping
    wrapped = [0.1, 0.2 + 2pi, 0.3, 0.4]
    unwrapped = [0.1, 0.2, 0.3, 0.4]
    wrapped = hcat(wrapped, wrapped)
    unwrapped = hcat(unwrapped, unwrapped)
    @test unwrap(wrapped) ≈ wrapped
    @test unwrap(wrapped, 1) ≈ unwrapped

    # test unwrapping with other ranges
    unwrapped = [1.0:100;]
    wrapped = Float64[i % 10 for i in unwrapped]
    @test unwrap(wrapped, range=10) ≈ unwrapped
end

@testset "hilbert" begin
    # Testing hilbert transform
    t = (0:1/256:2-1/256)
    a0 = sinpi.(t)
    a1 = cospi.(t)
    a2 = sinpi.(2*t)
    a3 = cospi.(2*t)
    a = hcat(a0, a1, a2, a3)

    h = hcat(hilbert(a0), hilbert(a1), hilbert(a2), hilbert(a3))
    h_abs = abs.(h)
    h_angle = angle.(h)
    h_real = real.(h)

    #The real part should be equal to the original signals:
    @test h_real ≈ a

    #The absolute value should be one everywhere, for this input:
    @test h_abs ≈ ones(size(a))

    #For the 'slow' sine - the phase should go from -pi/2 to pi/2 in
    #the first 256 bins:
    @test h_angle[1:256,1] ≈ -pi/2:pi/256:pi/2-pi/256

    #For the 'slow' cosine - the phase should go from 0 to pi in the
    #same interval:
    @test h_angle[1:256,2] ≈ 0:pi/256:pi-pi/256

    #The 'fast' sine should make this phase transition in half the time:
    @test h_angle[1:128,3] ≈ -pi/2:pi/128:pi/2-pi/128

    #Ditto for the 'fast' cosine:
    @test h_angle[1:128,4] ≈ 0:pi/128:pi-pi/128

    #The imaginary part of hilbert(cos(t)) = sin(t) Wikipedia
    @test imag(h[:,2]) ≈ a0

    #Sanity check with odd number of samples
    h2 = hilbert([ones(10); zeros(9)])
    @test real(h2) ≈ [ones(10); zeros(9)]

    #Sanity check with integer arguments
    r = rand(1:20, 128)
    @test hilbert(r) == hilbert(map(Float64, r))

    # Test hilbert with 2D input
    @test h ≈ hilbert(a)
end

@testset "fft helpers" begin
    ## FFTFREQ
    @test fftfreq(1) ≈ [0.]
    @test fftfreq(2) ≈ [0., -1/2]
    @test fftfreq(2, 1/2) ≈ [0., -1/4]
    @test fftfreq(3) ≈ [0., 1/3, -1/3]
    @test fftfreq(3, 1/2) ≈ [0., 1/6, -1/6]
    @test fftfreq(6) ≈ [0., 1/6, 1/3, -1/2, -1/3, -1/6]
    @test fftfreq(7) ≈ [0., 1/7, 2/7, 3/7, -3/7, -2/7, -1/7]

    @test rfftfreq(1) ≈ [0.]
    @test rfftfreq(2) ≈ [0., 1/2]
    @test rfftfreq(2, 1/2) ≈ [0., 1/4]
    @test rfftfreq(3) ≈ [0., 1/3]
    @test rfftfreq(3, 1/2) ≈ [0., 1/6]
    @test rfftfreq(6) ≈ [0., 1/6, 1/3, 1/2]
    @test rfftfreq(7) ≈ [0., 1/7, 2/7, 3/7]

    for n = 1:7
        @test fftshift(fftfreq(n)) ≈ fftshift([fftfreq(n);])
    end

    @test meanfreq(sin.(2*π*10*(0:1e-3:10*π)),1e3) ≈ 10.0 rtol=1e-3

    # nextfastfft
    @test nextfastfft(64) == 64
    @test nextfastfft(65) == 70
    @test nextfastfft(127) == 128
    @test nextfastfft((64,65,127)) == (64,70,128)
    @test nextfastfft(64,65,127) == nextfastfft((64,65,127))
end

## COMMON DSP TOOLS

@testset "dB conversion" begin
    # dB conversion
    @test 3dB ≈ db2pow(3)
    @test -3dB ≈ db2pow(-3)
    @test 3dBa ≈ db2amp(3)
    @test -3dBa ≈ db2amp(-3)
    @test isa(3e0dB, Float64)
    @test isa(3f0dB, Float32)
    test_num = convert(Float64, pi)
    @test pow2db(test_num) ≈ 10*log10(test_num)
    @test amp2db(test_num) ≈ 20*log10(test_num)
    @test test_num*dB ≈ db2pow(test_num)
    @test test_num*dBa ≈ db2amp(test_num)
    @test test_num ≈ db2pow(pow2db(test_num))
    @test test_num ≈ db2amp(amp2db(test_num))
end

@testset "rms" begin
    n = (5,6)
    for x in ( randn(n), randn(n)+randn(n)im )
        @test rms(x) ≈ sqrt(mean(abs.(x).^2))
        @test rmsfft(fft(x)) ≈ rms(x)
    end # for
end

@test shiftin!([1,2,3,4],[5,6]) ≈ [3,4,5,6]
