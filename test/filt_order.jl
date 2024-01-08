using DSP, Test


Δ = 1e-3
@testset "buttord" begin

    #
    # https://www.mathworks.com/help/signal/ref/buttord.html#d123e9937
    # Example 1: Lowpass Butterworth Example
    #
    # Wp = 40/500;
    # Ws = 150/500;
    # [n, Wn] = buttord(Wp, Ws, 3, 60)
    #

    # z-domain test
    (n, Wn) = buttord(40 / 500, 150 / 500, 3, 60, domain=:z)
    @test n == 5
    @test Wn ≈ 0.081038494957764

    # s-domain test
    (ns, Wns) = buttord(40 / 500, 150 / 500, 3, 60, domain=:s)
    @test ns == 6
    @test Wns ≈ 0.0948683377107

    #
    # Highpass filter example
    # Wp = 1200/2000; Ws=600/2000;
    # Rs = 3; Rp = 60;

    # z-domain (default)
    (nhpf, Wnhpf) = buttord(1200 / 2000, 600 / 2000, 3, 60, domain=:z)
    @test nhpf == 7
    @test Wnhpf ≈ 0.597905417809


    # s-domain test
    (nhpfs, Wnhpfs) = buttord(1200 / 2000, 600 / 2000, 3, 60, domain=:s)
    @test nhpfs == 10
    @test Wnhpfs ≈ 0.598578664562

    #
    # https://www.mathworks.com/help/signal/ref/buttord.html#d123e9937
    # Example 2: Bandpass Butterworth Example
    #
    # Wp = [100 200]/500;
    # Ws = [50 250]/500;
    # Rp = 3;
    # Rs = 40;
    # [n, Wn] = buttord(Wp, Ws, Rp, Rs)
    #

    (nbp, Wnbp) = buttord((100 / 500, 200 / 500), (50 / 500, 250 / 500), 3, 40, domain=:z)
    @test nbp == 8
    @test Wnbp[1] ≈ 0.195101359239
    @test Wnbp[2] ≈ 0.408043633382

    # s-domain
    (nbps, Wnbps) = buttord((100 / 500, 200 / 500), (50 / 500, 250 / 500), 3, 40, domain=:s)
    @test nbps == 9
    @test Wnbps[1] ≈ 0.198730150231
    @test Wnbps[2] ≈ 0.402555927759

    #
    # Bandstop Example, (44.1 kHz Nyquist)
    # Wp = [3200 7800]/22050
    # Ws = [4800 5600]/22050
    # Rp = 2
    # Rs = 60

    # this test may be more sensitive...MATLAB's implementation of bounded minimization
    # will yield different results in comparison to Optim.jl.
    (nbs, Wnbs) = buttord((3200 / 22050, 7800 / 22050), (4800 / 22050, 5600 / 22050), 2, 60, domain=:z)
    @test nbs == 5
    @test ≈(Wnbs[1], 0.172660908966, rtol=Δ)
    @test ≈(Wnbs[2], 0.314956388749, rtol=Δ)

    # s-domain
    (nbss, Wnbss) = buttord((3200 / 22050, 7800 / 22050), (4800 / 22050, 5600 / 22050), 2, 60, domain=:s)
    @test ≈(Wnbss[1], 0.173677826752, rtol=Δ)
    @test ≈(Wnbss[2], 0.318267164272, rtol=Δ)

end

#
# Elliptic/Chebyshev I/II Filter cases tested against SciPy signal results.
#
@testset "ellipord" begin

    Rp, Rs = 3, 40
    Wp = (0.2, 0.7)
    Ws = (0.1, 0.8)

    # Lowpass
    (nl, Wnl) = ellipord(0.1, 0.2, Rp, Rs, domain=:z)
    @test nl == 3
    @test Wnl == 0.1

    # Highpass
    (nh, Wnh) = ellipord(0.3, 0.1, Rp, Rs, domain=:z)
    @test nh == 3
    @test Wnh == 0.3

    # Bandpass (z-domain)
    (nbp, Wnbp) = ellipord(Wp, Ws, Rp, Rs, domain=:z)
    @test nbp == 4
    @test Wnbp == Wp

    # Bandpass (s-domain)
    (nbp, Wnbp) = ellipord(Wp, Ws, Rp, Rs, domain=:s)
    @test nbp == 5
    @test Wnbp == Wp

    # Bandstop (z-domain)
    (nbs, Wnbs) = ellipord(Ws, Wp, Rp, Rs, domain=:z)
    @test nbs == 4
    @test Wnbs == Ws

    # Bandstop (s-domain)
    # n, Wn = scipy.signal.ellipord([0.1, 0.8], [0.2, 0.7], 3, 40, True)
    (nbs, Wnbs) = ellipord(Ws, Wp, Rp, Rs, domain=:s)
    @test nbs == 5
    @test ≈(Wnbs[1], 0.17500000332788998, rtol=Δ)
    @test ≈(Wnbs[2], 0.799993389303865, rtol=Δ)

end

@testset "cheb1ord" begin
    Rp, Rs = 2, 70
    Wp = (0.2, 0.5)
    Ws = (0.1, 0.6)

    # Lowpass
    (nl, Wnl) = cheb1ord(0.1, 0.21, Rp, Rs, domain=:z)
    @test nl == 7
    @test Wnl == 0.1

    # Highpass
    (nh, Wnh) = cheb1ord(0.1, 0.04, Rp, Rs, domain=:z)
    @test nh == 6
    @test Wnh == 0.1

    # Bandpass (z-domain)
    (nbp, Wnbp) = cheb1ord(Wp, Ws, Rp, Rs, domain=:z)
    @test nbp == 9
    @test Wnbp == (0.2, 0.5)


    # Bandpass (s-domain)
    (nbp, Wnbp) = cheb1ord(Wp, Ws, Rp, Rs, domain=:s)
    @test nbp == 10
    @test Wnbp == (0.2, 0.5)

    # Bandstop (z-domain)
    (nbs, Wnbs) = cheb1ord(Ws, Wp, Rp, Rs, domain=:z)
    @test nbs == 9
    @test Wnbs == Ws

    # Bandstop (s-domain)
    (nbs, Wnbs) = cheb1ord(Ws, Wp, Rp, Rs, domain=:s)
    @test nbs == 10
    @test ≈(Wnbs[1], 0.166666612185443, rtol=Δ)
    @test ≈(Wnbs[2], 0.5999933893038649, rtol=Δ)
end

@testset "cheb2ord" begin
    Rp, Rs = 1.2, 80
    Wp = (0.22, 0.51)
    Ws = (0.14, 0.63)

    # Lowpass
    (nl, Wnl) = cheb2ord(0.1, 0.21, Rp, Rs, domain=:z)
    @test nl == 8
    @test Wnl == 0.19411478246577737
    (nl, Wnl) = cheb2ord(0.1, 0.21, Rp, Rs, domain=:s)
    @test nl == 8
    @test Wnl == 0.1987124302811051

    # Highpass
    (nh, Wnh) = cheb2ord(0.21, 0.1, Rp, Rs, domain=:z)
    @test nh == 8
    @test Wnh == 0.10862150541420543
    (nh, Wnh) = cheb2ord(0.21, 0.1, Rp, Rs, domain=:s)
    @test nh == 8
    @test Wnh == 0.10568035411923006

    # Bandpass
    (nbp, Wnbp) = cheb2ord(Wp, Ws, Rp, Rs, domain=:z)
    @test nbp == 9
    @test ≈(Wnbp[1], 0.1608459041132262)
    @test ≈(Wnbp[2], 0.6133747025904719)
    (nbp, Wnbp) = cheb2ord(Wp, Ws, Rp, Rs, domain=:s)
    @test nbp == 11
    @test ≈(Wnbp[1], 0.18262279523940905)
    @test ≈(Wnbp[2], 0.6143811338169016)

    # Bandstop
    (nbs, Wnbs) = cheb2ord(Ws, Wp, Rp, Rs, domain=:z)
    @test nbs == 9
    @test ≈(Wnbs[1], 0.21211425852327126, rtol=Δ)
    @test ≈(Wnbs[2], 0.5225427194473862, rtol=Δ)
    (nbs, Wnbs) = cheb2ord(Ws, Wp, Rp, Rs, domain=:s)
    @test nbs == 11
    @test ≈(Wnbs[1], 0.2159740591083134, rtol=Δ)
    @test ≈(Wnbs[2], 0.5195028184932494, rtol=Δ)

end

# using some simple examples for testing Brent's method.
@testset "brent" begin
    f1(x) = (x + 3) * ((x - 1)^2) # x³ + x² - 5x + 3
    @test ≈(f1(DSP.Filters.brent(f1, -4.0, 4.0)), 0.0, atol=1e-8)
    @test ≈(sin(DSP.Filters.brent(sin, 0.0, 2pi)), -1.0, atol=1e-8)
    @test ≈(cos(DSP.Filters.brent(cos, 0.0, 2pi)), -1.0, atol=1e-8)
end


@testset "remezord" begin
    #
    # Using the test-cases highlighted in [^Parks] Figures 8 and 15.
    #
    # [^Parks]: Rabiner, L. R., McClellan, J. H., & Parks, T. W. (1975). 
    #   FIR digital filter design techniques using weighted Chebyshev 
    #   approximation. Proceedings of the IEEE, 63(4), 595-610.
    #
    @test remezord(0.41665, 0.49417, 0.0116, 0.0001) == 39
    @test remezord(0.135, 0.203, 0.1, 0.1) == 10
    @test remezord(0.203, 0.135, 0.1, 0.1) == 10 # flipped HPF case.
    @test remezord(0.00963, 0.13271, 0.0116, 0.0001) == 24
end