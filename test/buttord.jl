using DSP, Test

@testset "buttord" begin
    
    #
    # https://www.mathworks.com/help/signal/ref/buttord.html#d123e9937
    # Example 1: Lowpass Butterworth Example
    #
    # Wp = 40/500;
    # Ws = 150/500;
    # [n, Wn] = buttord(Wp, Ws, 3, 60)
    #

    (n, Wn) = buttord(40/500, 150/500, 3, 60)
    @test n == 5
    @test isapprox(Wn, 0.081038494957764, rtol=1e-12)

    #
    # Highpass filter example
    # Wp = 600/2000; Ws = 1200/2000;
    # Rs = 3; Rp = 60;

    (nhpf, Wnhpf) = buttord(600/2000, 1200/2000, 3, 60)
    @test nhpf == 7
    @test isapprox(Wnhpf, 0.301783479785, rtol=1e-12)


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

    (nbp, Wnbp) = buttord([100/500, 200/500], [50/500, 250/500], 3, 40)
    @test nbp == 8
    @test isapprox(Wnbp[1], 0.195101359239, rtol=1e-12)
    @test isapprox(Wnbp[2], 0.408043633382, rtol=1e-12)

    #
    # Bandstop Example, (44.1 kHz Nyquist)
    # Wp = [3200 7800]/22050
    # Ws = [4800 5600]/22050
    # Rp = 2
    # Rs = 60

    # this test may be more sensitive...MATLAB's implementation of bounded minimization
    # will yield different results in comparison to Optim.jl.
    (nbs, Wnbs) = buttord([3200/22050, 7800/22050], [4800/22050, 5600/22050], 2, 60)
    @test nbs == 5
    @test isapprox(Wnbs[1], 0.172660908966, rtol=1e-12)
    @test isapprox(Wnbs[2], 0.314956388749, rtol=1e-12)

end