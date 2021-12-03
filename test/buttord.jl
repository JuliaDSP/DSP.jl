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

    (n, Wn) = buttord(40/500, 150/500, 3, 60, Lowpass)
    @test n == 5
    @test isapprox(Wn, 0.0810, rtol=1e-3)

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

    (nbp, Wnbp) = buttord([100/500, 200/500], [50/500, 250/500], 3, 40, Bandpass)
    @test nbp == 8
    @test isapprox(Wnbp[1], 0.1951, rtol=1e-3)
    @test isapprox(Wnbp[2], 0.4080, rtol=1e-3)
end