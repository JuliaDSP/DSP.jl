using DSP
using Test
using DelimitedFiles: readdlm

@testset "rational ratio" begin
    # AM Modulator
    # sig(t) = [(1 + sin(2π*0.005*t)) * sin(2π*.05*t) for t in t]
    x_ml   = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_x.txt"),'\t'))

    #
    # [y1,b1] = resample(x, 1, 2)
    #
    rate   = 1//2
    h1_ml  = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_taps_1_2.txt"),'\t'))
    y1_ml  = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_y_1_2.txt"),'\t'))
    y1_jl  = resample(x_ml, rate, h1_ml)
    @test y1_jl ≈ y1_ml


    #
    # [y2,b2] = resample(x, 2, 1)
    #
    rate   = 2//1
    h2_ml  = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_taps_2_1.txt"),'\t'))
    y2_ml  = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_y_2_1.txt"),'\t'))
    y2_jl  = resample(x_ml, rate, h2_ml)
    @test y2_jl ≈ y2_ml


    #
    # [y3,b3] = resample(x, 3, 2)
    #
    rate   = 3//2
    h3_ml  = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_taps_3_2.txt"),'\t'))
    y3_ml  = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_y_3_2.txt"),'\t'))
    y3_jl  = resample(x_ml, rate, h3_ml)
    @test y3_jl ≈ y3_ml


    #
    # [y4,b4] = resample(x, 2, 3)
    #
    rate  = 2//3
    h4_ml = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_taps_2_3.txt"),'\t'))
    y4_ml = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_y_2_3.txt"),'\t'))
    y4_jl = resample(x_ml, rate, h4_ml)
    @test y4_jl ≈ y4_ml
end

@testset "array signal" begin
    rate   = 1//2
    x_ml   = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_x.txt"),'\t'))
    h1_ml  = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_taps_1_2.txt"),'\t'))
    y1_ml  = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_y_1_2.txt"),'\t'))

    expected_result = [y1_ml 2y1_ml]
    X = [x_ml 2x_ml]
    y1_jl  = resample(X, rate, h1_ml, dims=1)
    @test y1_jl ≈ expected_result

    y1_jl  = resample(X', rate, h1_ml, dims=2)
    @test y1_jl ≈ expected_result'

    expected_result_3d = permutedims(reshape(expected_result, (size(expected_result, 1), size(expected_result, 2), 1)), (3, 1, 2))
    X_3d = permutedims(reshape(X, (size(X, 1), size(X, 2), 1)), (3, 1, 2))
    y1_jl  = resample(X_3d, rate, h1_ml, dims=2)
    @test y1_jl ≈ expected_result_3d

    expected_result_3d = permutedims(expected_result_3d, (1, 3, 2))
    X_3d = permutedims(X_3d, (1, 3, 2))
    y1_jl  = resample(X_3d, rate, h1_ml, dims=3)
    @test y1_jl ≈ expected_result_3d
end

@testset "irrational ratio" begin
    ratio    = 3.141592653589793
    cycles   = 2
    tx       = range(0, stop=cycles, length=1000)
    x        = sinpi.(2*tx)
    y        = resample(x, ratio)
    yLen     = length(y)
    ty       = range(0, length=yLen, step=step(tx)/ratio)
    yy       = sinpi.(2*ty)
    idxLower = round(Int, yLen/3)
    idxUpper = idxLower*2
    yDelta   = abs.(y[idxLower:idxUpper].-yy[idxLower:idxUpper])
    @test maximum(yDelta) < 0.00025

    # Test Float32 ratio (#302)
    f32_ratio = convert(Float32, ratio)
    f32_y     = resample(x, f32_ratio)
    ty        = range(0, length=yLen, step=step(tx)/f32_ratio)
    yy        = sinpi.(2*ty)
    idxLower  = round(Int, yLen/3)
    idxUpper  = idxLower*2
    yDelta    = abs.(f32_y[idxLower:idxUpper].-yy[idxLower:idxUpper])
    @test maximum(yDelta) < 0.00025
end

@testset "arbitrary ratio" begin
    # https://github.com/JuliaDSP/DSP.jl/issues/317
    @testset "Buffer length calculation" begin
        @test length(resample(sin.(1:1:35546),  1/55.55)) == 641
        @test length(resample(randn(1822), 0.9802414928649835)) == 1787
        @test length(resample(1:16_367_000*2, 10_000_000/16_367_000)) == 20_000_001
        @test resample(zeros(1000), 0.012) == zeros(13)
    end
end

@testset "resample_filter" begin
    @testset "decimation" begin
        ratio = 1//2
        h     = resample_filter(ratio)
        r0    = abs.(freqresp(PolynomialRatio(h, [1]), 0))
        rc    = abs.(freqresp(PolynomialRatio(h, [1]), ratio*π))
        @test isapprox(r0, 1.0)
        @test isapprox(rc, numerator(ratio)/2, rtol=0.001)

        ratio = 1//32
        h     = resample_filter(ratio)
        r0    = abs.(freqresp(PolynomialRatio(h, [1]), 0))
        rc    = abs.(freqresp(PolynomialRatio(h, [1]), ratio*π))
        @test isapprox(r0, numerator(ratio))
        @test isapprox(rc, numerator(ratio)/2, rtol=0.001)
    end

    @testset "interpolation" begin
        ratio = 2//1
        h     = resample_filter(ratio)
        r0    = abs.(freqresp(PolynomialRatio(h, [1]), 0))
        rc    = abs.(freqresp(PolynomialRatio(h, [1]), 1/ratio*π))
        @test isapprox(r0, numerator(ratio))
        @test isapprox(rc, numerator(ratio)/2, rtol=0.001)

        ratio = 32//1
        h     = resample_filter(ratio)
        r0    = abs.(freqresp(PolynomialRatio(h, [1]), 0))
        rc    = abs.(freqresp(PolynomialRatio(h, [1]), 1/ratio*π))
        @test isapprox(r0, numerator(ratio))
        @test isapprox(rc, numerator(ratio)/2, rtol=0.001)
    end

    @testset "arbitrary rate" begin
        ratio = 3.141592653589793
        Nϕ    = 32
        fc    = 1/Nϕ
        h     = resample_filter(ratio, Nϕ)
        r0    = abs.(freqresp(PolynomialRatio(h, [1]), 0))
        rc    = abs.(freqresp(PolynomialRatio(h, [1]), fc*π))
        @test isapprox(r0, Nϕ)
        @test isapprox(rc, Nϕ/2, rtol=0.001)
    end
end

@testset "inputlength" begin
    # FIRDecimator, FIRInterpolator, FIRRational
    for _ in 1:1000
        M=rand(1:10)//rand(1:10)
        H=FIRFilter(zeros(rand(1:100)), M)
        if M != 1
            setphase!(H, 10*rand())
        end
        yL = rand(1:100)
        @test outputlength(H, inputlength(H, yL)) <= yL < outputlength(H, inputlength(H, yL)+1)
        @test outputlength(H, inputlength(H, yL, RoundUp)-1) < yL <= outputlength(H, inputlength(H, yL, RoundUp))
    end

    # FIRArbitrary
    for _ in 1:1000
        M = 10*rand()
        H = FIRFilter(zeros(rand(1:100)), M)
        setphase!(H, 10*rand())
        yL = rand(1:100)
        @test outputlength(H, inputlength(H, yL)) <= yL < outputlength(H, inputlength(H, yL)+1)
        @test outputlength(H, inputlength(H, yL, RoundUp)-1) < yL <= outputlength(H, inputlength(H, yL, RoundUp))
    end
    let
        M = 2.0
        H = FIRFilter(resample_filter(M), M)
        setphase!(H, timedelay(H))
        yL = 200
        @test outputlength(H, inputlength(H, yL)) <= yL < outputlength(H, inputlength(H, yL)+1)
        @test outputlength(H, inputlength(H, yL, RoundUp)-1) < yL <= outputlength(H, inputlength(H, yL, RoundUp))
    end
end
