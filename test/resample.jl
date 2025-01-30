using DSP
using Test
using DelimitedFiles: readdlm

reference_data(s) = vec(readdlm(joinpath(dirname(@__FILE__), "data", s), '\t'))

@testset "rational ratio" begin
    # AM Modulator
    # sig(t) = [(1 + sin(2π*0.005*t)) * sin(2π*.05*t) for t in t]
    x_ml  = reference_data("resample_x.txt")

    #
    # [y,b] = resample(x, $num, $den)
    #
    for rate in (1//2, 2//1, 3//2, 2//3)
        num, den = numerator(rate), denominator(rate)
        h1_ml = reference_data("resample_taps_$(num)_$(den).txt")
        y1_ml = reference_data("resample_y_$(num)_$(den).txt")
        y1_jl = resample(x_ml, rate, h1_ml)
        @test y1_ml ≈ y1_jl
        @test y1_ml ≈ resample(x_ml, rate) rtol = 0.001     # check default taps are ok
    end
end

@testset "array signal" begin
    x, mat, rate = rand(10_000), rand(121, 212), 1.23456789
    h = resample_filter(rate)
    # check that filtering with taps from `resample_filter` give equal results
    @test resample(x, rate) == resample(x, rate, h)
    @test resample(x, rate, 64, 0.8) == resample(x, rate, resample_filter(rate, 64, 0.8), 64)
    @test resample(mat, rate; dims=1) == resample(mat, rate, h; dims=1)
    @test resample(mat, rate; dims=2) == resample(mat, rate, h; dims=2)

    rate  = 1//2
    x_ml  = reference_data("resample_x.txt")
    h1_ml = reference_data("resample_taps_1_2.txt")
    y1_ml = reference_data("resample_y_1_2.txt")
    expected_result = [y1_ml ℯ * y1_ml]
    X = [x_ml ℯ * x_ml]

    y1_jl  = resample(X, rate, h1_ml; dims=1)
    arb_y1 = resample(X, float(rate); dims=1)
    rat_y1 = resample(X, rate; dims=1)
    @test y1_jl ≈ expected_result
    @test arb_y1 ≈ expected_result  rtol = 0.002     # test that default taps are good enough
    @test rat_y1 ≈ expected_result  rtol = 0.0005

    y1_jl  = resample(X', rate, h1_ml; dims=2)
    arb_y1 = resample(X', float(rate); dims=2)
    rat_y1 = resample(X', rate; dims=2)
    @test y1_jl ≈ expected_result'
    @test arb_y1 ≈ expected_result' rtol = 0.002
    @test rat_y1 ≈ expected_result' rtol = 0.0005

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
        @test length(resample(sin.(1:1:35546),  1/55.55)) == 640
        @test length(resample(randn(1822), 0.9802414928649835)) == 1786
        @test length(resample(1:16_367_000*2, 10_000_000/16_367_000)) == 20_000_000
        @test resample(zeros(1000), 0.012) == zeros(12)
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

@testset "FIRFilter types" begin
    using DSP.Filters: FIRStandard, FIRInterpolator
    # test FIRRational(::Vector, ::Int(/Int32)) inferred result type
    FIRFutyp = Union{FIRFilter{<:FIRInterpolator},FIRFilter{<:FIRStandard}}
    @test only(Base.return_types(FIRFilter, (Vector, Int64))) <: FIRFutyp broken=VERSION<v"1.11"
    @test only(Base.return_types(FIRFilter, (Vector, Int32))) <: FIRFutyp broken=VERSION<v"1.11"
    FIRFutype{T} = Union{FIRFilter{FIRInterpolator{T}},FIRFilter{FIRStandard{T}}} where T
    @test only(Base.return_types(FIRFilter, (Vector{Float64}, Int64))) == FIRFutype{Float64}
    @test only(Base.return_types(FIRFilter, (Vector{Float64}, Int32))) == FIRFutype{Float64}

    # check that non-Int / Rational{Int} ratios get converted properly
    x = rand(200)
    @test resample(x, Int32(3)) == resample(x, Int64(3))
    @test resample(x, Rational{Int32}(1, 3)) == resample(x, Rational{Int64}(1, 3))
    @test resample(x, Rational{Int32}(2, 3)) == resample(x, Rational{Int64}(2, 3))
end
