using DSP, Compat, Compat.Test, Compat.DelimitedFiles

@testset "dspss" begin
    # Test dpss against dpss computed with MATLAB
    d1 = dpss(128, 4)
    d2 = readdlm(joinpath(dirname(@__FILE__), "data", "dpss128,4.txt"), '\t')
    @test d1 ≈ d2

    # Test dpsseig against dpss from MATLAB
    lambda = [0.9999999997159923,0.9999999731146645,0.9999988168667646,
              0.9999680890685374,0.9994167543397652,0.9925560207018469,
              0.9368556668429153]
    @test dpsseig(d1, 4) ≈ lambda
end

@testset "common windows" begin
    # Checking Hanning, Hamming, Triangular, Bartlett, Bartlett-Hann, Tukey,
    # and Blackman windows against values computed with MATLAB.
    # Lanczos and cosine are checked against values generated with DSP.jl v0.4.0
    # to test for regressions, as there's no reference MATLAB implementation
    # Gaussian is compared against DSP.jl from v0.5.0, when the implementation
    # was corrected (see GH issue #204)
    @test rect(128) == ones(128)

    hanning_jl = hanning(128)
    hanning_ml = readdlm(joinpath(dirname(@__FILE__), "data", "hanning128.txt"), '\t')
    @test hanning_jl ≈ hanning_ml

    hamming_jl = hamming(128)
    hamming_ml = readdlm(joinpath(dirname(@__FILE__), "data", "hamming128.txt"), '\t')
    @test hamming_jl ≈ hamming_ml

    triang_jl = triang(128)
    triang_ml = readdlm(joinpath(dirname(@__FILE__), "data", "triang128.txt"), '\t')
    @test triang_jl ≈ triang_ml

    bartlett_jl = bartlett(128)
    bartlett_ml = readdlm(joinpath(dirname(@__FILE__), "data", "bartlett128.txt"), '\t')
    @test bartlett_jl ≈ bartlett_ml

    barthann_jl = bartlett_hann(128)
    barthann_ml = readdlm(joinpath(dirname(@__FILE__), "data", "bartlett_hann128.txt"), '\t')
    @test bartlett_jl ≈ bartlett_ml

    blackman_jl = blackman(128)
    blackman_ml = readdlm(joinpath(dirname(@__FILE__), "data", "blackman128.txt"), '\t')
    @test blackman_jl ≈ blackman_ml

    kaiser_jl = kaiser(128, 0.4/π)
    kaiser_ml = readdlm(joinpath(dirname(@__FILE__), "data", "kaiser128,0.4.txt"), '\t')
    @test kaiser_jl ≈ kaiser_ml

    gaussian_jl = gaussian(128, 0.2)
    gaussian_ref = readdlm(joinpath(dirname(@__FILE__), "data", "gaussian128,0.2.txt"), '\t')
    @test gaussian_jl ≈ gaussian_ref

    tukey_jl = tukey(128, 0.4)
    tukey_ml = readdlm(joinpath(dirname(@__FILE__), "data", "tukey128,0.4.txt"), '\t')
    @test tukey_jl ≈ tukey_ml

    lanczos_jl = lanczos(128)
    lanczos_ref = readdlm(joinpath(dirname(@__FILE__), "data", "lanczos128.txt"), '\t')
    @test lanczos_jl ≈ lanczos_ref

    cosine_jl = cosine(128)
    cosine_ref = readdlm(joinpath(dirname(@__FILE__), "data", "cosine128.txt"), '\t')
    @test cosine_jl ≈ cosine_ref
end

@testset "window return types" begin
    # make sure return types are correct
    n = 12
    ft = typeof(1.0)
    for W in(:rect, :hanning, :hamming, :cosine, :lanczos, :triang, :bartlett, :bartlett_hann, :blackman)
        @eval @test Array{$ft,1} == typeof($W($n)) && length($W($n)) == $n
    end

    @test Array{ft,1} == typeof(tukey(n, 0.4)) && length(tukey(n, 0.4)) == n
    @test Array{ft,1} == typeof(gaussian(n, 0.4)) && length(gaussian(n, 0.4)) == n
    @test Array{ft,1} == typeof(kaiser(n, 0.4)) && length(kaiser(n, 0.4)) == n
    @test Array{ft,2} == typeof(dpss(n, 1.5)) && size(dpss(n, 1.5),1) == n  # size(,2) depends on the parameters
end

@testset "tensor product windows" begin
    # tensor product windows
    w = hamming(15)
    w2 = hamming(20)
    @test w*w2' ≈ hamming((15,20))
    w = tukey(10, 0.4)
    w2 = tukey(4, 0.4)
    @test w*w2' ≈ tukey((10,4), 0.4)
end
