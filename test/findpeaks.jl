using DSP
using Test

@testset "Findpeaks: Findpeaks: Single Gaussian peak" begin
    gaussian(x,μ,σ) = 1/sqrt(2*π)/σ*exp(-((x-μ)^2)/2/σ^2)
    x = collect(range(1,stop=1000,length=1000))
    μ = rand(x)
    data = gaussian.(x,μ,rand(x)/3)

    p = findpeaks(data)
    @test length(p) == 1
    @test p[1] == μ

    p = findpeaks(data,x)
    @test length(p) == 1
    @test p[1] == μ
end


@testset "Findpeaks: Prominence" begin
    gaussian(x,μ,σ) = exp(-((x-μ)^2)/2/σ^2)

    x = range(0., stop = 1., length = 1000) |> collect

    positions = [1, 2, 3, 5, 7, 8]/10;
    amps = [3, 7, 5, 5, 4, 5];
    widths = [1, 3, 3, 4, 2, 3]/100;

    base = 4 * cos.(2π * x);
    signal = ( a * gaussian.(x, μ, σ) for (a, μ, σ) in zip(amps, positions, widths) )
    data = sum(signal) + base

    peaks = findpeaks(data, x, min_prom=4.)
    
    @test length(peaks) == 2

    expected_peak_1 = argmin(abs.(x .- positions[2])) # global highest peak
    expected_peak_2 = argmin(abs.(x .- positions[4])) # lowest peak

    # check peaks are around expected - may be shifted because of background
    @test abs(peaks[1] - expected_peak_1) < 20
    @test abs(peaks[2] - expected_peak_2) < 20
end

@testset "Findpeaks: Threshold" begin
    y = [0., 1., 3., 1., 4., 5., 3., 0., 2., 5., 4., 0.]
    
    @test Set(findpeaks(y, threshold=3.1)) == Set([])
    @test Set(findpeaks(y, threshold=2.1)) == Set([])
    @test Set(findpeaks(y, threshold=1.1)) == Set([3])
    @test Set(findpeaks(y, threshold=0.1)) == Set([3, 6, 10])
    @test Set(findpeaks(y, threshold=0.0)) == Set([3, 6, 10])
end

@testset "Findpeaks: Min. Distance" begin
    y = [0., 1., 3., 1., 4., 5., 3., 0., 2., 5., 4., 0.]
    
    @test Set(findpeaks(y, min_dist=4)) == Set([6])
    @test Set(findpeaks(y, min_dist=3)) == Set([6, 10])
    @test Set(findpeaks(y, min_dist=2)) == Set([3, 6, 10])
    @test Set(findpeaks(y, min_dist=0)) == Set([3, 6, 10])
    @test Set(findpeaks(y, min_dist=0)) == Set([3, 6, 10])
end

@testset "Findpeaks: Min. Height" begin
    y = [0., 1., 3., 1., 4., 5., 3., 0., 2., 5., 4., 0.]
    
    @test Set(findpeaks(y, min_height=6.)) == Set([])
    @test Set(findpeaks(y, min_height=4.9)) == Set([6, 10])
    @test Set(findpeaks(y, min_height=2.9)) == Set([3, 6, 10])
end



