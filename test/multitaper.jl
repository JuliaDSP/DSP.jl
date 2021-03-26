
@testset "`MTConfig`" begin
    @test_throws ArgumentError MTConfig{Float64}(100, 1; n_for_fft = 99)
    @test_throws ArgumentError MTConfig{Float64}(100, -1)
    @test_throws ArgumentError MTConfig{Float64}(100, 1; n_tapers = -1)
    @test_throws DimensionMismatch MTConfig{Float64}(100, 1; window = rand(1, 200))
    @test_throws ArgumentError MTConfig{Complex{Float64}}(100, 1; onesided=true)
    
end
