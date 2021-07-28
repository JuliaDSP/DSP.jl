# diric.jl Dirichlet kernel tests

@testset "diric" begin
    @test_throws ArgumentError diric(0, -2)
    @test @inferred(diric(0, 4)) ≈ 1
    @test @inferred(diric(0 // 1, 5)) == 1
    @test @inferred(diric(4π, 4)) ≈ 1
    @test @inferred(diric(2π, 4)) ≈ -1
    @test @inferred(diric(0, 5)) ≈ 1
    @test @inferred(diric(2π, 5)) ≈ 1
    @test @inferred(diric(π, 5)) == 1 // 5

    # check type inference
    @inferred diric(2, 4)
    @inferred diric(Float16(1), 4)
    @inferred diric(0.1f0, 4)
    @inferred diric(0., 4)
    @inferred diric(BigFloat(1), 4)
end
