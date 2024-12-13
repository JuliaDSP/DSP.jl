# diric.jl Dirichlet kernel tests

@testset "diric" begin
    @test_throws DomainError diric(0, -2)
    @test @inferred(diric(0, 4)) ≈ 1
    @test @inferred(diric(0 // 1, 5)) == 1
    @test @inferred(diric(4π, 4)) ≈ 1
    @test @inferred(diric(2π, 4)) ≈ -1
    @test @inferred(diric(0, 5)) ≈ 1
    @test @inferred(diric(2π, 5)) ≈ 1
    @test @inferred(diric(π, 5)) ≈ 1 // 5
    @test @inferred(diric(π/2, 5)) ≈ diric(-π/2, 5) ≈ -0.2
    @test abs(@inferred(diric(π/2, 4))) < eps()
    @test abs(@inferred(diric(3π/2, 4))) < eps()

    # check type inference
    @inferred diric(2, 4)
    @inferred diric(Float16(1), 4)
    @inferred diric(0.1f0, 4)
    @inferred diric(0., 4)
    @inferred diric(BigFloat(1), 4)

    # accuracy test near 2π for even n
    dirics(Ω,n) = real(exp(-1im * (n-1) * Ω/2)/n * sum(exp.(1im*(0:(n-1))*Ω)))
    Ω = 2π .+ 4π * (-101:2:101)/100 * 1e-9
    n = 400
    @test dirics.(Ω,n) ≈ diric.(Ω,n)
end
