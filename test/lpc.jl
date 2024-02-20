# Filter some noise, try to learn the coefficients
using Statistics: std
@testset "$method" for method in (LPCBurg(), LPCLevinson())
    @testset "$T" for T in (Float64, ComplexF64)
        coeffs = T <: Complex ? [1, 0.5 + 0.2im, 0.3 - 0.5im, 0.2] : [1, 0.5, 0.3, 0.2]
        s = randn(T, 50000)
        x = filt(1, coeffs, s)

        # Analyze the filtered noise, attempt to reconstruct the coefficients above
        ar, e = lpc(x[1000:end], length(coeffs) - 1, method)

        @test all(<(0.01), abs.(ar .- coeffs[2:end]))
        @test isapprox(e, std(s); rtol=0.01)
    end
end

# test dotu, levinson with complex Int coefficients
@test isapprox(levinson(complex(1:10), 3)[1] |> real, -[1.25, 0, 0.25])

# test that lpc defaults to Burg
@test let v = rand(1000)
    lpc(v, 20) == lpc(v, 20, LPCBurg())
end
