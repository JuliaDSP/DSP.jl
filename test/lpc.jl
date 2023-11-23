# Filter some noise, try to learn the coefficients
@testset "$method" for method in (LPCBurg(), LPCLevinson())
    @testset "$T" for T in (Float64, ComplexF64)
        coeffs = T <: Complex ? [1, 0.5 + 0.2im, 0.3 - 0.5im, 0.2] : [1, .5, .3, .2]
        x = filt(1, coeffs, randn(T, 50000))

        # Analyze the filtered noise, attempt to reconstruct the coefficients above
        ar, e = lpc(x[1000:end], length(coeffs)-1, method)

        @test all(<(0.01), abs.(ar - coeffs[2:end]))
    end
end
