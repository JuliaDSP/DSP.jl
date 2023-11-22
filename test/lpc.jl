# Filter some noise, try to learn the coefficients
@testset "$method, Float64" for method in (LPCBurg(), LPCLevinson())
    coeffs = [1, .5, .3, .2]
    x = filt(1, coeffs, randn(50000))

        # Analyze the filtered noise, attempt to reconstruct the coefficients above
    ar, e = lpc(x[1000:end], length(coeffs)-1, method)

    @test isapprox(ar, coeffs[2:end]; atol=0.01)
end

@testset "$method, ComplexF64" for method in (LPCBurg(), LPCLevinson())
    coeffs = [1, 0.5 + 0.2im, 0.3 - 0.5im, 0.2]
    x = filt(1, coeffs, randn(ComplexF64, 50000))

    # Analyze the filtered noise, attempt to reconstruct the coefficients above
    ar, e = lpc(x[1000:end], length(coeffs) - 1, method)

    @test isapprox(ar, coeffs[2:end]; atol=0.01)
end
