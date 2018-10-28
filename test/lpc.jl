# Filter some noise, try to learn the coefficients
@testset "$method" for method in (LPCBurg(), LPCLevinson())
    coeffs = [1, .5, .3, .2]
    x = filt([1], coeffs, randn(50000))

    # Analyze the filtered noise, attempt to reconstruct the coefficients above
    ar, e = lpc(x[1000:end], length(coeffs)-1, method)

    @test all(abs.([1; ar] .- coeffs) .<= 1e-2)
end
