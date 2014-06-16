using DSP, Base.Test

# Test dpss against dpss computed with MATLAB
d1 = dpss(128, 4)
d2 = readdlm(joinpath(dirname(@__FILE__), "data", "dpss128,4.txt"), '\t')
@test_approx_eq d1 d2


# make sure return types are correct
n = 12
ft = typeof(1.0)
for W in(:rect, :hanning, :hamming, :cosine, :lanczos, :triang, :bartlett, :bartlett_hann, :blackman)
    @eval @test Array{ft,1} == typeof($W(n)) && length($W(n)) == n
end

@test Array{ft,1} == typeof(gaussian(n, 0.4)) && length(gaussian(n, 0.4)) == n
@test Array{ft,1} == typeof(kaiser(n, 0.4)) && length(kaiser(n, 0.4)) == n
@test Array{ft,2} == typeof(dpss(n, 1.5)) && size(dpss(n, 1.5),1) == n  # size(,2) depends on the parameters

