using DSP
using Base.Test

# AM Modulator
# sig(t) = [(1 + sin(2π*0.005*t)) * sin(2π*.05*t) for t in t]
x_ml   = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_x.txt"),'\t'))

#
# [y1,b1] = resample(x, 1, 2)
#
rate   = 1//2
h1_ml  = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_taps_1_2.txt"),'\t'))
y1_ml  = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_y_1_2.txt"),'\t'))
y1_jl  = resample(x_ml, rate, h1_ml)
@test_approx_eq y1_jl y1_ml


#
# [y2,b2] = resample(x, 2, 1)
#
rate   = 2//1
h2_ml  = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_taps_2_1.txt"),'\t'))
y2_ml  = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_y_2_1.txt"),'\t'))
y2_jl  = resample(x_ml, rate, h2_ml)
@test_approx_eq y2_jl y2_ml


#
# [y3,b3] = resample(x, 3, 2)
#
rate   = 3//2
h3_ml  = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_taps_3_2.txt"),'\t'))
y3_ml  = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_y_3_2.txt"),'\t'))
y3_jl  = resample(x_ml, rate, h3_ml)
@test_approx_eq y3_jl y3_ml


#
# [y4,b4] = resample(x, 2, 3)
#
rate  = 2//3
h4_ml = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_taps_2_3.txt"),'\t'))
y4_ml = vec(readdlm(joinpath(dirname(@__FILE__), "data", "resample_y_2_3.txt"),'\t'))
y4_jl = resample(x_ml, rate, h4_ml)
@test_approx_eq y4_jl y4_ml

#
# Test resample with a rational factor
#

ratio    = 3.141592653589793
cycles   = 2
tx       = linspace(0, cycles, 1000)
x        = sinpi(2*tx)
y        = resample(x, ratio)
yLen     = length(y)
ty       = linspace(0, cycles, yLen)
yy       = sinpi(2*ty)
idxLower = round(Int, yLen/3)
idxUpper = idxLower*2
yDelta   = abs(y[idxLower:idxUpper].-yy[idxLower:idxUpper])
@test all(map(delta -> abs(delta) < 0.005, yDelta))


#
# Test response of resample_filter
#

# Decimation
ratio = 1//2
h     = resample_filter(ratio)
r0    = abs(freqz(PolynomialRatio(h, [1]), 0))
rc    = abs(freqz(PolynomialRatio(h, [1]), ratio*π))
@test isapprox(r0, 1.0)
@test isapprox(rc, num(ratio)/2, rtol=0.001)

ratio = 1//32
h     = resample_filter(ratio)
r0    = abs(freqz(PolynomialRatio(h, [1]), 0))
rc    = abs(freqz(PolynomialRatio(h, [1]), ratio*π))
@test isapprox(r0, num(ratio))
@test isapprox(rc, num(ratio)/2, rtol=0.001)


# Interpolation
ratio = 2//1
h     = resample_filter(ratio)
r0    = abs(freqz(PolynomialRatio(h, [1]), 0))
rc    = abs(freqz(PolynomialRatio(h, [1]), 1/ratio*π))
@test isapprox(r0, num(ratio))
@test isapprox(rc, num(ratio)/2, rtol=0.001)

ratio = 32//1
h     = resample_filter(ratio)
r0    = abs(freqz(PolynomialRatio(h, [1]), 0))
rc    = abs(freqz(PolynomialRatio(h, [1]), 1/ratio*π))
@test isapprox(r0, num(ratio))
@test isapprox(rc, num(ratio)/2, rtol=0.001)

# Arbitrary rate resampling
ratio = 3.141592653589793
Nϕ    = 32
fc    = 1/Nϕ
h     = resample_filter(ratio, Nϕ)
r0    = abs(freqz(PolynomialRatio(h, [1]), 0))
rc    = abs(freqz(PolynomialRatio(h, [1]), fc*π))
@test isapprox(r0, Nϕ)
@test isapprox(rc, Nϕ/2, rtol=0.001)
