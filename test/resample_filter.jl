using DSP
using MAT
using Base.Test

function Filters.freqz(h::Vector, w = linspace(0, π, 250))
    pr = PolynomialRatio(h, [one(eltype(h))])
    freqz(pr, w)
end

Util.amp2db(x::Complex) = amp2db(abs(x))
Util.amp2db(x::AbstractVector) = [amp2db(x) for x in x]

# AM Modulator
sig(t) = [(1 + sin(2π*0.005*t)) * sin(2π*.05*t) for t in t]


mlData = matread(joinpath(dirname(@__FILE__), "data","resample.mat"))
tx_jl  = [0:100]
x_jl   = sig(tx_jl) # Cretae signal vector
x_ml   = vec(mlData["x"])
@test_approx_eq x_jl x_ml


#
# [y1,b1] = resample(x, 1, 2)
#
rate   = 1//2
h1_ml  = vec(mlData["b1"])
h1_jl  = resample_filter(rate)
y1_ml  = vec(mlData["y1"])
ty1_ml = [0:length(y1_ml)-1] / float64(rate)
y1_jl  = resample(x_ml, rate, h1_jl)
ty1_jl = [0:length(y1_jl)-1] / float64(rate)
# @test_approx_eq y1_jl y1_ml


#
# [y2,b2] = resample(x, 2, 1)
#
rate   = 2//1
h2_ml  = vec(mlData["b2"])
h2_jl  = resample_filter(rate)
y2_ml  = vec(mlData["y2"])
ty2_ml = [0:length(y2_ml)-1] / float64(rate)
y2_jl  = resample(x_jl, rate, h2_ml)
ty2_jl = [0:length(y2_jl)-1] / float64(rate)
@test_approx_eq y2_jl y2_ml


#
# [y3,b3] = resample(x, 3, 2)
#
rate   = 3//2
h3_ml  = vec(mlData["b3"])
h3_jl  = resample_filter(rate)
y3_ml  = vec(mlData["y3"])
ty3_ml = [0:length(y3_ml)-1] / float64(rate)
y3_jl  = resample(x_jl, rate, h3_ml)
ty3_jl = [0:length(y3_jl)-1] / float64(rate)
@test_approx_eq y3_jl y3_ml


#
# [y4,b4] = resample(x, 2, 3)
#
rate   = 2//3
h4_ml  = vec(mlData["b4"])
h4_jl  = resample_filter(rate)
y4_ml  = vec(mlData["y4"])
ty4_ml = [0:length(y4_ml)-1] / float64(rate)
y4_jl  = resample(x_jl, rate, h4_ml)
ty4_jl = [0:length(y4_jl)-1] / float64(rate)
@test_approx_eq y4_jl y4_ml



using PyPlot


ml_lf = "b-" # Line format for matlab results
ml_mf = "bo" # Marker format for matlab results
jl_lf = "r-" # Line format for julia results
jl_mf = "r." # Marker format for julia results


figure(1)
clf()
subplot(4,2,1)
title("ratio = 1//2")
hold(true)
plot(tx_jl, x_jl, "k.-")
stem(ty1_ml, y1_ml, markerfmt = ml_mf, linefmt = ml_lf)
stem(ty1_jl, y1_jl, markerfmt = jl_mf, linefmt = jl_lf)
hold(false)

subplot(4,2,3)
title("Error for ratio = 1//2")
hold(true)
plot(ty1_jl, y1_jl.-y1_ml)
hold(false)


subplot(4,2,2)
title("ratio = 2//1")
hold(true)
plot(tx_jl, x_jl, "k.-")
stem(ty2_ml, y2_ml, markerfmt = ml_mf, linefmt = ml_lf)
stem(ty2_jl, y2_jl, markerfmt = jl_mf, linefmt = jl_lf)
hold(false)

subplot(4,2,4)
title("Error for ratio = 2//1")
hold(true)
plot(ty2_jl, y2_jl.-y2_ml)
hold(false)


subplot(4,2,5)
title("ratio = 3//2")
hold(true)
plot(tx_jl, x_jl, "k.-")
stem(ty3_ml, y3_ml, markerfmt = ml_mf, linefmt = ml_lf)
stem(ty3_jl, y3_jl, markerfmt = jl_mf, linefmt = jl_lf)
hold(false)

subplot(4,2,7)
title("Error for ratio = 3//2")
hold(true)
plot(ty3_jl, y3_jl.-y3_ml)
hold(false)


subplot(4,2,6)
title("ratio = 2//3")
hold(true)
plot(tx_jl, x_jl, "k.-")
stem(ty4_ml, y4_ml, markerfmt = ml_mf, linefmt = ml_lf)
stem(ty4_jl, y4_jl, markerfmt = jl_mf, linefmt = jl_lf)
hold(false)

subplot(4,2,8)
title("Error for ratio = 2//3")
hold(true)
plot(ty4_jl, y4_jl.-y4_ml)
hold(false)



figure(2)
clf()
subplot(2,2,1)
title("ratio = 1/2")
hold(true)
plot(amp2db(freqz(h1_ml)), ml_lf)
plot(amp2db(freqz(h1_jl)), jl_lf)
hold(false)

subplot(2,2,2)
title("ratio = 2/1")
hold(true)
plot(amp2db(freqz(h2_ml)), ml_lf)
plot(amp2db(freqz(h2_jl)), jl_lf)
hold(false)

subplot(2,2,3)
title("ratio = 3/2")
hold(true)
plot(amp2db(freqz(h3_ml)), ml_lf)
plot(amp2db(freqz(h3_jl)), jl_lf)
hold(false)

subplot(2,2,4)
title("ratio = 2/3")
hold(true)
plot(amp2db(freqz(h4_ml)), ml_lf)
plot(amp2db(freqz(h4_jl)), jl_lf)
hold(false)


# #
# # Resample with arbitrary factor
# #
# y5_jl  = resample(x_jl, 3/2, h3_ml)
# ty5_jl = [0.0:length(y5_jl)-1]/float64(3//2)
#
# figure(5)
# suptitle("ratio = $(3/2)")
# hold(true)
# plot(tx_jl, x_jl, "k.-")
# stem(ty5_jl, y5_jl, "r.-")
# hold(false)
