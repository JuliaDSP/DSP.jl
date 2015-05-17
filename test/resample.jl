using DSP
using MAT
using Base.Test

# AM Modulator
sig(t) = [(1 + sin(2π*0.005*t)) * sin(2π*.05*t) for t in t]

mlData = matread(joinpath(dirname(@__FILE__), "data","resample.mat"))
tx_jl  = [0:100]
ty_jl  = [0:2:100]
x_jl   = sig(tx_jl) # Cretae signal vector
x_ml   = vec(mlData["x"])
@test_approx_eq x_jl x_ml


#
# [y1,b1] = resample(x, 1, 2)
#
rate   = 1//2
h1_ml  = vec(mlData["b1"])
y1_ml  = vec(mlData["y1"])
ty1_ml = [0:length(y1_ml)-1] / float64(rate)
y1_jl  = resample(x_ml, rate, h1_ml)
ty1_jl = [0:length(y1_jl)-1] / float64(rate)
# h1_jl  = resample_filter(rate)
@test_approx_eq y1_jl y1_ml


#
# [y2,b2] = resample(x, 2, 1)
#
rate   = 2//1
h2_ml  = vec(mlData["b2"])
y2_ml  = vec(mlData["y2"])
ty2_ml = [0:length(y2_ml)-1] / float64(rate)
y2_jl  = resample(x_jl, rate, h2_ml)
ty2_jl = [0:length(y2_jl)-1] / float64(rate)
# h2_jl  = resample_filter(rate)
@test_approx_eq y2_jl y2_ml


#
# [y3,b3] = resample(x, 3, 2)
#
rate   = 3//2
h3_ml  = vec(mlData["b3"])
y3_ml  = vec(mlData["y3"])
ty3_ml = [0:length(y3_ml)-1] / float64(rate)
y3_jl  = resample(x_jl, rate, h3_ml)
ty3_jl = [0:length(y3_jl)-1] / float64(rate)
# h3_jl  = resample_filter(rate)
@test_approx_eq y3_jl y3_ml


#
# [y4,b4] = resample(x, 2, 3)
#
rate   = 2//3
h4_ml  = vec(mlData["b4"])
y4_ml  = vec(mlData["y4"])
ty4_ml = [0:length(y4_ml)-1] / float64(rate)
y4_jl  = resample(x_jl, rate, h4_ml)
ty4_jl = [0:length(y4_jl)-1] / float64(rate)
# h4_jl  = resample_filter(rate)
@test_approx_eq y4_jl y4_ml
