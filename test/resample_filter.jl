using DSP
using MAT
using Base.Test
using PyPlot

# AM Modulator
sig(t) = [(1 + sin(2π*0.005*t)) * sin(2π*.05*t) for t in t]

mlData = matread(joinpath(dirname(@__FILE__), "data","results.mat"))
tx_jl  = [0:100]
x_jl   = sig(tx_jl) # Cretae signal vector
x_ml   = mlData["x"]' 
@test_approx_eq x_jl x_ml


#
# [y1,b1] = resample(x, 1, 2)
#
y1_ml = vec(mlData["y1"])
h1_ml = vec(mlData["b1"])
ty1_ml = [0:length(y1_ml)-1] / float64(1//2)
# y1_jl = resample(x_jl, 1//2)
# ty1_jl = [0:length(y1_jl)-1] / float64(1//2)
# @test_approx_eq y1_jl y1_ml

figure(1)
suptitle("ratio = 1/2")
hold(true)
plot(tx_jl, x_jl, "k.-")
stem(ty1_ml, y1_ml, "b")
hold(false)


#
# [y2,b2] = resample(x, 2, 1)
#
y2_ml  = vec(mlData["y2"])
h2_ml  = vec(mlData["b2"])
ty2_ml = [0:length(y2_ml)-1] / float64(2//1)
y2_jl  = resample(x_jl, 2//1)
ty2_jl  = [0:length(y2_jl)-1] / float64(2//1)
# @test_approx_eq y2_jl y2_ml

figure(2)
suptitle("ratio = 2/1")
hold(true)
plot(tx_jl, x_jl, "k.-")
stem(ty2_ml, y2_ml, "b")
stem(ty2_jl, y2_jl, "r")
hold(false)


#
# [y3,b3] = resample(x, 3, 2)
#
h3_ml  = vec(mlData["b3"])

y3_ml  = vec(mlData["y3"])
ty3_ml = [0:length(y3_ml)-1] / float64(3//2)

y3_jl  = resample(x_jl, 3//2, h3_ml)
ty3_jl = [0:length(y3_jl)-1] / float64(3//2)

@test_approx_eq y3_jl y3_ml

figure(3)
suptitle("ratio = 3/2")
hold(true)
plot(tx_jl, x_jl, "k.-")
stem(ty3_ml, y3_ml, "b.-")
stem(ty3_jl, y3_jl, "r.-")
hold(false)


#
# [y4,b4] = resample(x, 2, 3)
#
h4_ml  = vec(mlData["b4"])

y4_ml  = vec(mlData["y4"])
ty4_ml = [0:length(y4_ml)-1] / float64(2//3)

y4_jl  = resample(x_jl, 2//3, h4_ml)
ty4_jl = [0:length(y4_jl)-1] / float64(2//3)

@test_approx_eq y4_jl y4_ml

figure(4)
suptitle("ratio = 2/3")
hold(true)
plot(tx_jl, x_jl, "k.-")
stem(ty4_ml, y4_ml, "b.-")
stem(ty4_jl, y4_jl, "r.-")
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
