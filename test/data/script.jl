# scratch work, delete me before marking PR as ready to review

# fs = 1000;
# t = 0:1/fs:2-1/fs;
# x = cos(2*pi*100*t)+randn(size(t));
# nfft = 2^nextpow2(length(x));
# [pxx,fx] = pmtm(x,4,nfft,fs);
# y = sin(2*pi*100*t)+randn(size(t));
# z = x + 1i*y;
# [pzz,fz] = pmtm(z,4,nfft,fs);

fs = 1000;
nw=4
nfft = nextpow(2, length(x))

t = 0:1/fs:2-1/fs;
x = cos.((2*pi*100) .* t) .+ randn(length(t))
y = sin.((2*pi*100) .* t) .+ randn(length(t))
nfft = nextpow(2, length(x))

mt_pgram(x; fs, nw, nfft)
z = x + im*y
mt_pgram(z; fs, nw, nfft)


using PyMNE
fs = 1000.0
n_samples = 1024
t = (0:1023) ./ fs
sin_1 = sin.(2 * π * 12.0 * t)  # 12 Hz sinusoid signal
sin_2 = sin.(2 * π * 12.0 * t .+ π)
noise = rand(1024) * 2 .- 1

data = Array{Float64, 3}(undef, 1, 2, n_samples)
data[1,1,:] = sin_1
data[1,2,:] = sin_2

result = PyMNE.time_frequency.csd_array_multitaper(data, fs, n_fft = nextpow(2, n_samples), low_bias=true, adaptive=false)
out_f = result.frequencies
out= reduce((x,y) -> cat(x,y; dims=3), [ r.get_data() for r in result])

signal = dropdims(data; dims=1)
jl_result = mt_cross_spectral(signal; fs, n_tapers = 8)
@test freq(jl_result)[2:end] ≈ out_f
@assert !any(isnan, jl_result.values)
findmax(abs.(jl_result.values[:,:,2:end]./fs .- out))
@test jl_result.values[:,:,2:end] ./ fs ≈ out atol=1e-3
jl = [jl_result.values[:,:,i] for i in axes(jl_result.values,3)]
findmax(norm.(jl[2:end] .- [ r.get_data() for r in result]))

using MAT, DelimitedFiles
mat_results = matread("test/data/JuliaMatlabCompare.mat")

writedlm("pmtm_x.txt", vec(mat_results["x"]))
writedlm("pmtm_y.txt", vec(mat_results["y"]))
writedlm("pmtm_fx.txt", vec(mat_results["fx"]))
writedlm("pmtm_fz.txt", vec(mat_results["fz"]))
writedlm("pmtm_pxx.txt", vec(mat_results["pxx"]))
writedlm("pmtm_pzz.txt", vec(mat_results["pzz"]))
