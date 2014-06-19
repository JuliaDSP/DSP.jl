#in matlab:
#x=rand(512,1);
#[s,f,t,p]=spectrogram(x,ones(1,256),128,256,10);
#save
#
#in julia:
#using MAT
#
#matdata=matread("matlab.mat")
#
#for i in ("x", "f", "t", "p")
#  fid=open("spectrogram_$i.txt","w")
#  print(fid,matdata["$i"])
#  close(fid)
#end

using DSP, Base.Test

x0 = vec(readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_x.txt"),'\t'))
f0 = vec(readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_f.txt"),'\t'))
t0 = vec(readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_t.txt"),'\t'))
p0 = readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_p.txt"),'\t')
spec = spectrogram(x0, 256, 128; fs=10)
p, f, t = power(spec), freq(spec), time(spec)

# with real input matlab outputs a 1-sided PSD
@test_approx_eq p0 p
@test_approx_eq f0 f
@test_approx_eq t0 t

#Matlab: p = pwelch(0:7, [1, 1, 1, 1, 1, 1, 1, 1], 0, 8, 1, 'twosided')
data = Float64[0:7]
data0 = Float64[98.0,
                13.656854249492380,
                 4.0,
                 2.343145750507620,
                 2.0,
                 2.343145750507620,
                 4.0,
                13.656854249492380]
@test_approx_eq power(periodogram(data, onesided=false)) data0
@test_approx_eq power(welch_pgram(data, length(data), 0, onesided=false)) data0

# # ~~~~~~~~ Tests with no window ~~~~~~~~~~~~~~~~~~~
# Matlab: p = pwelch(0:7, [1, 1], 0, 2, 1, 'twosided')
@test_approx_eq power(welch_pgram(data, 2, 0; onesided=false)) Float64[34.5, 0.5]

# Matlab: p = pwelch(0:7, [1, 1, 1], 0, 3, 1, 'twosided')
@test_approx_eq power(welch_pgram(data, 3, 0; onesided=false)) Float64[25.5, 1.0, 1.0]

# Matlab: p = pwelch(0:7, [1, 1, 1], 1, 3, 1, 'twosided')
@test_approx_eq power(welch_pgram(data, 3, 1; onesided=false)) Float64[35.0, 1.0, 1.0]

# Matlab: p = pwelch(0:7, [1, 1, 1, 1], 1, 4, 1, 'twosided')
@test_approx_eq power(welch_pgram(data, 4, 1; onesided=false)) Float64[45, 2, 1, 2]

# ~~~~~~~~~~~ This one tests periodogram ~~~~~~~~~~~~
# ~ If functionality of the other arguments has been
# ~ tested above, we only test here that the correct
# ~ value of the spectral density is obtained when
# ~ using a window. More tests to be added if needed
#Matlab: p = pwelch(0:7, window_func(8), 0, 8, 1, 'twosided')
cases = [
    hamming => Float64[65.461623986801527,
                       20.556791795515764,
                        0.369313143650544,
                        0.022167446610882,
                        0.025502985564107,
                        0.022167446610882,
                        0.369313143650544,
                        20.556791795515764],
    bartlett => Float64[62.999999999999993,
                        21.981076052592442,
                         0.285714285714286,
                         0.161781090264695,
                         0.142857142857143,
                         0.161781090264695,
                         0.285714285714286,
                        21.981076052592442]
]

for (window1, expected) in cases
    @test_approx_eq power(periodogram(data; window=window1, onesided=false)) expected
    @test_approx_eq power(welch_pgram(data, length(data), 0; window=window1, onesided=false)) expected
end
