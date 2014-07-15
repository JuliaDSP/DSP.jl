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
data = 0:7
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
@test_approx_eq power(spectrogram(data, length(data), 0, onesided=false)) data0

# # ~~~~~~~~ Tests with no window ~~~~~~~~~~~~~~~~~~~
# Matlab: p = pwelch(0:7, [1, 1], 0, 2, 1, 'twosided')
expected = Float64[34.5, 0.5]
@test_approx_eq power(welch_pgram(data, 2, 0; onesided=false)) expected
@test_approx_eq mean(power(spectrogram(data, 2, 0; onesided=false)), 2) expected

# Matlab: p = pwelch(0:7, [1, 1, 1], 0, 3, 1, 'twosided')
expected = Float64[25.5, 1.0, 1.0]
@test_approx_eq power(welch_pgram(data, 3, 0; onesided=false)) expected
@test_approx_eq mean(power(spectrogram(data, 3, 0; onesided=false)), 2) expected

# Matlab: p = pwelch(0:7, [1, 1, 1], 1, 3, 1, 'twosided')
expected = Float64[35.0, 1.0, 1.0]
@test_approx_eq power(welch_pgram(data, 3, 1; onesided=false)) expected
@test_approx_eq mean(power(spectrogram(data, 3, 1; onesided=false)), 2) expected

# Matlab: p = pwelch(0:7, [1, 1, 1, 1], 1, 4, 1, 'twosided')
expected = Float64[45, 2, 1, 2]
@test_approx_eq power(welch_pgram(data, 4, 1; onesided=false)) expected
@test_approx_eq mean(power(spectrogram(data, 4, 1; onesided=false)), 2) expected

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
    @test_approx_eq power(spectrogram(data, length(data), 0; window=window1, onesided=false)) expected
    @test_approx_eq power(periodogram(data; window=window1(length(data)), onesided=false)) expected
    @test_approx_eq power(welch_pgram(data, length(data), 0; window=window1(length(data)), onesided=false)) expected
    @test_approx_eq power(spectrogram(data, length(data), 0; window=window1(length(data)), onesided=false)) expected
end

# Padded periodogram
# MATLAB: a = periodogram(0:7, [], 32);
expected = [
                  98
    174.463067389405
    121.968086934209
    65.4971744936088
    27.3137084989848
    12.1737815028909
    10.3755170959439
    10.4034038628775
                   8
    5.25810953219633
    4.47015397150535
    4.89522578856669
    4.68629150101524
    3.69370284475603
     3.1862419983415
    3.61553458569862
                   2
]
@test_approx_eq power(periodogram(data; nfft=32)) expected
@test_approx_eq power(welch_pgram(data, length(data), 0; nfft=32)) expected
@test_approx_eq power(spectrogram(data, length(data), 0; nfft=32)) expected

# Padded periodogram with window
# MATLAB: a = periodogram(0:7, hamming(8), 32, 1)
expected = [
      65.4616239868015
      122.101693164395
      98.8444689598445
       69.020252632913
      41.1135835910315
      20.5496474310966
      8.43291449161938
      2.78001620362588
     0.738626287301088
     0.174995741770789
    0.0501563022944516
    0.0327357460012861
    0.0443348932217643
    0.0553999745503552
    0.0561319901616643
    0.0526025934871384
    0.0255029855641069
]
@test_approx_eq power(periodogram(data; window=hamming, nfft=32)) expected
@test_approx_eq power(welch_pgram(data, length(data), 0; window=hamming, nfft=32)) expected
@test_approx_eq power(spectrogram(data, length(data), 0; window=hamming, nfft=32)) expected

# Test fftshift
p = periodogram(data)
@test power(p) == power(fftshift(p))
@test_approx_eq freq(p) freq(fftshift(p))

p = periodogram(data; onesided=false)
@test fftshift(power(p)) == power(fftshift(p))
@test fftshift(freq(p)) == freq(fftshift(p))

data = 1:100

p = spectrogram(data)
@test power(p) == power(fftshift(p))
@test_approx_eq freq(p) freq(fftshift(p))

p = spectrogram(data; onesided=false)
@test fftshift(power(p), 1) == power(fftshift(p))
@test fftshift(freq(p)) == freq(fftshift(p))
