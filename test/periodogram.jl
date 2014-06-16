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

x0 = readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_x.txt"),'\t')
f0 = readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_f.txt"),'\t')
t0 = readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_t.txt"),'\t')
p0 = readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_p.txt"),'\t')
p, t, f = spectrogram(x0, n=256, r=10, m=128)

# with real input matlab outputs a 1-sided PSD
@test_approx_eq p0[[1,129],:] p[[1,129],:]
@test_approx_eq p0[2:128,:]./2 p[2:128,:]
@test_approx_eq f0 f[1:129]
@test_approx_eq vec(t0) t


# # ~~~~ TESTS FOR DSP.Periodogram.welch_pgram ~~~~
data = Float64[0:7]

# ~~~~~~~~~~~ This one tests periodogram ~~~~~~~~~~
#Matlab: p = pwelch(0:7, [1, 1, 1, 1, 1, 1, 1, 1], 0, 8, 1, 'twosided')
@test_approx_eq welch_pgram(data, length(data), 0)  Float64[ 98.0,
                                                             13.656854249492380,
                                                              4.0,
                                                              2.343145750507620,
                                                              2.0,
                                                              2.343145750507620,
                                                              4.0,
                                                             13.656854249492380]

# ~~~~~~~~ Tests with no window ~~~~~~~~~~~~~~~~~~~
# Matlab: p = pwelch(0:7, [1, 1], 0, 2, 1, 'twosided')
@test_approx_eq welch_pgram(data, 2, 0)  Float64[34.5, 0.5]

# Matlab: p = pwelch(0:7, [1, 1, 1], 0, 3, 1, 'twosided')
@test_approx_eq welch_pgram(data, 3, 0)  Float64[25.5, 1.0, 1.0]

# Matlab: p = pwelch(0:7, [1, 1, 1], 1, 3, 1, 'twosided')
@test_approx_eq welch_pgram(data, 3, 1)  Float64[35.0, 1.0, 1.0]

# Matlab: p = pwelch(0:7, [1, 1, 1, 1], 1, 4, 1, 'twosided')
@test_approx_eq welch_pgram(data, 4, 1)  Float64[45, 2, 1, 2]

# ~~~~~~~~~ Tests with window ~~~~~~~~~~~~~~~
data = Float64[0:7]

# ~~~~~~~~~~~ This one tests periodogram ~~~~~~~~~~~~
# ~ If functionality of the other arguments has been 
# ~ tested above, we only test here that the correct 
# ~ value of the spectral density is obtained when 
# ~ using a window. More tests to be added if needed
#Matlab: p = pwelch(0:7, window_func(8), 0, 8, 1, 'twosided')
cases = {hamming => Float64[65.461623986801527,
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
                            21.981076052592442]					 }

for (window, expected) in cases
	@test_approx_eq welch_pgram(data, length(data), 0, window) expected
end
