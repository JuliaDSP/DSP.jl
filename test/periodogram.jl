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
