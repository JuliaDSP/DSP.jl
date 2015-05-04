using DSP
using Winston

h = resample_filter(147//160, 0.9)
f = linspace(-1,1,length(h))
H = abs2(fftshift(fft(h)))
p = plot(f, H)

display(p)

x = [cos(2*pi*0.25*x) for x in 0:99]
p = plot(x)
display(p)

y = resample(x, 4.0)
p = plot(y)
display(p)