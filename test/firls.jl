using DSP

# Lowpass
N    = 41
b    = [(0.0, 0.4, 1), (0.6, 1, 0)]
h_lp = firls(N, b)

# Highpass
N    = 41
b    = [(0.0, 0.4, 0), (0.6, 1, 1)]
h_hp = firls(N, b)

# Bandpass
N    = 41
b    = [(0, .25, -50dBa), (.4, .6, 0dBa), (.75, 1, -50dBa)]
h_bp = firls(N, b)

# Bandstop
N    = 41
b    = [(0, .2, 0dBa), (.4, .6, -10dBa), (.8, 1, 0dBa)]
h_bs = firls(N, b)

using PyPlot

function Filters.freqz(h::Vector, w = linspace(0, π, 250))
    pr = PolynomialRatio(h, [one(eltype(h))])
    freqz(pr, w)
end
Util.amp2db(x::Complex) = amp2db(abs(x))
Util.amp2db(x::AbstractVector) = [amp2db(x) for x in x]

clf()
subplot(2,2,1)
title("Lowpass")
ylabel("dB")
xlabel("Normalized ƒ (π rad/s)")
H = amp2db(freqz(h_lp))
f = linspace(0,1, length(H))
plot(f, H)
grid()

subplot(2,2,2)
title("Highpass")
ylabel("dB")
xlabel("Normalized ƒ (π rad/s)")
H = amp2db(freqz(h_hp))
f = linspace(0,1, length(H))
plot(f, H)
grid()

subplot(2,2,3)
title("Bandpass")
ylabel("dB")
xlabel("Normalized ƒ (π rad/s)")
H = amp2db(freqz(h_bp))
f = linspace(0,1, length(H))
plot(f, H)
grid()

subplot(2,2,4)
title("Bandstop")
ylabel("dB")
xlabel("Normalized ƒ (π rad/s)")
H = amp2db(freqz(h_bs))
f = linspace(0,1, length(H))
plot(f, H)
grid()
