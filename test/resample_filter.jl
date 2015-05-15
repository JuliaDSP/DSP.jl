using DSP
# MATLAB
# t = [0:100]
# x = (1 + sin(2*pi*0.005*t)) .* sin(2*pi*.05*t)

# AM Modulator

sig(t) = [(1 + sin(2π*0.005*t)) * sin(2π*.05*t) for t in t]

rerate = 3//2       # Input sample rate
xLen   = 100        # Number of signal samples
xTime  = [0:100]    # Create time vector    
x      = sig(xTime) # Cretae signal vector

# Resample with arbitrary factor
yArb    = resample(x, float64(rerate))
arbTime = [0.0:length(yArb)-1]/float64(rerate)

# Resample with rational factor
yRat    = resample(x, rerate)
ratTime = [0.0:length(yRat)-1]/float64(rerate)


using PyPlot

figure(num=1, figsize=(10, 10/golden), dpi=100, facecolor="w", edgecolor="k")
clf()
rc("font", size=10)

hold(true)

subplot(211)
plt.title("Arbitrary Polphase Resampling, rate = $(float64(rerate))")
plot(xTime, x, "g.-")
plot(arbTime, yArb, "b.-")

subplot(212)
plt.title("Rational Polphase Resampling, rate = $rerate")
plot(xTime, x, "g.-")
plot(ratTime, yRat, "b.-")

hold(false)
