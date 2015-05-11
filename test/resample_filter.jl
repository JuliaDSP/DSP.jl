using DSP

rerate = 11//3      # Input sample rate
xƒ1    = 0.125      # First singal frequency
xƒ2    = 0.3        # Second signal frequency
xLen   = 80         # Number of signal samples
xTime  = [0:xLen-1] # Time vector
x      = cos(2*pi*0.125*xTime) + 0.5sin(2*pi*0.3*xTime*pi) + cos(0.1*xTime)

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
