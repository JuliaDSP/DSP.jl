using DSP
using PyPlot

ƒsIn   = 1.0                           # Input sample rate
ƒsOut  = π                             # Input sample rate
rerate = ƒsOut/ƒsIn
xƒ1    = 0.125                         # First singal frequency
xƒ2    = 0.3                           # Second signal frequency
xLen   = 100                           # Number of signal samples
xTime  = [0:xLen-1]                    # Signal time vector
x      = cos(2*pi*xƒ1*xTime)
x      = x + 0.5sin(2*pi*xƒ2*xTime*pi) # Create the two signals and add them
x      = x + cos(0.1*xTime)

pf      = FIRFilter(rerate)
δfilter = (pf.kernel.tapsPerϕ-1)/2     # ~(hLen-1)/(2*Nϕ) The time delay (at output rate) caused by the filtering process

(phase,throwaway) = modf(δfilter) 
pf.kernel.inputDeficit += throwaway

y     = filt(pf, x)
yTime = [0.0:length(y)-1]/rerate       # Create y time vector. Accout for filter delay so the plots line up


figure(num=1, figsize=(10, 10/golden), dpi=100, facecolor="w", edgecolor="k" )
clf()
rc( "font", size=10 )
hold(true)

plt.suptitle( "Arbitrary Polphase Resampling, ƒsIn = $(ƒsIn), ƒsOut = $(ƒsOut)" )

subplot( 311 )

plot( xTime, x, "b")
stem( xTime, x, linefmt = "b-", markerfmt = "b." )
xlabel( "Time" )
ylabel( "Amplitude" )

subplot( 312 )
stem( yTime, y, linefmt = "r-", markerfmt = "r." )
plot( yTime, y, "r-" )
xlabel( "Time" )
ylabel( "Amplitude" )
xlim( yTime[1], xTime[end]+yTime[1] )

subplot( 313 )
plot( xTime, x, "b.-")
plot( yTime, y, "r.-" )

hold(false)
