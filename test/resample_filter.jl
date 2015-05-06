using DSP
using PyPlot

ƒsIn   = 1.0                               # Input sample rate
ƒsOut  = π                                 # Input sample rate
rerate = ƒsOut/ƒsIn                        # Resampling ratio
xƒ1    = 0.125                             # First singal frequency
xƒ2    = 0.3                               # Second signal frequency
xLen   = 30                                # Number of signal samples
xTime  = [0:xLen-1]                        # Time vector
x      = cos(2*pi*xƒ1*xTime)               # Create a signal vector
x      = x + 0.5sin(2*pi*xƒ2*xTime*pi)   
x      = x + cos(0.1*xTime)

tapsPerϕ = 11                              # Taps per polyphase filter, or phase
Nϕ       = 32                              # Number of filters, or phases
hLen     = tapsPerϕ*Nϕ                     # Total number of taps on coef vector
relBW    = 0.8                             # Relative bandwidth. Or where to set the cutoff frequency in realtion to the nyquist frequency
ƒnyq     = rerate > 1.0 ? 0.5/Nϕ : rate/Nϕ # Nyquist frequency. We must filter out everything below this point to avoid aliases
ƒc       = ƒnyq * relBW                    # Cutoff frequency
tw       = (1.0-relBW) * ƒnyq              # Transistion with

h = digitalfilter(Lowpass(ƒc), FIRWindow(hamming(hLen))) # Create filter taps
scale!(h, Nϕ)                                            # Scale h to Nϕ to account for interpolaton 

pf      = FIRFilter(h, rerate)
δfilter = pf.kernel.tapsPerϕ/2   # The time delay (at output rate) caused by the filtering process

(phase,throwaway) = modf(δfilter)
pf.kernel.inputDeficit += throwaway
setphase!(pf, phase)

y     = filt(pf, x)
yTime = [0.0:length(y)-1]/rerate # Create y time vector. Accout for filter delay so the plots line up


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