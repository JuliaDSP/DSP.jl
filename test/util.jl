using DSP, Base.Test, Compat

@test_approx_eq(unwrap([0.1, 0.2, 0.3, 0.4]), [0.1, 0.2, 0.3, 0.4])
@test_approx_eq(unwrap([0.1, 0.2 + 2pi, 0.3, 0.4]), [0.1, 0.2, 0.3, 0.4])
@test_approx_eq(unwrap([0.1, 0.2 - 2pi, 0.3, 0.4]), [0.1, 0.2, 0.3, 0.4])
@test_approx_eq(unwrap([0.1, 0.2 - 2pi, 0.3 - 2pi, 0.4]), [0.1, 0.2, 0.3, 0.4])
@test_approx_eq(unwrap([0.1 + 2pi, 0.2, 0.3, 0.4]),
                [0.1 + 2pi, 0.2 + 2pi, 0.3 + 2pi, 0.4 + 2pi])
@test_approx_eq(unwrap([0.1, 0.2 + 6pi, 0.3, 0.4]), [0.1, 0.2, 0.3, 0.4])

test_v = [0.1, 0.2, 0.3 + 2pi, 0.4]
res_v = unwrap(test_v)
@test_approx_eq(test_v, [0.1, 0.2, 0.3 + 2pi, 0.4])
unwrap!(test_v)
@test_approx_eq(test_v, [0.1, 0.2, 0.3, 0.4])

# test multi-dimensional unwrapping
wrapped = [0.1, 0.2 + 2pi, 0.3, 0.4]
unwrapped = [0.1, 0.2, 0.3, 0.4]
wrapped = hcat(wrapped, wrapped)
unwrapped = hcat(unwrapped, unwrapped)
@test_approx_eq(unwrap(wrapped), wrapped)
@test_approx_eq(unwrap(wrapped, 1), unwrapped)

# test unwrapping with other ranges
unwrapped = [1.0:100;]
wrapped = unwrapped % 10
@test_approx_eq(unwrap(wrapped, range=10), unwrapped)

# Testing hilbert transform
t = (0:1/256:2-1/256)
a0 = sinpi(t)
a1 = cospi(t)
a2 = sinpi(2*t)
a3 = cospi(2*t)
a = hcat(a0, a1, a2, a3)

h = hcat(hilbert(a0), hilbert(a1), hilbert(a2), hilbert(a3))
h_abs = abs(h)
h_angle = angle(h)
h_real = real(h)

#The real part should be equal to the original signals:
@test_approx_eq h_real a

#The absolute value should be one everywhere, for this input:
@test_approx_eq h_abs ones(size(a))

#For the 'slow' sine - the phase should go from -pi/2 to pi/2 in
#the first 256 bins:
@test_approx_eq h_angle[1:256,1] -pi/2:pi/256:pi/2-pi/256

#For the 'slow' cosine - the phase should go from 0 to pi in the
#same interval:
@test_approx_eq h_angle[1:256,2] 0:pi/256:pi-pi/256

#The 'fast' sine should make this phase transition in half the time:
@test_approx_eq h_angle[1:128,3] -pi/2:pi/128:pi/2-pi/128
                    
#Ditto for the 'fast' cosine:
@test_approx_eq h_angle[1:128,4] 0:pi/128:pi-pi/128

#The imaginary part of hilbert(cos(t)) = sin(t) Wikipedia
@test_approx_eq imag(h[:,2]) a0

#Sanity check with odd number of samples
h2 = hilbert([ones(10); zeros(9)]) 
@test_approx_eq real(h2) [ones(10); zeros(9)]

#Sanity check with integer arguments
r = round(Int, rand(128)*20)
@compat @test hilbert(r) == hilbert(map(Float64, r))

# Test hilbert with 2D input
@test_approx_eq h hilbert(a)

## FFTFREQ

@test_approx_eq fftfreq(1) [0.]
@test_approx_eq fftfreq(2) [0., -1/2]
@test_approx_eq fftfreq(2, 1/2) [0., -1/4]
@test_approx_eq fftfreq(3) [0., 1/3, -1/3]
@test_approx_eq fftfreq(3, 1/2) [0., 1/6, -1/6]
@test_approx_eq fftfreq(6) [0., 1/6, 1/3, -1/2, -1/3, -1/6]
@test_approx_eq fftfreq(7) [0., 1/7, 2/7, 3/7, -3/7, -2/7, -1/7]

@test_approx_eq rfftfreq(1) [0.]
@test_approx_eq rfftfreq(2) [0., 1/2]
@test_approx_eq rfftfreq(2, 1/2) [0., 1/4]
@test_approx_eq rfftfreq(3) [0., 1/3]
@test_approx_eq rfftfreq(3, 1/2) [0., 1/6]
@test_approx_eq rfftfreq(6) [0., 1/6, 1/3, 1/2]
@test_approx_eq rfftfreq(7) [0., 1/7, 2/7, 3/7]

for n = 1:7
    @test_approx_eq fftshift(fftfreq(n)) fftshift([fftfreq(n);])
end

# nextfastfft
@test nextfastfft(64) == 64
@test nextfastfft(65) == 70
@test nextfastfft(127) == 128
@test nextfastfft((64,65,127)) == (64,70,128)
@test nextfastfft(64,65,127) == nextfastfft((64,65,127))


## COMMON DSP TOOLS

# dB conversion
@test_approx_eq 3dB db2pow(3)
@test_approx_eq -3dB db2pow(-3)
@test_approx_eq 3dBa db2amp(3)
@test_approx_eq -3dBa db2amp(-3)
@test isa(3e0dB, Float64)
@test isa(3f0dB, Float32)
test_num = convert(Float64, pi)
@test_approx_eq pow2db(test_num) 10*log10(test_num)
@test_approx_eq amp2db(test_num) 20*log10(test_num)
@test_approx_eq test_num*dB db2pow(test_num)
@test_approx_eq test_num*dBa db2amp(test_num)
@test_approx_eq test_num db2pow(pow2db(test_num))
@test_approx_eq test_num db2amp(amp2db(test_num))

n = (5,6)
for x in ( randn(n), randn(n)+randn(n)im )
    @test_approx_eq rms(x) sqrt(mean(abs(x).^2))
    @test_approx_eq rmsfft(fft(x)) rms(x)
end # for

@test_approx_eq shiftin!([1,2,3,4],[5,6]) [3,4,5,6]

## DELAY TOOLS

# Test Variables
max_delay = 1025    # Maximum delay for copied shifted signals
max_order = 1025    # Maximum order for FIR shifting filters

# Test delays that do not require to erase data to align the signals (delay
# shorter than 0-pads around data)
x = [zeros(max_delay); 2 * (rand(4 * max_delay) - (1 / 2)); zeros(max_delay)]
delays = 1:(max_delay - 1)
x_delay = zeros(x)
x_advance = zeros(x)

for d in delays

    fill!(x_delay, 0)
    fill!(x_advance, 0)

    x_delay[(d + 1):end] = x[1:(end - d)]
    x_advance[1:(end - d)] = x[(d + 1):end]

    @test shiftsignals(x, d) == x_advance
    @test d == finddelay(x, x_delay)
    @test (x, d) == alignsignals(x, x_delay)
    @test shiftsignals(x, -d) == x_delay
    @test -d == finddelay(x, x_advance)
    @test (x, -d) == alignsignals(x, x_advance)
    shiftsignals!(x_advance, -d)
    @test x == x_advance
    shiftsignals!(x_delay, d)
    @test x == x_delay

end

# Test data erasing delays (when the delay is longer than the 0-pads): 
x = [zeros(fld(max_delay, 2)); 2 * (rand(4 * max_delay) - (1 / 2)); 
    zeros(fld(max_delay, 2))]
delays = (fld(max_delay, 2)):(max_delay - 1)
x_delay = zeros(x)
x_advance = zeros(x)

for d in delays

    fill!(x_delay, 0)
    fill!(x_advance, 0)

    x_delay[(d + 1):end] = x[1:(end - d)]
    x_advance[1:(end - d)] = x[(d + 1):end]

    @test shiftsignals(x, d) == x_advance
    @test d == finddelay(x, x_delay)
    @test ([x[1:(end - d)]; zeros(d)], d) == alignsignals(x, x_delay)
    @test shiftsignals(x, -d) == x_delay
    @test -d == finddelay(x, x_advance)
    @test ([zeros(d); x[(d + 1):end]], -d) == alignsignals(x, x_advance)
    shiftsignals!(x_advance, -d)
    @test [zeros(d); x[(d + 1):end]] == x_advance
    shiftsignals!(x_delay, d)
    @test [x[1:(end - d)]; zeros(d)] == x_delay

end

# Test delay calculation for not identical signals. A FIR filter will supply
# a known delay easily:
orders = 2 * (4:(fld(max_order, 2)))    # Make even to remove delay ambiguity
fir_delays = round(Int, orders / 2)
fir_taps = orders + 1

x = [zeros(maximum(fir_delays)); 2 * (rand(4 * maximum(fir_delays)) - (1 / 2)); 
    zeros(maximum(fir_delays))]
y = zeros(x)

for n = 1:length(orders)

    w = hamming(fir_taps[n])
    dsm = FIRWindow(w)
    rty = Lowpass(0.1) # Too small bandwidth can create errors of few samples
    fir = digitalfilter(rty, dsm)
    
    y = fftfilt(fir, x)
    
    @test fir_delays[n] == finddelay(x, y)

end
