using DSP, Base.Test

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
r = int(rand(128)*20)
@test hilbert(r) == hilbert(float64(r))

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
num = float64(pi)
@test_approx_eq pow2db(num) 10*log10(num)
@test_approx_eq amp2db(num) 20*log10(num)
@test_approx_eq num*dB db2pow(num)
@test_approx_eq num*dBa db2amp(num)
@test_approx_eq num db2pow(pow2db(num))
@test_approx_eq num db2amp(amp2db(num))

n = (5,6)
for x in ( randn(n), randn(n)+randn(n)im )
    @test_approx_eq rms(x) sqrt(mean(abs(x).^2))
    @test_approx_eq rmsfft(fft(x)) rms(x)
end # for


