using DSP, Base.Test, Polynomial
import Base.Sort.Lexicographic
import DSP.FilterDesign: coeffs

function lt(a, b)
    if abs(real(a) - real(b)) > 1e-10
        isless(real(a), real(b))
    else
        isless(imag(a), imag(b))
    end
end

function tffilter_eq(f1, f2)
    b1, a1 = (coeffs(f1.b), coeffs(f1.a))
    b2, a2 = (coeffs(f2.b), coeffs(f2.a))
    @test_approx_eq float64(b1) float64(b2)
    @test_approx_eq float64(a1) float64(a2)
end

function zpkfilter_eq(f1, f2)
    @test_approx_eq complex128(sort(f1.z, lt=lt)) complex128(sort(f2.z, lt=lt))
    @test_approx_eq complex128(sort(f1.p, lt=lt)) complex128(sort(f2.p, lt=lt))
    @test_approx_eq float64(f1.k) float64(f2.k)
end

function zpkfilter_eq(f1, f2, eps)
    @test_approx_eq_eps complex128(sort(f1.z, lt=lt)) complex128(sort(f2.z, lt=lt)) eps
    @test_approx_eq_eps complex128(sort(f1.p, lt=lt)) complex128(sort(f2.p, lt=lt)) eps
    @test_approx_eq_eps float64(f1.k) float64(f2.k) eps
end

function tffilter_accuracy(f1, f2, accurate_f)
    b1, a1 = (coeffs(f1.b), coeffs(f1.a))
    b2, a2 = (coeffs(f2.b), coeffs(f2.a))
    accurate_b, accurate_a = (coeffs(accurate_f.b), coeffs(accurate_f.a))
    @test sum(abs(b1 - accurate_b)) <= sum(abs(b2 - accurate_b))
    @test sum(abs(a1 - accurate_a)) <= sum(abs(a2 - accurate_a))
end

#
# Butterworth filter prototype
#

# Poles of 20 pole Butterworth filter prototype from MATLAB 2013b (buttap(20))
m_prototype = [-0.07845909572784487+0.996917333733128im,-0.07845909572784487-0.996917333733128im,-0.2334453638559053+0.9723699203976767im,-0.2334453638559053-0.9723699203976767im,-0.3826834323650897+0.9238795325112867im,-0.3826834323650897-0.9238795325112867im,-0.5224985647159488+0.8526401643540923im,-0.5224985647159488-0.8526401643540923im,-0.6494480483301835+0.760405965600031im,-0.6494480483301835-0.760405965600031im,-0.7604059656000308+0.6494480483301838im,-0.7604059656000308-0.6494480483301838im,-0.8526401643540922+0.5224985647159489im,-0.8526401643540922-0.5224985647159489im,-0.9238795325112867+0.3826834323650899im,-0.9238795325112867-0.3826834323650899im,-0.9723699203976766+0.2334453638559055im,-0.9723699203976766-0.2334453638559055im,-0.996917333733128+0.07845909572784507im,-0.996917333733128-0.07845909572784507im]

# Test that our answers are close to MATLAB's
f = Butterworth(20)
prototype = sort(f.p, lt=lt)
m_prototype = sort(m_prototype, lt=lt)
@test isempty(f.z)
@test_approx_eq prototype m_prototype
@test f.k == 1

# Test that our answers are more accurate than MATLAB's
x = [BigFloat(2*i-1) for i = 1:20]/(2*20)
accurate_prototype = sort(complex(-sinpi(x), cospi(x)), lt=lt)
accurate_f = ZPKFilter(BigFloat[], accurate_prototype, BigFloat(1))
@test_approx_eq prototype complex128(accurate_prototype)
@test sum(abs(prototype - accurate_prototype)) <= sum(abs(m_prototype - accurate_prototype))

#
# Conversion between zpk and tf
#

# Output of [z, p, k] = buttap(20); [b, a] = zp2tf(z, p, k)
m_a = [1,12.74549484318237,81.22381939879423,343.6513712403923,1081.352361133001,2687.409807920676,5468.931438945091,9326.061201886809,13528.36656744904,16852.27707949905,18122.54155403868,16852.27707949905,13528.36656744904,9326.061201886809,5468.931438945092,2687.409807920676,1081.352361133001,343.6513712403923,81.22381939879423,12.74549484318237,1]

f = convert(TFFilter, Butterworth(20))
@test_approx_eq coeffs(f.b) [1]
@test_approx_eq coeffs(f.a) m_a

# Test that our answers are more accurate than MATLAB's
accurate_a = coeffs(convert(TFFilter, accurate_f).a)
@test_approx_eq coeffs(f.a) float64(accurate_a)
@test sum(abs(coeffs(f.a) - accurate_a)) <= sum(abs(m_a - accurate_a))

#
# Conversion between tf and zpk
#

f = convert(ZPKFilter, convert(TFFilter, Butterworth(20)))
@test isempty(f.z)
@test_approx_eq_eps sort(f.p, lt=lt) prototype 1e-7

# Test that our answers are more accurate than MATLAB's
# Output of [z, p, k] = buttap(20); [b, a] = zp2tf(z, p, k); tf2zpk(b, a)
m_p = [-0.07845909573254482+0.9969173337335029im,-0.07845909573254482-0.9969173337335029im,-0.2334453637958131+0.9723699203918822im,-0.2334453637958131-0.9723699203918822im,-0.3826834327796701+0.9238795325396184im,-0.3826834327796701-0.9238795325396184im,-0.5224985628488221+0.8526401643454914im,-0.5224985628488221-0.8526401643454914im,-0.6494480541398985+0.7604059651905597im,-0.6494480541398985-0.7604059651905597im,-0.760405952587916+0.6494480502272874im,-0.760405952587916-0.6494480502272874im,-0.8526401859847815+0.5224985598169277im,-0.8526401859847815-0.5224985598169277im,-0.9238795057196649+0.3826834415328767im,-0.9238795057196649-0.3826834415328767im,-0.9969173244841298+0.07845911266921719im,-0.9969173244841298-0.07845911266921719im,-0.97236994351794+0.2334453500964366im,-0.97236994351794-0.2334453500964366im]
@test sum(abs(sort(f.p, lt=lt) - accurate_prototype)) .<= sum(abs(m_p - accurate_prototype))

#
# Frequency scaling
#

# Output of [b, a] = butter(20, 0.5, 's')
m_b = [9.536743164062593e-07]
m_a = [1,6.372747421591193,20.30595484969859,42.95642140504915,67.58452257081277,83.98155649752141,85.4520537335174,72.85985313974101,52.84518190409806,32.91460367089675,17.697794486366,8.228650917724192,3.302823869006134,1.138435205308456,0.3337970848965534,0.08201323876711111,0.01650012758076492,0.002621851892398036,0.0003098442817642021,2.43101021636629e-05,9.536743164062593e-07]
m_f = TFFilter(m_b, m_a)

f = analogfilter(Lowpass(0.5), Butterworth(20))
tffilter_eq(f, m_f)

# Test that our answers are more accurate than MATLAB's
accurate_ft = analogfilter(Lowpass(0.5), accurate_f)
tffilter_eq(f, accurate_ft)
tffilter_accuracy(f, m_f, accurate_ft)

#
# High pass filter construction
#

# Output of [b, a] = butter(20, 0.5, 'high', 's')
m_b = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
m_a = [1,6.372747421591192,20.30595484969859,42.95642140504915,67.58452257081278,83.98155649752141,85.4520537335174,72.85985313974103,52.84518190409808,32.91460367089676,17.697794486366,8.228650917724194,3.302823869006133,1.138435205308455,0.333797084896553,0.08201323876711097,0.01650012758076489,0.002621851892398029,0.000309844281764201,2.431010216366281e-05,9.536743164062553e-07]
m_f = TFFilter(m_b, m_a)

f = analogfilter(Highpass(0.5), Butterworth(20))
tffilter_eq(f, m_f)

# Test that our answers are more accurate than MATLAB's
accurate_ft = analogfilter(Highpass(0.5), accurate_f)
tffilter_eq(f, accurate_ft)
tffilter_accuracy(f, m_f, accurate_ft)

#
# Band pass filter construction
#

# Output of [b, a] = butter(10, 0.5, 'high', 's')
m_b = [0.0009765624999999457,-0,-0,-0,-0,-0,-0,-0,-0,-0,-0]
m_a = [1,3.19622661074983,6.982932273632673,10.74389003924683,13.29907942734037,13.38723211262245,11.39496738633279,8.228979126528134,5.12370722185132,2.74825868033936,1.279689995802157,0.5152985025636301,0.1801303320182105,0.05424375889068843,0.01408374570149165,0.003102395442359219,0.000577868753822514,8.753272710678296e-05,1.066714028066583e-05,9.154803174808179e-07,5.370475264498965e-08]
m_f = TFFilter(m_b, m_a)

f = analogfilter(Bandpass(0.25, 0.75), Butterworth(10))
tffilter_eq(f, m_f)

# Test that our answers are more accurate than MATLAB's
x = [BigFloat(2*i-1) for i = 1:10]/(2*10)
accurate_prototype = sort(complex(-sinpi(x), cospi(x)), lt=lt)
accurate_f10 = DSP.FilterDesign.ZPKFilter(BigFloat[], accurate_prototype, BigFloat(1))
accurate_ft = analogfilter(Bandpass(0.25, 0.75), accurate_f10)
tffilter_eq(f, accurate_ft)
tffilter_accuracy(f, m_f, accurate_ft)

#
# Band stop filter construction
#

# Output of [b, a] = butter(10, 0.5, 'stop', 's')
m_b = [0.9999999999999911,0,1.874999999999983,0,1.582031249999986,0,0.7910156249999928,0,0.2595520019531226,0,0.05839920043945258,0,0.009124875068664466,0,0.0009776651859283356,0,6.87420833855861e-05,0,2.86425347439942e-06,0,5.370475264498911e-08]
m_a = [1,3.196226610749829,6.982932273632667,10.74389003924681,13.29907942734034,13.38723211262241,11.39496738633274,8.228979126528094,5.12370722185129,2.748258680339342,1.279689995802148,0.5152985025636261,0.180130332018209,0.05424375889068796,0.01408374570149152,0.003102395442359189,0.0005778687538225081,8.753272710678211e-05,1.066714028066572e-05,9.154803174808088e-07,5.370475264498911e-08]
m_f = TFFilter(m_b, m_a)

f = analogfilter(Bandstop(0.25, 0.75), Butterworth(10))
tffilter_eq(f, m_f)

# Test that our answers are more accurate than MATLAB's
accurate_ft = analogfilter(Bandstop(0.25, 0.75), accurate_f10)
tffilter_eq(f, accurate_ft)
tffilter_accuracy(f, m_f, accurate_ft)

#
# Digital filter creation
#

# Output of [b, a] = butter(20, [0.5], 'low')
m_b = [8.555155699386467e-06,0.0001711031139877293,0.001625479582883429,0.009752877497300572,0.04144972936352743,0.1326391339632878,0.3315978349082195,0.663195669816439,1.077692963451713,1.436923951268951,1.580616346395846,1.436923951268951,1.077692963451713,0.663195669816439,0.3315978349082195,0.1326391339632878,0.04144972936352743,0.009752877497300572,0.001625479582883429,0.0001711031139877293,8.555155699386467e-06]
m_a = [1,-9.216626485446244e-15,2.719352210996272,-2.276898740746355e-14,2.954502383565172,-2.181204855972978e-14,1.663404391798304,-9.886771142150971e-15,0.5272563634534938,-1.578399203426284e-15,0.09588816575455775,3.201156931050434e-16,0.009788250615318749,1.651392612127069e-16,0.0005259945756568433,3.034673947897289e-17,1.306795536212839e-05,3.205637858266449e-18,1.137794123801438e-07,9.233765069631315e-20,1.463813780864114e-10]
m_f = TFFilter(m_b, m_a)

f = convert(TFFilter, digitalfilter(Lowpass(0.5), Butterworth(20)))
tffilter_eq(f, m_f)

# Cannot test accuracy since the roots function in Polynomial ends up
# calling BLAS and so cannot preserve accuracy for BigFloats

#
# Conversion to second-order sections
#

x = randn(100)

# Test that biquad filt() works
for ftype in (Lowpass(0.5), Highpass(0.5))
    f = digitalfilter(ftype, Butterworth(1))
    @test_approx_eq filt(f, x) filt(SOSFilter([convert(BiquadFilter, convert(TFFilter, f))], 1.0), x)
    f = digitalfilter(ftype, Butterworth(2))
    @test_approx_eq filt(f, x) filt(SOSFilter([convert(BiquadFilter, convert(TFFilter, f))], 1.0), x)
end

# Test that filters converted to SOS are equivalent to Butterworth filters
for ftype in (Lowpass(0.5), Highpass(0.5), Bandpass(0.25, 0.75), Bandstop(0.25, 0.75))
    for order in 1:4
        f = digitalfilter(ftype, Butterworth(order))
        @test_approx_eq filt(f, x) filt(convert(SOSFilter, f), x)
    end
end

# Test that a numerically challenging filter (high order, clustered roots) has acceptable errors in its coefficients after conversion to SOS
f = ZPKFilter(ones(100), 0.99*ones(100), 1)
g = convert(SOSFilter, f)
tffilter_eq(convert(TFFilter, f), convert(TFFilter, g))

# Test designing filters as SOS
for ftype in (Lowpass(0.5), Highpass(0.5), Bandpass(0.25, 0.75), Bandstop(0.25, 0.75))
    for order in 1:(isa(ftype, Lowpass) || isa(ftype, Highpass) ? 4 : 2)
        f1 = digitalfilter(ftype, Butterworth(order))
        f2 = digitalfilter(ftype, convert(SOSFilter, Butterworth(order)))
        @assert isa(f2, SOSFilter)
        @test_approx_eq filt(f1, x) filt(f2, x)
    end
end

#
# Filter conversion tests
#

for f in (digitalfilter(Lowpass(0.5), Butterworth(1)), digitalfilter(Lowpass(0.5), Butterworth(2)),
          digitalfilter(Bandpass(0.25, 0.75), Butterworth(1)))
    for ftype1 in (ZPKFilter, TFFilter, BiquadFilter, SOSFilter)
        f2 = convert(ftype1, f)
        for ftype2 in (ZPKFilter, TFFilter, BiquadFilter, SOSFilter)
            f3 = convert(ftype2, f)
            zpkfilter_eq(f, convert(ZPKFilter, f3), eps())
        end
    end
end

f = digitalfilter(Lowpass(0.5), Butterworth(3))
for ftype1 in (ZPKFilter, TFFilter, SOSFilter)
    f2 = convert(ftype1, f)
    for ftype2 in (ZPKFilter, TFFilter, SOSFilter)
        f3 = convert(ftype2, f2)
        zpkfilter_eq(f, convert(ZPKFilter, f3), 1e-5)
    end
end
