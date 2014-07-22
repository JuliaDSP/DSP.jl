using DSP, Base.Test, Polynomials
import Base.Sort.Lexicographic

function lt(a, b)
    if abs(real(a) - real(b)) > 1e-10
        isless(real(a), real(b))
    else
        isless(imag(a), imag(b))
    end
end

function tffilter_eq(f1, f2)
    b1, a1 = (coefb(f1), coefa(f1))
    b2, a2 = (coefb(f2), coefa(f2))
    @test_approx_eq float64(b1) float64(b2)
    @test_approx_eq float64(a1) float64(a2)
end

function zpkfilter_eq(f1, f2)
    if !isempty(f1.z) || !isempty(f2.z)
        @test_approx_eq complex128(sort(f1.z, lt=lt)) complex128(sort(f2.z, lt=lt))
    end
    @test_approx_eq complex128(sort(f1.p, lt=lt)) complex128(sort(f2.p, lt=lt))
    @test_approx_eq float64(f1.k) float64(f2.k)
end

function zpkfilter_eq(f1, f2, eps)
    if !isempty(f1.z) || !isempty(f2.z)
        @test_approx_eq_eps complex128(sort(f1.z, lt=lt)) complex128(sort(f2.z, lt=lt)) eps
    end
    @test_approx_eq_eps complex128(sort(f1.p, lt=lt)) complex128(sort(f2.p, lt=lt)) eps
    @test_approx_eq_eps float64(f1.k) float64(f2.k) eps
end

loss(x::Real, y::Real) = abs(float(x) - float(y))/eps(float(x))
loss(x::Union(Real,Complex), y::Union(Real,Complex)) = loss(real(x), real(y)) + loss(imag(x), imag(y))
loss(x::AbstractVector, y::AbstractVector) = sum(map(loss, x, y))

function accuracy_check(err1, err2, part, relerr=1)
    try
        @test err1 <= relerr*err2
    catch e
        println("Filter 1 $part error (ULP): ", float64(err1))
        println("Filter 2 $part error (ULP): ", float64(err2))
        println("Ratio: ", float64(err1/err2))
        rethrow(e)
    end
end

function tffilter_accuracy(f1, f2, accurate_f)
    b1, a1 = (coefb(f1), coefa(f1))
    b2, a2 = (coefb(f2), coefa(f2))
    accurate_b, accurate_a = (coefb(accurate_f), coefa(accurate_f))
    accuracy_check(loss(b1, accurate_b), loss(b2, accurate_b), "b")
    accuracy_check(loss(a1, accurate_a), loss(a2, accurate_a), "a")
end

function zpkfilter_accuracy(f1, f2, accurate_f; relerr=1, eps=nothing)
    z1, p1 = sort(f1.z, lt=lt), sort(f1.p, lt=lt)
    z2, p2 = sort(f2.z, lt=lt), sort(f2.p, lt=lt)
    accurate_z, accurate_p = sort(accurate_f.z, lt=lt), sort(accurate_f.p, lt=lt)
    if !isempty(z1) || !isempty(z2) || !isempty(accurate_z)
        if eps != nothing
            @test_approx_eq_eps z1 accurate_z eps
            @test_approx_eq_eps z2 accurate_z eps
        else
            @test_approx_eq z1 accurate_z
            @test_approx_eq z2 accurate_z
        end
        accuracy_check(loss(z1, accurate_z), loss(z2, accurate_z), "z", relerr)
    end
    if eps != nothing
        @test_approx_eq_eps p1 accurate_p eps
        @test_approx_eq_eps p2 accurate_p eps
        @test_approx_eq_eps f1.k accurate_f.k eps
        @test_approx_eq_eps f2.k accurate_f.k eps
    else
        @test_approx_eq p1 accurate_p
        @test_approx_eq p2 accurate_p
        @test_approx_eq f1.k accurate_f.k
        @test_approx_eq f2.k accurate_f.k
    end
    accuracy_check(loss(p1, accurate_p), loss(p2, accurate_p), "p", relerr)
    accuracy_check(loss(f1.k, accurate_f.k), loss(f2.k, accurate_f.k), "k", relerr)
end

#
# Butterworth filter prototype
#

# Poles of 20 pole Butterworth filter prototype from MATLAB 2013b (buttap(20))
matlab_p = [-0.07845909572784487+0.996917333733128im,-0.07845909572784487-0.996917333733128im,-0.2334453638559053+0.9723699203976767im,-0.2334453638559053-0.9723699203976767im,-0.3826834323650897+0.9238795325112867im,-0.3826834323650897-0.9238795325112867im,-0.5224985647159488+0.8526401643540923im,-0.5224985647159488-0.8526401643540923im,-0.6494480483301835+0.760405965600031im,-0.6494480483301835-0.760405965600031im,-0.7604059656000308+0.6494480483301838im,-0.7604059656000308-0.6494480483301838im,-0.8526401643540922+0.5224985647159489im,-0.8526401643540922-0.5224985647159489im,-0.9238795325112867+0.3826834323650899im,-0.9238795325112867-0.3826834323650899im,-0.9723699203976766+0.2334453638559055im,-0.9723699203976766-0.2334453638559055im,-0.996917333733128+0.07845909572784507im,-0.996917333733128-0.07845909572784507im]
matlab_butter = ZPKFilter([], matlab_p, 1)

# Test that our answers are close to MATLAB's and at least as accurate
butter = Butterworth(20)
zpkfilter_eq(butter, matlab_butter)
zpkfilter_accuracy(butter, matlab_butter, Butterworth(BigFloat, 20))

# Poles of 19 pole Butterworth filter prototype from MATLAB 2013b (buttap(19))
matlab_p = [-0.08257934547233227+0.9965844930066698im,-0.08257934547233227-0.9965844930066698im,-0.2454854871407991+0.9694002659393304im,-0.2454854871407991-0.9694002659393304im,-0.4016954246529694+0.9157733266550574im,-0.4016954246529694-0.9157733266550574im,-0.5469481581224267+0.8371664782625287im,-0.5469481581224267-0.8371664782625287im,-0.6772815716257409+0.7357239106731317im,-0.6772815716257409-0.7357239106731317im,-0.7891405093963935+0.6142127126896679im,-0.7891405093963935-0.6142127126896679im,-0.8794737512064892+0.4759473930370733im,-0.8794737512064892-0.4759473930370733im,-0.9458172417006346+0.3246994692046836im,-0.9458172417006346-0.3246994692046836im,-0.9863613034027224+0.1645945902807336im,-0.9863613034027224-0.1645945902807336im,-1+0im]

# Test that our answers are close to MATLAB's and at least as accurate
zpkfilter_eq(Butterworth(19), ZPKFilter([], matlab_p, 1))
zpkfilter_accuracy(Butterworth(19), ZPKFilter([], matlab_p, 1), Butterworth(BigFloat, 19))

#
# Chebyshev type I filter
#

# Poles of 20 pole Butterworth filter prototype with 1 dB passband ripple from MATLAB 2013b:
#=
    [z, p, k] = cheb1ap(20, 1)
    sprintf('%.16g%+.16gim,', [real(p) imag(p)]')
=#
matlab_p = [-0.005606643513655412+0.9994594480354867im,-0.01668187637027696+0.9748494394091176im,-0.02734634606668882+0.9262354022446964im,-0.03733745796561948+0.854814375477951im,-0.04640919773351807+0.7623449818206758im,-0.05433818902992066+0.6511041246330369im,-0.06092919377429644+0.5238309230062319im,-0.06601991955554824+0.3836592655598524im,-0.06948501580979141+0.2340406437032877im,-0.07123916036725134+0.07865916446238384im,-0.07123916036725134-0.07865916446238384im,-0.06948501580979141-0.2340406437032877im,-0.06601991955554824-0.3836592655598524im,-0.06092919377429644-0.5238309230062319im,-0.05433818902992066-0.6511041246330369im,-0.04640919773351807-0.7623449818206758im,-0.03733745796561948-0.854814375477951im,-0.02734634606668882-0.9262354022446964im,-0.01668187637027696-0.9748494394091176im,-0.005606643513655412-0.9994594480354867im]
matlab_k = 3.748372513504540e-06

# Test that our answers are close to MATLAB's and at least as accurate
zpkfilter_eq(Chebyshev1(20, 1), ZPKFilter([], matlab_p, matlab_k))
zpkfilter_accuracy(Chebyshev1(20, 1), ZPKFilter([], matlab_p, matlab_k), Chebyshev1(BigFloat, 20, 1))

# Poles of 19 pole Butterworth filter prototype with 1 dB passband ripple from MATLAB 2013b:
#=
    [z, p, k] = cheb1ap(19, 1)
    sprintf('%.16g%+.16gim,', [real(p) imag(p)]')
    sprintf('%.16g', k)
=#
matlab_p = [-0.006212227114275604+0.9994004289494353im,-0.0184672279812169+0.9721393904901743im,-0.03021849100930104+0.9183609236365976im,-0.04114547237637834+0.8395319647744703im,-0.05095011251526942+0.7378027624097928im,-0.0593649664017747+0.6159482239948341im,-0.06616049875775722+0.4772922236864799im,-0.07115134517517649+0.3256169357239265im,-0.07420136837393095+0.1650596665748249im,-0.07522737167197566+0im,-0.07420136837393095-0.1650596665748249im,-0.07115134517517649-0.3256169357239265im,-0.06616049875775722-0.4772922236864799im,-0.0593649664017747-0.6159482239948341im,-0.05095011251526942-0.7378027624097928im,-0.04114547237637834-0.8395319647744703im,-0.03021849100930104-0.9183609236365976im,-0.0184672279812169-0.9721393904901743im,-0.006212227114275604-0.9994004289494353im]
matlab_k = 7.496745027009062e-06

# Test that our answers are close to MATLAB's and at least as accurate
zpkfilter_eq(Chebyshev1(19, 1), ZPKFilter([], matlab_p, matlab_k))
zpkfilter_accuracy(Chebyshev1(19, 1), ZPKFilter([], matlab_p, matlab_k), Chebyshev1(BigFloat, 19, 1))

#
# Chebyshev type II filter
#

# 20 pole Butterworth filter prototype with 1 dB passband ripple from MATLAB 2013b:
#=
    [z, p, k] = cheb2ap(20, 1)
    sprintf('%.16g%+.16gim,', [real(z) imag(z)]')
    sprintf('%.16g%+.16gim,', [real(p) imag(p)]')
    sprintf('%.16g', k)
=#
matlab_z = [0+1.003092198482826im,0-1.003092198482826im,0+1.028415193665209im,0-1.028415193665209im,0+1.082392200292394im,0-1.082392200292394im,0+1.172827696614009im,0-1.172827696614009im,0+1.315086999890785im,0-1.315086999890785im,0+1.539769043222366im,0-1.539769043222366im,0+1.913880855430943im,0-1.913880855430943im,0+2.613125929752753im,0-2.613125929752753im,0+4.283657569731185im,0-4.283657569731185im,0+12.74549484318238im,0-12.74549484318238im]
matlab_p = [-0.001929675700544753-1.002788598382701im,-0.006034873757627809-1.028072311003105im,-0.010957850142167-1.081957626858823im,-0.01756373888305481-1.172213900920771im,-0.02744258171352171-1.314120759169777im,-0.04403133250736632-1.538048178699878im,-0.07621942506983267-1.910267531994771im,-0.1536690174908018-2.6032737531287im,-0.4316587213861887-4.238414904224276im,-3.610083779585366-11.62012078608171im,-3.610083779585366+11.62012078608171im,-0.4316587213861887+4.238414904224276im,-0.1536690174908018+2.6032737531287im,-0.07621942506983267+1.910267531994771im,-0.04403133250736632+1.538048178699878im,-0.02744258171352171+1.314120759169777im,-0.01756373888305481+1.172213900920771im,-0.010957850142167+1.081957626858823im,-0.006034873757627809+1.028072311003105im,-0.001929675700544753+1.002788598382701im]
matlab_k = 0.8912509381337452
matlab_f = ZPKFilter(matlab_z, matlab_p, matlab_k)

# The gain shows a greater discrepancy with the gain of the enhanced
# precision filter than MATLAB's, but AFAICT our filter is still more
# accurate: the poles and zeros are each more accurate, and computing
# the gain in enhanced precision yields the same value.

f = Chebyshev2(20, 1)
zpkfilter_eq(f, matlab_f)
zpkfilter_accuracy(f, matlab_f, Chebyshev2(BigFloat, 20, 1), relerr=4)

# 19 pole Butterworth filter prototype with 1 dB passband ripple from MATLAB 2013b:
#=
    [z, p, k] = cheb2ap(19, 1)
    sprintf('%.16g%+.16gim,', [real(z) imag(z)]')
    sprintf('%.16g%+.16gim,', [real(p) imag(p)]')
    sprintf('%.16g', k)
=#
matlab_z = [0+1.003427212662145im,0-1.003427212662145im,0+1.031565634068626im,0-1.031565634068626im,0+1.091973276457601im,0-1.091973276457601im,0+1.194505544554793im,0-1.194505544554793im,0+1.359205519207709im,0-1.359205519207709im,0+1.628100459889458im,0-1.628100459889458im,0+2.101072544213108im,0-2.101072544213108im,0+3.079770972368366im,0-3.079770972368366im,0+6.075533820974263im,0-6.075533820974263im]
matlab_p = [-0.002139218568677465-1.003090263923361im,-0.006720708432694556-1.031180124433753im,-0.01232195186345456-1.091472452757043im,-0.02007306444044342-1.193772338103026im,-0.03217434756636293-1.357992987527419im,-0.05375960613144559-1.625783486413104im,-0.09966387580015423-2.09563676642751im,-0.2295218019875045-3.061543704346429im,-0.9149742173123645-5.932401750359961im,-38.84170469908301+0im,-0.9149742173123645+5.932401750359961im,-0.2295218019875045+3.061543704346429im,-0.09966387580015423+2.09563676642751im,-0.05375960613144559+1.625783486413104im,-0.03217434756636293+1.357992987527419im,-0.02007306444044342+1.193772338103026im,-0.01232195186345456+1.091472452757043im,-0.006720708432694556+1.031180124433753im,-0.002139218568677465+1.003090263923361im]
matlab_k = 37.33930783884512
matlab_f = ZPKFilter(matlab_z, matlab_p, matlab_k)

f = Chebyshev2(19, 1)
zpkfilter_eq(f, matlab_f)
zpkfilter_accuracy(f, matlab_f, Chebyshev2(BigFloat, 19, 1), relerr=2)

#
# Elliptic filter
#

# 20 pole elliptic filter prototype with 0.1 dB passband ripple and 10
# dB stopband ripple from MATLAB 2013b:
#=
    [z, p, k] = ellipap(20, 0.1, 10)
    sprintf('%.16g%+.16gim,', [real(z) imag(z)]')
    sprintf('%.16g%+.16gim,', [real(p) imag(p)]')
    sprintf('%.16g', k)
=#
matlab_z = [0-1.953252853757711im,0+1.953252853757711im,0-1.069599937693626im,0+1.069599937693626im,0-1.007032209276402im,0+1.007032209276402im,0-1.000730367252677im,0+1.000730367252677im,0-1.000076070678288im,0+1.000076070678288im,0-1.000007925871874im,0+1.000007925871874im,0-1.000000826312751im,0+1.000000826312751im,0-1.000000086632384im,0+1.000000086632384im,0-1.000000009576065im,0+1.000000009576065im,0-1.000000001633984im,0+1.000000001633984im]
matlab_p = [-0.6665552331000151-0.9586923119289487im,-0.6665552331000151+0.9586923119289487im,-0.06598422536270995-1.015686365686964im,-0.06598422536270995+1.015686365686964im,-0.006773291804375061-1.001821704927991im,-0.006773291804375061+1.001821704927991im,-0.0007045345890457818-1.000191780070487im,-0.0007045345890457818+1.000191780070487im,-7.339097059946433e-05-1.000020003037497im,-7.339097059946433e-05+1.000020003037497im,-7.646273200288023e-06-1.000002084836051im,-7.646273200288023e-06+1.000002084836051im,-7.966432964216136e-07-1.000000217756124im,-7.966432964216136e-07+1.000000217756124im,-8.299912436196738e-08-1.000000023227775im,-8.299912436196738e-08+1.000000023227775im,-8.637784231770464e-09-1.000000002962719im,-8.637784231770464e-09+1.000000002962719im,-8.070974080610083e-10-1.000000000874029im,-8.070974080610083e-10+1.000000000874029im]
matlab_k = 0.3162277662398871
matlab_f = ZPKFilter(matlab_z, matlab_p, matlab_k)

f = Elliptic(20, 0.1, 10)
zpkfilter_eq(f, matlab_f)
zpkfilter_accuracy(f, matlab_f, Elliptic(BigFloat, 20, 0.1, 10), eps=1e-9)

# 19 pole elliptic filter prototype with 0.1 dB passband ripple and 10
# dB stopband ripple from MATLAB 2013b:
#=
    [z, p, k] = ellipap(19, 0.1, 10)
    sprintf('%.16g%+.16gim,', [real(z) imag(z)]')
    sprintf('%.16g%+.16gim,', [real(p) imag(p)]')
    sprintf('%.16g', k)
=#
matlab_z = [0-1.232609672486912im,0+1.232609672486912im,0-1.021948255265862im,0+1.021948255265862im,0-1.002264470627064im,0+1.002264470627064im,0-1.000235691798532im,0+1.000235691798532im,0-1.000024555183434im,0+1.000024555183434im,0-1.000002559984621im,0+1.000002559984621im,0-1.000000268394018im,0+1.000000268394018im,0-1.00000002966741im,0+1.00000002966741im,0-1.000000005062212im,0+1.000000005062212im]
matlab_p = [-0.2102804655411482-1.034611529620449im,-0.2102804655411482+1.034611529620449im,-0.0210635609742992-1.005504639644632im,-0.0210635609742992+1.005504639644632im,-0.002183581313750937-1.00059265523267im,-0.002183581313750937+1.00059265523267im,-0.0002273805846473938-1.000061954788962im,-0.0002273805846473938+1.000061954788962im,-2.368886404726722e-05-1.000006458811389im,-2.368886404726722e-05+1.000006458811389im,-2.468065209676509e-06-1.00000067462382im,-2.468065209676509e-06+1.00000067462382im,-2.571378798158274e-07-1.000000071961468im,-2.571378798158274e-07+1.000000071961468im,-2.676054057433227e-08-1.000000009178736im,-2.676054057433227e-08+1.000000009178736im,-2.500451772444974e-09-1.00000000270781im,-2.500451772444974e-09+1.00000000270781im,-1.309071549907812+0im]
matlab_k = 0.9266824319626478
matlab_f = ZPKFilter(matlab_z, matlab_p, matlab_k)

f = Elliptic(19, 0.1, 10)
zpkfilter_eq(f, matlab_f)
zpkfilter_accuracy(f, matlab_f, Elliptic(BigFloat, 19, 0.1, 10), eps=4e-9)

#
# Conversion between zpk and tf
#

# Output of [z, p, k] = buttap(20); [b, a] = zp2tf(z, p, k)
m_a = [1,12.74549484318237,81.22381939879423,343.6513712403923,1081.352361133001,2687.409807920676,5468.931438945091,9326.061201886809,13528.36656744904,16852.27707949905,18122.54155403868,16852.27707949905,13528.36656744904,9326.061201886809,5468.931438945092,2687.409807920676,1081.352361133001,343.6513712403923,81.22381939879423,12.74549484318237,1]

f = convert(TFFilter, Butterworth(20))
@test_approx_eq coefb(f) [1]
@test_approx_eq coefa(f) m_a

# Test that our answers are more accurate than MATLAB's
accurate_a = coefa(convert(TFFilter, Butterworth(BigFloat, 20)))
@test_approx_eq coefa(f) float64(accurate_a)
@test sum(abs(coefa(f) - accurate_a)) <= sum(abs(m_a - accurate_a))

#
# Conversion between tf and zpk
#

f = convert(ZPKFilter, convert(TFFilter, Butterworth(20)))
zpkfilter_eq(f, butter, 1e-6)

# TODO: Make this more accurate!

# Output of [z, p, k] = buttap(20); [b, a] = zp2tf(z, p, k); tf2zpk(b, a)
m_p = [-0.07845909573254482+0.9969173337335029im,-0.07845909573254482-0.9969173337335029im,-0.2334453637958131+0.9723699203918822im,-0.2334453637958131-0.9723699203918822im,-0.3826834327796701+0.9238795325396184im,-0.3826834327796701-0.9238795325396184im,-0.5224985628488221+0.8526401643454914im,-0.5224985628488221-0.8526401643454914im,-0.6494480541398985+0.7604059651905597im,-0.6494480541398985-0.7604059651905597im,-0.760405952587916+0.6494480502272874im,-0.760405952587916-0.6494480502272874im,-0.8526401859847815+0.5224985598169277im,-0.8526401859847815-0.5224985598169277im,-0.9238795057196649+0.3826834415328767im,-0.9238795057196649-0.3826834415328767im,-0.9969173244841298+0.07845911266921719im,-0.9969173244841298-0.07845911266921719im,-0.97236994351794+0.2334453500964366im,-0.97236994351794-0.2334453500964366im]
# println(complex128([sort(f.p, lt=lt) - sort(Butterworth(BigFloat, 20).p, lt=lt) sort(m_p, lt=lt) - sort(Butterworth(BigFloat, 20).p, lt=lt)]))
zpkfilter_accuracy(f, ZPKFilter([], m_p, 1), Butterworth(BigFloat, 20); eps=1e-6, relerr=9)

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
accurate_ft = analogfilter(Lowpass(0.5), Butterworth(BigFloat, 20))
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
accurate_ft = analogfilter(Highpass(0.5), Butterworth(BigFloat, 20))
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
accurate_f10 = ZPKFilter(BigFloat[], accurate_prototype, BigFloat(1))
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

# Test that filters converted to SOS are equivalent
for ftype in (Lowpass(0.5), Highpass(0.5), Bandpass(0.25, 0.75), Bandstop(0.25, 0.75))
    for order in 1:4
        for proto in (Butterworth(order), Chebyshev1(order, 1), Chebyshev2(order, 1))
            f = digitalfilter(ftype, proto)
            @test_approx_eq filt(f, x) filt(convert(SOSFilter, f), x)
        end
    end
end

# Test that a numerically challenging filter (high order, clustered roots) has acceptable errors in its coefficients after conversion to SOS
f = ZPKFilter(ones(100), 0.99*ones(100), 1)
g = convert(SOSFilter, f)
tffilter_eq(convert(TFFilter, f), convert(TFFilter, g))

# Test designing filters as SOS
for ftype in (Lowpass(0.5), Highpass(0.5), Bandpass(0.25, 0.75), Bandstop(0.25, 0.75))
    for order in 1:(isa(ftype, Lowpass) || isa(ftype, Highpass) ? 4 : 2)
        for proto in (Butterworth(order), Chebyshev1(order, 1), Chebyshev2(order, 1))
            f1 = digitalfilter(ftype, proto)
            f2 = digitalfilter(ftype, convert(SOSFilter, proto))
            @assert isa(f2, SOSFilter)

            # Test out-of-place filtering
            @test_approx_eq filt(f1, x) filt(f2, x)
            # Test in-place filtering
            y = copy(x)
            @test_approx_eq filt(f1, x) filt!(y, f2, y)
        end
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
            try
                zpkfilter_eq(f, convert(ZPKFilter, f3), sqrt(eps()))
            catch e
                println("Conversion from $ftype1 to $ftype2 failed:")
                rethrow(e)
            end
        end
    end
end

for proto in (Butterworth(3), Chebyshev1(3, 1), Chebyshev2(3, 1))
    f = digitalfilter(Lowpass(0.5), proto)
    for ftype1 in (ZPKFilter, TFFilter, SOSFilter)
        f2 = convert(ftype1, f)
        for ftype2 in (ZPKFilter, TFFilter, SOSFilter)
            f3 = convert(ftype2, f2)
            zpkfilter_eq(f, convert(ZPKFilter, f3), 2e-5)
        end
    end
end
