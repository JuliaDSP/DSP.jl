!(dirname(@__FILE__) in LOAD_PATH) && push!(LOAD_PATH, dirname(@__FILE__))
using DSP, Test, FilterTestHelpers, Polynomials

@testset "convert to SOS" begin
    # Test conversion to SOS against MATLAB
    # Test poles/zeros generated with:
    #=
    mag = rand(10)
    arg = rand(10)
    z = mag.*cospi(arg) + im*mag.*sinpi(arg)
    z = [z, conj(z)]
    mag = rand(10)
    arg = rand(10)
    p = mag.*cospi(arg) + im*mag.*sinpi(arg)
    p = [p, conj(p)]
    k = real(prod(p))/real(prod(z))
    =#
    z = [0.07674942314081117 + 0.5605934468331276im,
        -0.10637764563083413 + 0.00938830970435945im,
        0.19723688182100613 + 0.20500254264692958im,
        0.07019769777809694 + 0.00040337356287483566im,
        -0.029225401438832663 + 0.35500551559734517im,
        -0.35980697033738923 + 0.21595798755003964im,
        0.3712083544916094 + 0.41767928564585416im,
        -0.17585090520154495 + 0.32300027988596314im,
        -0.2308322620393009 + 0.3539248310694154im,
        -0.008498685920569043 + 0.028356975487226484im,
        0.07674942314081117 - 0.5605934468331276im,
        -0.10637764563083413 - 0.00938830970435945im,
        0.19723688182100613 - 0.20500254264692958im,
        0.07019769777809694 - 0.00040337356287483566im,
        -0.029225401438832663 - 0.35500551559734517im,
        -0.35980697033738923 - 0.21595798755003964im,
        0.3712083544916094 - 0.41767928564585416im,
        -0.17585090520154495 - 0.32300027988596314im,
        -0.2308322620393009 - 0.3539248310694154im,
        -0.008498685920569043 - 0.028356975487226484im]
    p = [-0.946179900464128 + 0.23187351222922784im,
        0.05191136844411161 + 0.8713748123758278im,
        -0.05109307872385063 + 0.7440280322767342im,
        -0.032892467440199066 + 0.024218620496805687im,
        -0.1266287653888985 + 0.33150956246323654im,
        -0.15051989876024102 + 0.9373469058653078im,
        -0.6672740143157547 + 0.5034138963468052im,
        -0.790336466064852 + 0.109069102200402im,
        -0.009317017740249631 + 0.013158400271660778im,
        0.5075897927503011 + 0.02902816288546107im,
        -0.946179900464128 - 0.23187351222922784im,
        0.05191136844411161 - 0.8713748123758278im,
        -0.05109307872385063 - 0.7440280322767342im,
        -0.032892467440199066 - 0.024218620496805687im,
        -0.1266287653888985 - 0.33150956246323654im,
        -0.15051989876024102 - 0.9373469058653078im,
        -0.6672740143157547 - 0.5034138963468052im,
        -0.790336466064852 - 0.109069102200402im,
        -0.009317017740249631 - 0.013158400271660778im,
        0.5075897927503011 - 0.02902816288546107im]
    k = 10586.686805294861

    # MATLAB:
    #=
    [sos, g] = zp2sos(z, p, k)
    save('...', 'sos', '-ascii', '-double')
    =#
    m_sos_full = [
    1.0000000000000000e+00  -1.4039539555619387e-01   4.9278794835762620e-03   1.0000000000000000e+00   1.8634035480499262e-02   2.5995031728136877e-04
    1.0000000000000000e+00   1.6997371841138087e-02   8.7634572115964215e-04   1.0000000000000000e+00   6.5784934880398133e-02   1.6684559930728524e-03
    1.0000000000000000e+00   2.1275529126166826e-01   1.1404343849064294e-02   1.0000000000000000e+00   2.5325753077779700e-01   1.2593343422848324e-01
    1.0000000000000000e+00  -3.9447376364201225e-01   8.0928430042179728e-02   1.0000000000000000e+00  -1.0151795855006023e+00   2.5849003194479853e-01
    1.0000000000000000e+00  -7.4241670898321876e-01   3.1225162810199936e-01   1.0000000000000000e+00   1.0218615744770126e-01   5.5618821550707065e-01
    1.0000000000000000e+00   3.5170181040308990e-01   1.3525272166661328e-01   1.0000000000000000e+00   1.5806729321297039e+00   6.3652779864668074e-01
    1.0000000000000000e+00   4.6166452407860181e-01   1.7854631924569470e-01   1.0000000000000000e+00   1.3345480286315095e+00   6.9868016121613397e-01
    1.0000000000000000e+00   5.8450802877665325e-02   1.2688304019379779e-01   1.0000000000000000e+00  -1.0382273688822322e-01   7.6198885381674941e-01
    1.0000000000000000e+00  -1.5349884628162233e-01   3.2015548658469395e-01   1.0000000000000000e+00   3.0103979752048204e-01   9.0127546185805940e-01
    1.0000000000000000e+00   7.1961394067477846e-01   1.7609890829003397e-01   1.0000000000000000e+00   1.8923598009282561e+00   9.4902172971582510e-01
    ]
    @test m_sos_full ≈ sosfilter_to_matrix(convert(SecondOrderSections, ZeroPoleGain(z, p, k)))

    # And with half of the zeros removed
    # MATLAB:
    # [sos, g] = zp2sos(z([1:5, 11:15]), p, k)
    zp = z[[1:5; 11:15]]
    m_sos_half = [
    0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   1.0000000000000000e+00   1.8634035480499262e-02   2.5995031728136877e-04
    0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   1.0000000000000000e+00   6.5784934880398133e-02   1.6684559930728524e-03
    0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   1.0000000000000000e+00   2.5325753077779700e-01   1.2593343422848324e-01
    0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   1.0000000000000000e+00  -1.0151795855006023e+00   2.5849003194479853e-01
    0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   1.0000000000000000e+00   1.0218615744770126e-01   5.5618821550707065e-01
    1.0000000000000000e+00  -3.9447376364201225e-01   8.0928430042179728e-02   1.0000000000000000e+00   1.5806729321297039e+00   6.3652779864668074e-01
    1.0000000000000000e+00  -1.4039539555619387e-01   4.9278794835762620e-03   1.0000000000000000e+00   1.3345480286315095e+00   6.9868016121613397e-01
    1.0000000000000000e+00   5.8450802877665325e-02   1.2688304019379779e-01   1.0000000000000000e+00  -1.0382273688822322e-01   7.6198885381674941e-01
    1.0000000000000000e+00  -1.5349884628162233e-01   3.2015548658469395e-01   1.0000000000000000e+00   3.0103979752048204e-01   9.0127546185805940e-01
    1.0000000000000000e+00   2.1275529126166826e-01   1.1404343849064294e-02   1.0000000000000000e+00   1.8923598009282561e+00   9.4902172971582510e-01
    ]
    @test m_sos_half ≈ sosfilter_to_matrix(convert(SecondOrderSections, ZeroPoleGain(zp, p, k)))

    # And with an extra real pole
    pp = [p; 0.7]
    # MATLAB:
    #=
    pp = [p'; 0.7]
    [sos, g] = zp2sos(z, pp, k)
    =#
    m_sos_extra = [
    0.0000000000000000e+00   1.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00  -6.9999999999999996e-01   0.0000000000000000e+00
    1.0000000000000000e+00  -1.4039539555619387e-01   4.9278794835762620e-03   1.0000000000000000e+00   1.8634035480499262e-02   2.5995031728136877e-04
    1.0000000000000000e+00   1.6997371841138087e-02   8.7634572115964215e-04   1.0000000000000000e+00   6.5784934880398133e-02   1.6684559930728524e-03
    1.0000000000000000e+00   2.1275529126166826e-01   1.1404343849064294e-02   1.0000000000000000e+00   2.5325753077779700e-01   1.2593343422848324e-01
    1.0000000000000000e+00  -3.9447376364201225e-01   8.0928430042179728e-02   1.0000000000000000e+00  -1.0151795855006023e+00   2.5849003194479853e-01
    1.0000000000000000e+00  -7.4241670898321876e-01   3.1225162810199936e-01   1.0000000000000000e+00   1.0218615744770126e-01   5.5618821550707065e-01
    1.0000000000000000e+00   3.5170181040308990e-01   1.3525272166661328e-01   1.0000000000000000e+00   1.5806729321297039e+00   6.3652779864668074e-01
    1.0000000000000000e+00   4.6166452407860181e-01   1.7854631924569470e-01   1.0000000000000000e+00   1.3345480286315095e+00   6.9868016121613397e-01
    1.0000000000000000e+00   5.8450802877665325e-02   1.2688304019379779e-01   1.0000000000000000e+00  -1.0382273688822322e-01   7.6198885381674941e-01
    1.0000000000000000e+00  -1.5349884628162233e-01   3.2015548658469395e-01   1.0000000000000000e+00   3.0103979752048204e-01   9.0127546185805940e-01
    1.0000000000000000e+00   7.1961394067477846e-01   1.7609890829003397e-01   1.0000000000000000e+00   1.8923598009282561e+00   9.4902172971582510e-01
    ]
    @test m_sos_extra ≈ sosfilter_to_matrix(convert(SecondOrderSections, ZeroPoleGain(z, pp, k)))

    # And with only poles (no zeros)
    m_sos_only_poles = copy(m_sos_full)
    m_sos_only_poles[:, 1:2] .= 0
    m_sos_only_poles[:, 3] .= 1
    @test m_sos_only_poles ≈ sosfilter_to_matrix(convert(SecondOrderSections, ZeroPoleGain(Float64[], p, k)))

    # Test that a filter with repeated zeros is handled properly
    # MATLAB:
    #=
    [z,p,k] = butter(2, [49.5 50.5]/500, 'stop')
    [sos,g] = zp2sos(z, p, k)
    =#
    m_sos_butterworth_bs = [
    1.0000000000000000e+00  -1.9021224191804869e+00   1.0000000000000000e+00   1.0000000000000000e+00  -1.8964983429993663e+00   9.9553672990017417e-01
    1.0000000000000000e+00  -1.9021224191804869e+00   1.0000000000000000e+00   1.0000000000000000e+00  -1.8992956433548462e+00   9.9559721515078736e-01
    ]
    f = convert(SecondOrderSections, digitalfilter(Bandstop(49.5, 50.5; fs=1000), DSP.Butterworth(2)))
    @test m_sos_butterworth_bs ≈ sosfilter_to_matrix(f)
    @test f.g ≈ 0.995566972017647

    # Test that a numerically challenging filter (high order, clustered
    # roots) has acceptable errors in its coefficients after conversion to
    # SOS
    f = ZeroPoleGain(ones(100), 0.99*ones(100), 1)
    g = convert(SecondOrderSections, f)
    tffilter_eq(convert(PolynomialRatio, f), convert(PolynomialRatio, g))
end

@testset "conversions" begin
    for f in (digitalfilter(Lowpass(0.5), Butterworth(1)), digitalfilter(Lowpass(0.5), Butterworth(2)),
              digitalfilter(Bandpass(0.25, 0.75), Butterworth(1)))
        for ftype1 in (ZeroPoleGain, PolynomialRatio, Biquad, SecondOrderSections)
            f2 = convert(ftype1, f)
            for ftype2 in (ZeroPoleGain, PolynomialRatio, Biquad, SecondOrderSections)
                f3 = convert(ftype2, f)
                try
                    zpkfilter_eq(convert(ZeroPoleGain, f), convert(ZeroPoleGain, f3), sqrt(eps()))
                catch e
                    println("Conversion from $ftype1 to $ftype2 failed:")
                    rethrow(e)
                end
            end
        end
    end

    for proto in (Butterworth(3), Chebyshev1(3, 1), Chebyshev2(3, 1))
        f = digitalfilter(Lowpass(0.5), proto)
        for ftype1 in (ZeroPoleGain, PolynomialRatio, SecondOrderSections)
            f2 = convert(ftype1, f)
            for ftype2 in (ZeroPoleGain, PolynomialRatio, SecondOrderSections)
                f3 = convert(ftype2, f2)
                zpkfilter_eq(convert(ZeroPoleGain, f), convert(ZeroPoleGain, f3), 2e-5)
            end
        end
    end
end

@testset "filter gain and composition" begin
    # Test setting filter gain
    x = randn(100)
    f1 = digitalfilter(Lowpass(0.3), Butterworth(2))
    y = filt(f1, x)
    for ty in (ZPKFilter, PolynomialRatio, Biquad, SecondOrderSections)
        @test filt(3*convert(ty, f1), x) ≈ 3*y
        @test filt(convert(ty, f1)*3, x) ≈ 3*y
    end

    # Test composing filters
    f2 = digitalfilter(Highpass(0.5), Butterworth(1))
    y = filt(f2, y)
    for ty in (ZPKFilter, PolynomialRatio, Biquad, SecondOrderSections)
        @test filt(convert(ty, f1)*convert(ty, f2), x) ≈ y
    end

    f3 = digitalfilter(Bandstop(0.35, 0.4), Butterworth(1))
    y = filt(f3, y)
    for ty in (ZPKFilter, PolynomialRatio, Biquad, SecondOrderSections)
        @test filt(convert(ty, f1)*convert(ty, f2)*convert(ty, f3), x) ≈ y
    end
    @test filt(convert(Biquad, f1)*(convert(Biquad, f2)*convert(Biquad, f3)), x) ≈ y
    @test filt((convert(Biquad, f1)*convert(Biquad, f2))*convert(Biquad, f3), x) ≈ y
end

@testset "misc" begin
    f1 = digitalfilter(Lowpass(0.3), Butterworth(2))
    f2 = digitalfilter(Highpass(0.5), Butterworth(1))
    # Test some otherwise untested code paths
    @test promote_type(ZeroPoleGain{ComplexF32,ComplexF64,Float32}, ZeroPoleGain{ComplexF64,ComplexF32,Float64}) == ZeroPoleGain{ComplexF64,ComplexF64,Float64}
    @test convert(ZeroPoleGain{ComplexF64,ComplexF64,Float64}, f1) === f1
    f1f = convert(ZeroPoleGain{ComplexF32,ComplexF32,Float32}, f1)
    @test f1f.z == convert(Vector{ComplexF32}, f1.z)
    @test f1f.p == convert(Vector{ComplexF32}, f1.p)
    @test f1f.k == convert(Float32, f1.k)

    @test_throws ArgumentError PolynomialRatio(Float64[], Float64[])
    @test promote_type(PolynomialRatio{Float32}, PolynomialRatio{Int}) == PolynomialRatio{Float32}
    f1f = convert(PolynomialRatio{Float32}, f1)
    f1p = convert(PolynomialRatio, f1)
    @test f1f.b == convert(Poly{Float32}, f1p.b)
    @test f1f.a == convert(Poly{Float32}, f1p.a)

    b = Biquad(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 2)
    @test b.b0 === 0.0
    @test b.b1 === 2*1/3
    @test b.b2 === 2*2/3
    @test b.a1 === 4/3
    @test b.a2 === 5/3
    @test promote_type(Biquad{Float32}, Biquad{Int}) == Biquad{Float32}
    @test convert(Biquad{Float32}, b) == Biquad{Float32}(0.0, 2*1/3, 2*2/3, 4/3, 5/3)
    zpkfilter_eq(convert(ZeroPoleGain{ComplexF64,ComplexF64,Float64}, b), convert(ZeroPoleGain, b), 0.0)
    f = convert(PolynomialRatio, Biquad(2.0, 0.0, 0.0, 0.0, 0.0))
    @test coefb(f) == [2.0]
    @test coefa(f) == [1.0]
    @test convert(Biquad, PolynomialRatio([4.0], [2.0])) == Biquad(2.0, 0.0, 0.0, 0.0, 0.0)
    @test Biquad(2.0, 0.0, 0.0, 0.0, 0.0)*2 == Biquad(4.0, 0.0, 0.0, 0.0, 0.0)
    @test convert(Biquad{Float64}, f1) == convert(Biquad, f1)
    f = PolynomialRatio(Float64[1.0], Float64[1.0])
    empty!(f.b.a)
    empty!(f.a.a)
    @test_throws ArgumentError convert(Biquad, f)
    @test_throws ArgumentError convert(SOSFilter, ZPKFilter([0.5 + 0.5im, 0.5 + 0.5im], [0.5 + 0.5im, 0.5 - 0.5im], 1))
    @test_throws ArgumentError convert(SOSFilter, ZPKFilter([0.5 + 0.5im, 0.5 - 0.5im], [0.5 + 0.5im, 0.5 + 0.5im], 1))

    @test promote_type(SecondOrderSections{Float64,Float32}, SecondOrderSections{Float32,Float64}) == SecondOrderSections{Float64,Float64}
    @test convert(SecondOrderSections{Float32,Float32}, convert(SecondOrderSections, b)).biquads == convert(SecondOrderSections, convert(Biquad{Float32}, b)).biquads
    @test convert(Biquad, convert(SecondOrderSections, b)) == b
    @test_throws ArgumentError convert(Biquad, convert(SecondOrderSections, f1*f2))
    f1p = convert(PolynomialRatio{Float32}, convert(PolynomialRatio, convert(SecondOrderSections, f1)))
    f1f = convert(PolynomialRatio{Float32}, convert(SecondOrderSections, f1))
    @test f1p.b == f1f.b
    @test f1p.a == f1f.a
end

@testset "coef() testing" begin
    # relies on conversion for non-polynomial ratio filter coeff objects
    # so no need to test conversion correctness

    @testset "Biquad" begin
        bs = [1, 2, 3, 4, 5]
        B = Biquad(bs...)
        @test coefa(B) == [1, 4, 5]
        @test coefb(B) == [1, 2, 3]

        bs = Float64[20, 16, 13, 31, 33]
        B = Biquad(bs...)
        @test coefa(B) == [1.0, 31, 33]
        @test coefb(B) == [20.0, 16, 13]
    end

    @testset "SecondOrderSections" begin
        bs = [2.0, 0, 0, 0, 0]
        B = SecondOrderSections(repeat([Biquad(bs...)], 2), 0.25)
        @test coefb(B) == [1.0]
        @test coefa(B) == [1.0]

        bs = [0, 1, 0, 0, 0]
        B = SecondOrderSections(repeat([Biquad(bs...)], 2), 1)
        @test coefb(B) == [1]
        @test coefa(B) == [1, 0, 0]
    end

    @testset "ZeroPoleGain" begin
        f = ZeroPoleGain([0], [-1, 1], 1)
        @test coefa(f) == [1, 0, -1]
        @test coefb(f) == [1, 0]

        f = ZeroPoleGain(Int[], [-0.25, 0.25], 1)
        @test coefa(f) == [1.0, 0, -1/16]
        @test coefb(f) == [1.0]
    end
end
