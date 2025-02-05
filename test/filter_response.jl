using DSP, Test
using FilterTestHelpers: read_reference_data
using Polynomials: Polynomial


#######################################
#
#  http://www.mathworks.com.au/help/signal/ref/freqz.html
#  Example 1: Frequency response from transfer function
#
#  dlmwrite('freqz-eg1.txt',[w, abs(h)], 'delimiter', '\t', 'precision', '%.12f')
#
#######################################
@testset "freqresp from TF" begin
    # Matlab
    freqz_eg1_w_abs = read_reference_data("freqz-eg1.txt")
    matlab_abs = freqz_eg1_w_abs[:, 2]

    # Julia
    b0 = 0.05634
    b1 = [1,  1]
    b2 = [1, -1.0166, 1]
    a1 = [1, -0.683]
    a2 = [1, -1.4461, 0.7957]
    b = b0*conv(b1, b2)
    a = conv(a1, a2)

    w     = range(0, stop=6.280045284537, length=2001)
    h     = freqresp(PolynomialRatio(b, a), w)
    h_abs = convert(Array{Float64}, abs.(h))
    h_ref = [b0 * Polynomial(b1)(z) * Polynomial(b2)(z) / (Polynomial(a1)(z) * Polynomial(a2)(z)) for z ∈ exp.(-im*w)]

    # Test
    @test h_abs ≈ matlab_abs
    @test h ≈ h_ref
end

#######################################
#
# freqresp with conversion to PolynomialRatio may lead to undesirable roundoff errors
# check by using the same poles and zeros such that they should cancel
#
#######################################
@testset "freqresp pole/zero cancellation" begin
    rs = [1-10.0^-n for n in 1:15]
    @test freqresp(ZeroPoleGain(rs, reverse(rs), 42.0), range(0, stop=2π, length=50)) ≈ fill(42., 50)

    biq = Biquad(42.0, 84.0*real(.999999*exp(im)), 42.0*.999999^2, 2.0*real(.999999*exp(im)), .999999^2)
    @test freqresp(biq, range(0, stop=2π, length=50)) ≈ fill(42.0, 50)

    pr = [1-10.0^-n for n in 1:10]
    zr = reverse(pr)
    bs = [Biquad(42.0, 84.0*real(r[1]*exp(im)), 42.0*r[1]^2, 2.0*real(r[2]*exp(im)), r[2]^2) for r in zip(zr, pr)]
    sos = SecondOrderSections(bs, 42.0^-9)
    @test freqresp(sos, range(0, stop=2π, length=50)) ≈ fill(42.0, 50)

    @test freqresp(ZeroPoleGain(ComplexF64[], ComplexF64[], 42.0), range(0, stop=2π, length=50)) == fill(42.0, 50)
    @test freqresp(SecondOrderSections(Biquad{:z,Float64}[], 42.0), range(0, stop=2π, length=50)) == fill(42.0, 50)
end

#######################################
#
#  Test frequency, phase, impulse and step response
#
#  Data from Matlab using b,a and a from above:
#  ir = impz(b,a, 512);
#  stepr = stepz(b,a, 512);
#  h = abs(freqz(b,a));
#  [phi,f] =  phasez(b, a);
#  all = [f, ir, stepr, h, phi];
#  dlmwrite('responses-eg1.txt',all, 'delimiter', '\t', 'precision', '%.12f')
#
#######################################
@testset "impresp, stepresp, freqresp, phaseresp" begin
    Hdelay = PolynomialRatio{:z}([0, 1], [1])
    W = range(0, stop=2π, length=100)
    for Tf in (PolynomialRatio, ZeroPoleGain, Biquad, SecondOrderSections)
        H = convert(Tf, Hdelay)
        @test freqresp(H, W) ≈ exp.(-W*im)
        @test phaseresp(H, W) ≈ -W
        @test grpdelay(H, W) ≈ fill(1.0, length(W))
        @test impresp(H, 10) == [0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        @test stepresp(H, 10) == [0, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    end

    matlab_resp = read_reference_data("responses-eg1.txt")
    b0 = 0.05634
    b1 = [1,  1]
    b2 = [1, -1.0166, 1]
    a1 = [1, -0.683]
    a2 = [1, -1.4461, 0.7957]
    b = b0*conv(b1, b2)
    a = conv(a1, a2)
    df = PolynomialRatio(b, a)
    w = matlab_resp[:,1]

    #Impulse response
    impz_matlab = matlab_resp[:,2]
    @test impresp(df, 512) ≈ impz_matlab

    #Step response
    stepz_matlab = matlab_resp[:,3]
    @test stepresp(df, 512) ≈ stepz_matlab

    h_matlab = matlab_resp[:,4]
    @test abs.(freqresp(df, w)) ≈ h_matlab
    @test abs.(freqresp(SecondOrderSections(df), w)) ≈ h_matlab
    @test abs.(freqresp(ZeroPoleGain(df), w)) ≈ h_matlab

    phi_matlab = matlab_resp[:,5]
    @test phaseresp(df, w) ≈ phi_matlab


    # Test diffent versions of the functions
    H, w = freqresp(df)
    @test H == freqresp(df, w)
    @test minimum(w) ≤ 0                 # w should cover frequencies from 0 ...
    @test maximum(w) ≥ float(π)          # ... to π ...
    @test maximum(abs.(diff(w))) < 0.03  # and be sufficiently dense

    phi, w = phaseresp(df)
    @test phi == phaseresp(df, w)
    @test minimum(w) ≤ 0                 # w should cover frequencies from 0 ...
    @test maximum(w) ≥ float(π)          # ... to π ...
    @test maximum(abs.(diff(w))) < 0.03  # and be sufficiently dense

    @test cumsum(impresp(df)) ≈ stepresp(df)
end

#######################################
#
#  http://www.mathworks.com.au/help/signal/ref/freqs.html
#  Example 1: Frequency response from the transfer function
#
#  dlmwrite('freqs-eg1.txt',[w; mag; phasedeg]', 'delimiter', '\t', 'precision', '%.10f')
#
#######################################
@testset "freqresp(::FilterCoefficients{:s})" begin
    w = 10 .^ range(-1, stop=1, length=50)
    for (b, a, H) in (
            ([1, 0], [1], w*im), # s
            ([1], [1, 0], 1 ./ (w*im)), # 1/s
            ([0, 0, 1], [2, 10], 1 ./ (2*w*im .+ 10)), # 1/(2s + 10)
            ([1, 0, 1], [2, 10], (1 .- w.^2) ./ (2*w*im .+ 10)), # (s^2 + 1)/(2s + 10)
        )
        @test freqresp(PolynomialRatio{:s}(b, a), w) ≈ H
        @test freqresp(ZeroPoleGain(PolynomialRatio{:s}(b, a)), w) ≈ H
        if isone(a[1]) && lastindex(Polynomial(reverse(b))) <= lastindex(Polynomial(reverse(a)))
            @test freqresp(Biquad(PolynomialRatio{:s}(b, a)), w) ≈ H
        else
            # cannot represent this PolynomialRatio as a Biquad due to implied a0=1
            @test_broken freqresp(Biquad(PolynomialRatio{:s}(b, a)), w) ≈ H
        end
        if lastindex(Polynomial(reverse(b))) <= lastindex(Polynomial(reverse(a)))
            @test freqresp(SecondOrderSections(PolynomialRatio{:s}(b, a)), w) ≈ H
        else
            # cannot represent this PolynomialRatio as a Biquad due to implied a0=1
            @test_broken freqresp(SecondOrderSections(PolynomialRatio{:s}(b, a)), w) ≈ H
        end
    end

    # Julia
    a = [1.0, 0.4, 1.0]
    b = [0.2, 0.3, 1.0]
    w = 10 .^ range(-1, stop=1, length=50)

    h        = freqresp(PolynomialRatio{:s}(b, a), w)
    mag      = convert(Array{Float64}, abs.(h))
    phasedeg = (180/pi)*phaseresp(PolynomialRatio{:s}(b, a), w)

    # Matlab
    freqs_eg1_w_mag_phasedeg = read_reference_data("freqs-eg1.txt")
    matlab_w        = freqs_eg1_w_mag_phasedeg[:,1]
    matlab_mag      = freqs_eg1_w_mag_phasedeg[:,2]
    matlab_phasedeg = freqs_eg1_w_mag_phasedeg[:,3]
    h_ref = [Polynomial(reverse(b))(s) / Polynomial(reverse(a))(s) for s ∈ im*w]

    # Test
    @test w ≈ matlab_w
    @test h ≈ h_ref
    @test mag ≈ matlab_mag
    @test phasedeg ≈ matlab_phasedeg

    @test h ≈ freqresp(ZeroPoleGain(PolynomialRatio{:s}(b, a)), w)
    @test h ≈ freqresp(SecondOrderSections(PolynomialRatio{:s}(b, a)), w)
end

# ######################################
#
#  Test grpdelay
#
#  Data from Matlab using b,a and a from above:
#  [gd, w] = grpdelay(b, a, 512)
#  all = [w gd]
#  dlmwrite('grpdelay_eg1.txt', all, 'delimiter', '\t', 'precision', '%.12f')
#
# ######################################
@testset "grpdelay" begin
    matlab_delay = read_reference_data("grpdelay_eg1.txt")
    b0 = 0.05634
    b1 = [1,  1]
    b2 = [1, -1.0166, 1]
    a1 = [1, -0.683]
    a2 = [1, -1.4461, 0.7957]
    b = b0*conv(b1, b2)
    a = conv(a1, a2)
    df = PolynomialRatio(b, a)
    w = matlab_delay[:, 1]

    grpdelay_matlab = matlab_delay[:, 2]
    @test grpdelay(df, w) ≈ grpdelay_matlab

    # Test with IIR filters types I-IV
    tau, w = grpdelay(PolynomialRatio([1, 1, 1, 1, 1], [1]))
    @test tau ≈ fill(2.0, length(w))
    tau, w = grpdelay(PolynomialRatio([1, 1, 1, 1, 1, 1], [1]))
    @test tau ≈ fill(2.5, length(w))
    tau, w = grpdelay(PolynomialRatio([1, 0, -1], [1]))
    @test tau ≈ fill(1.0, length(w))
    tau, w = grpdelay(PolynomialRatio([1, -1], [1]))
    @test tau ≈ fill(0.5, length(w))
    @test minimum(w) ≤ 0                 # w should cover frequencies from 0 ...
    @test maximum(w) ≥ float(π)          # ... to π ...
    @test maximum(abs.(diff(w))) < 0.03  # and be sufficiently dense

    # group delay of analog filters
    w = matlab_delay[:, 1]
    @test grpdelay(PolynomialRatio{:s}([1, 2], [123]), w) ≈ 0.5 ./ (1 .+ w.^2/4)
    @test grpdelay(PolynomialRatio{:s}([123], [1, 2]), w) ≈ -0.5 ./ (1 .+ w.^2/4)
    @test grpdelay(PolynomialRatio{:s}([1, 2], [3, 9]), w) ≈ 0.5 ./ (1 .+ w.^2/4) - 1/3 ./ (1 .+ w.^2/9)
    tau, w = grpdelay(PolynomialRatio{:s}([1, 2], [123]))
    @test tau ≈ 0.5 ./ (1 .+ w.^2/4)
    @test minimum(w) < 2 < maximum(w) # position of zero covered
    tau, w = grpdelay(PolynomialRatio{:s}([123], [1, 2]))
    @test tau ≈ -0.5 ./ (1 .+ w.^2/4)
    @test minimum(w) < 2 < maximum(w) # position of pole covered
    tau, w = grpdelay(PolynomialRatio{:s}([1, 2], [3, 9]))
    @test tau ≈ 0.5 ./ (1 .+ w.^2/4) - 1/3 ./ (1 .+ w.^2/9)
    @test minimum(w) < 2 && 3 < maximum(w) # position of pole and zero covered
end
