using DSP, Compat, Compat.Test


#######################################
#
#  http://www.mathworks.com.au/help/signal/ref/freqz.html
#  Example 1: Frequency response from transfer function
#
#  dlmwrite('freqz-eg1.txt',[w, abs(h)], 'delimiter', '\t', 'precision', '%.12f')
#
#######################################
@testset "freqz fro TF" begin
    # Matlab
    freqz_eg1_w_abs = readdlm(joinpath(dirname(@__FILE__), "data", "freqz-eg1.txt"),'\t')
    matlab_w     = freqz_eg1_w_abs[:,1]
    matlab_abs   = freqz_eg1_w_abs[:,2]

    # Julia
    b0 = 0.05634
    b1 = [1,  1]
    b2 = [1, -1.0166, 1]
    a1 = [1, -0.683]
    a2 = [1, -1.4461, 0.7957]
    b = b0*conv(b1, b2)
    a = conv(a1, a2)

    #=w     = range(0, stop=2*pi, length=200)=#        # Does not produce same values as matlab
    h     = freqz(PolynomialRatio(b, a), matlab_w)   # So use frequencies from matlab
    h_abs = convert(Array{Float64}, abs.(h))

    # Test
    @test h_abs ≈ matlab_abs

    #=using Winston=#
    #=figure = plot(matlab_w/pi, 20*log10(h_abs))=#
    #=figure = oplot(matlab_w/pi, 20*log10(matlab_abs), "r--")=#
    #=ylim(-100, 20)=#
    #=ylabel("Magnitude (dB)")=#
    #=xlabel("Normalised Frequency (x pi rad/s)")=#
    #=file(figure, "MATLAB-freqz.png", width=1200, height=800)=#
end

#######################################
#
# freqz with conversion to PolynomialRatio may lead to undesirable roundoff errors
# check by using the same poles and zeros such that they should cancel
#
#######################################
@testset "freqz ple/zero cancellation" begin
    rs = [1-10.0^-n for n in 1:15]
    @test freqz(ZeroPoleGain(rs, reverse(rs), 42.0), Compat.range(0, stop=2π, length=50)) ≈ fill(42., 50)

    biq = Biquad(42.0, 84.0*real(.999999*exp(im)), 42.0*.999999^2, 2.0*real(.999999*exp(im)), .999999^2)
    @test freqz(biq, Compat.range(0, stop=2π, length=50)) ≈ fill(42.0, 50)

    pr = [1-10.0^-n for n in 1:10]
    zr = reverse(pr)
    bs = [Biquad(42.0, 84.0*real(r[1]*exp(im)), 42.0*r[1]^2, 2.0*real(r[2]*exp(im)), r[2]^2) for r in zip(zr, pr)]
    sos = SecondOrderSections(bs, 42.0^-9)
    @test freqz(sos, Compat.range(0, stop=2π, length=50)) ≈ fill(42.0, 50)

    @test freqz(ZeroPoleGain(ComplexF64[], ComplexF64[], 42.0), Compat.range(0, stop=2π, length=50)) == fill(42.0, 50)
    @test freqz(SecondOrderSections(Biquad{Float64}[], 42.0), Compat.range(0, stop=2π, length=50)) == fill(42.0, 50)
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
@testset "impz, stepz, freqz, phasez" begin
    matlab_resp = readdlm(joinpath(dirname(@__FILE__), "data", "responses-eg1.txt"),'\t')
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
    @test impz(df, 512) ≈ impz_matlab

    #Step response
    stepz_matlab = matlab_resp[:,3]
    @test stepz(df, 512) ≈ stepz_matlab

    h_matlab = matlab_resp[:,4]
    @test abs.(freqz(df, w)) ≈ h_matlab
    @test abs.(freqz(SecondOrderSections(df), w)) ≈ h_matlab
    @test abs.(freqz(ZeroPoleGain(df), w)) ≈ h_matlab

    phi_matlab = matlab_resp[:,5]
    @test phasez(df, w) ≈ phi_matlab


    # Test diffent versions of the functions
    @test freqz(df) == freqz(df, Compat.range(0, stop=pi, length=250))
    @test phasez(df) == phasez(df, Compat.range(0, stop=pi, length=250))
    @test cumsum(impz(df)) ≈ stepz(df)
end

#######################################
#
#  Test digital filter with frequency specified in hz
#
#  TODO: Create a MATLAB ground truth to compare against
#        Currently this just checks it runs without error
#
#######################################
@testset "freqz with fs" begin
    # Julia
    b0 = 0.05634
    b1 = [1,  1]
    b2 = [1, -1.0166, 1]
    a1 = [1, -0.683]
    a2 = [1, -1.4461, 0.7957]
    b = b0*conv(b1, b2)
    a = conv(a1, a2)

    fs    = 8192
    hz    = Compat.range(0, stop=fs, length=200)
    h     = freqz(PolynomialRatio(b, a), hz, fs)
    h_abs = convert(Array{Float64}, abs.(h))

    #=using Winston=#
    #=figure = plot(hz, 20*log10(h_abs))=#
    #=ylim(-100, 20)=#
    #=ylabel("Magnitude (dB)")=#
    #=xlabel("Frequency (Hz)")=#
    #=file(figure, "MATLAB-freqz-hz.png", width=1200, height=800)=#
end

#######################################
#
#  http://www.mathworks.com.au/help/signal/ref/freqs.html
#  Example 1: Frequency response from the transfer function
#
#  dlmwrite('freqs-eg1.txt',[w; mag; phasedeg]', 'delimiter', '\t', 'precision', '%.10f')
#
#######################################
@testset "freqs" begin
    # Julia
    a = [1.0, 0.4, 1.0]
    b = [0.2, 0.3, 1.0]
    w = 10 .^ Compat.range(-1, stop=1, length=50)

    h        = freqs(PolynomialRatio(b, a), w)
    mag      = convert(Array{Float64}, abs.(h))
    phasedeg = (180/pi)*convert(Array{Float64}, angle.(h))

    # Matlab
    freqs_eg1_w_mag_phasedeg = readdlm(joinpath(dirname(@__FILE__), "data", "freqs-eg1.txt"),'\t')
    matlab_w        = freqs_eg1_w_mag_phasedeg[:,1]
    matlab_mag      = freqs_eg1_w_mag_phasedeg[:,2]
    matlab_phasedeg = freqs_eg1_w_mag_phasedeg[:,3]

    # Test
    @test w ≈ matlab_w
    @test mag ≈ matlab_mag
    @test phasedeg ≈ matlab_phasedeg

    @test h ≈ freqs(ZeroPoleGain(PolynomialRatio(b, a)), w)
    @test h ≈ freqs(SecondOrderSections(PolynomialRatio(b, a)), w)

    #=using Winston=#
    #=figure = loglog(w, mag)=#
    #=ylabel("Magnitude")=#
    #=xlabel("Frequency (rad/s)")=#
    #=file(figure, "MATLAB-freqs-mag.png", width=1200, height=800)=#

    #=figure = semilogx(w, phasedeg)=#
    #=ylim(-150, 0)=#
    #=setattr(figure, draw_grid=true, tickdir=1)=#
    #=ylabel("Phase (degrees)")=#
    #=xlabel("Frequency (rad/s)")=#
    #=file(figure, "MATLAB-freqs-phase.png", width=1200, height=800)=#
end

#######################################
#
#  Test analog filter with frequency specified in hz
#
#  TODO: Create a MATLAB ground truth to compare against
#        Currently this just checks it runs without error
#
#######################################
@testset "freqs with fs" begin
    # Julia
    a  = [1.0, 0.4, 1.0]
    b  = [0.2, 0.3, 1.0]
    fs = 8192
    hz = Compat.range(0, stop=fs, length=50)

    h        = freqs(PolynomialRatio(b, a), hz, fs)
    mag      = convert(Array{Float64}, abs.(h))
    phasedeg = (180/pi)*convert(Array{Float64}, angle.(h))

    #=using Winston=#
    #=figure = semilogx(hz, mag)=#
    #=ylabel("Magnitude")=#
    #=xlabel("Frequency (Hz)")=#
    #=file(figure, "MATLAB-freqs-hz.png", width=1200, height=800)=#
end
