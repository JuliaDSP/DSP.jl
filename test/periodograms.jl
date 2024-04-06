#in matlab:
#x=rand(512,1);
#[s,f,t,p]=spectrogram(x,ones(1,256),128,256,10);
#save
#
#in julia:
#using MAT
#
#matdata=matread("matlab.mat")
#
#for i in ("x", "f", "t", "p")
#  fid=open("spectrogram_$i.txt","w")
#  print(fid,matdata["$i"])
#  close(fid)
#end

using DSP, Test
using Statistics: mean
using FFTW: fftfreq

@testset "matlab ref" begin
    x0 = vec(readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_x.txt"),'\t'))
    f0 = vec(readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_f.txt"),'\t'))
    t0 = vec(readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_t.txt"),'\t'))
    p0 = readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_p.txt"),'\t')
    spec = spectrogram(x0, 256, 128; fs=10)
    p, f, t = power(spec), freq(spec), time(spec)

    # with real input matlab outputs a 1-sided PSD
    @test p0 ≈ p
    @test f0 ≈ f
    @test t0 ≈ t


    mt_spec = mt_spectrogram(x0, 256, 128; fs=10)
    @test freq(mt_spec) == f
    @test time(mt_spec) == t
    @test power(mt_spec)[:,1] ≈ power(mt_pgram(x0[1:256]; fs=10))

    # in-place full precision version:
    spec_config = MTSpectrogramConfig{Float64}(length(x0), 256, 128; fs=10)
    out = allocate_output(spec_config)
    mt_spec2 = mt_spectrogram!(out, x0, spec_config)
    @test power(mt_spec2) ≈ power(mt_spec)
    @test freq(mt_spec2) == freq(mt_spec)
    @test time(mt_spec2) == time(mt_spec)

    # out-of-place with config:
    r = mt_spectrogram(x0, spec_config)
    @test power(r) ≈ power(mt_spec)
    @test freq(r) == freq(mt_spec)
    @test time(r) == time(mt_spec)

    # in-place without config:
    r = mt_spectrogram!(out, x0, 256, 128; fs=10)
    @test power(r) ≈ power(mt_spec)
    @test freq(r) == freq(mt_spec)
    @test time(r) == time(mt_spec)

    # in-place half precision version:
    spec_config = MTSpectrogramConfig{Float32}(length(x0), 256, 128; fs=10)
    out = allocate_output(spec_config)
    mt_spec3 = mt_spectrogram!(out, x0, spec_config)
    @test power(mt_spec3) ≈ power(mt_spec)
    @test freq(mt_spec3) == freq(mt_spec)
    @test time(mt_spec3) == time(mt_spec)

    mt_spec4 = mt_spectrogram!(out, Float32.(x0), spec_config)
    @test power(mt_spec4) ≈ power(mt_spec)
    @test freq(mt_spec4) == freq(mt_spec)
    @test time(mt_spec4) == time(mt_spec)

    # We can also only pass the window config. Full precision:
    config = MTConfig{Float64}(256; fs=10)
    mt_spec5 = mt_spectrogram(x0, config, 128)
    @test power(mt_spec5) ≈ power(mt_spec)
    @test freq(mt_spec5) == freq(mt_spec)
    @test time(mt_spec5) == time(mt_spec)

    # We can also only pass the window config. Half precision:
    config = MTConfig{Float32}(256; fs=10)
    mt_spec6 = mt_spectrogram(x0, config, 128)
    @test power(mt_spec6) ≈ power(mt_spec)
    @test freq(mt_spec6) == freq(mt_spec)
    @test time(mt_spec6) == time(mt_spec)

    # with Float32 input:
    mt_spec7 = mt_spectrogram(Float32.(x0), config, 128)
    @test power(mt_spec7) ≈ power(mt_spec)
    @test freq(mt_spec7) == freq(mt_spec)
    @test time(mt_spec7) == time(mt_spec)

    @test_throws DimensionMismatch mt_spectrogram!(similar(out, size(out, 1), size(out, 2)+1), x0, spec_config)
    @test_throws DimensionMismatch mt_spectrogram!(out, vcat(x0, x0), spec_config)
end

@testset "0:7" begin
    #Matlab: p = pwelch(0:7, [1, 1, 1, 1, 1, 1, 1, 1], 0, 8, 1, 'twosided')
    data = 0:7
    data0 = Float64[98.0,
                    13.656854249492380,
                    4.0,
                    2.343145750507620,
                    2.0,
                    2.343145750507620,
                    4.0,
                    13.656854249492380]
    @test power(periodogram(data, onesided=false)) ≈ data0
    @test power(welch_pgram(data, length(data), 0, onesided=false)) ≈ data0
    @test power(spectrogram(data, length(data), 0, onesided=false)) ≈ data0
    @test power(periodogram(complex.([data;], [data;]), onesided=false)) ≈ data0*2
    @test power(welch_pgram(complex.([data;], [data;]), length(data), 0, onesided=false)) ≈ data0*2
    @test power(spectrogram(complex.([data;], [data;]), length(data), 0, onesided=false)) ≈ data0*2

    # # ~~~~~~~~ Tests with no window ~~~~~~~~~~~~~~~~~~~
    # Matlab: p = pwelch(0:7, [1, 1], 0, 2, 1, 'twosided')
    expected = Float64[34.5, 0.5]
    @test power(welch_pgram(data, 2, 0; onesided=false)) ≈ expected
    @test mean(power(spectrogram(data, 2, 0; onesided=false)), dims=2) ≈ expected

    # Matlab: p = pwelch(0:7, [1, 1, 1], 0, 3, 1, 'twosided')
    expected = Float64[25.5, 1.0, 1.0]
    @test power(welch_pgram(data, 3, 0; onesided=false)) ≈ expected
    @test mean(power(spectrogram(data, 3, 0; onesided=false)), dims=2) ≈ expected

    # Matlab: p = pwelch(0:7, [1, 1, 1], 1, 3, 1, 'twosided')
    expected = Float64[35.0, 1.0, 1.0]
    @test power(welch_pgram(data, 3, 1; onesided=false)) ≈ expected
    @test mean(power(spectrogram(data, 3, 1; onesided=false)), dims=2) ≈ expected

    # Matlab: p = pwelch(0:7, [1, 1, 1, 1], 1, 4, 1, 'twosided')
    expected = Float64[45, 2, 1, 2]
    @test power(welch_pgram(data, 4, 1; onesided=false)) ≈ expected
    @test mean(power(spectrogram(data, 4, 1; onesided=false)), dims=2) ≈ expected

    # ~~~~~~~~~~~ This one tests periodogram ~~~~~~~~~~~~
    # ~ If functionality of the other arguments has been
    # ~ tested above, we only test here that the correct
    # ~ value of the spectral density is obtained when
    # ~ using a window. More tests to be added if needed
    #Matlab: p = pwelch(0:7, window_func(8), 0, 8, 1, 'twosided')
    cases = (
        (hamming,  Float64[65.461623986801527,
                           20.556791795515764,
                            0.369313143650544,
                            0.022167446610882,
                            0.025502985564107,
                            0.022167446610882,
                            0.369313143650544,
                           20.556791795515764]),
        (bartlett, Float64[62.999999999999993,
                           21.981076052592442,
                            0.285714285714286,
                            0.161781090264695,
                            0.142857142857143,
                            0.161781090264695,
                            0.285714285714286,
                           21.981076052592442])
    )

    @testset for (window1, expected) in cases
        @test power(periodogram(data; window=window1, onesided=false)) ≈ expected
        @test power(welch_pgram(data, length(data), 0; window=window1, onesided=false)) ≈ expected
        @test power(spectrogram(data, length(data), 0; window=window1, onesided=false)) ≈ expected
        @test power(periodogram(data; window=window1(length(data)), onesided=false)) ≈ expected
        @test power(welch_pgram(data, length(data), 0; window=window1(length(data)), onesided=false)) ≈ expected
        @test power(spectrogram(data, length(data), 0; window=window1(length(data)), onesided=false)) ≈ expected
    end

    # Padded periodogram
    # MATLAB: a = periodogram(0:7, [], 32);
    expected = [
                      98
        174.463067389405
        121.968086934209
        65.4971744936088
        27.3137084989848
        12.1737815028909
        10.3755170959439
        10.4034038628775
                       8
        5.25810953219633
        4.47015397150535
        4.89522578856669
        4.68629150101524
        3.69370284475603
         3.1862419983415
        3.61553458569862
                       2
    ]
    @test power(periodogram(data; nfft=32)) ≈ expected
    @test power(welch_pgram(data, length(data), 0; nfft=32)) ≈ expected
    @test power(spectrogram(data, length(data), 0; nfft=32)) ≈ expected

    # Padded periodogram with window
    # MATLAB: a = periodogram(0:7, hamming(8), 32, 1)
    expected = [
          65.4616239868015
          122.101693164395
          98.8444689598445
           69.020252632913
          41.1135835910315
          20.5496474310966
          8.43291449161938
          2.78001620362588
         0.738626287301088
         0.174995741770789
        0.0501563022944516
        0.0327357460012861
        0.0443348932217643
        0.0553999745503552
        0.0561319901616643
        0.0526025934871384
        0.0255029855641069
    ]
    @test power(periodogram(data; window=hamming, nfft=32)) ≈ expected
    @test power(welch_pgram(data, length(data), 0; window=hamming, nfft=32)) ≈ expected
    @test power(spectrogram(data, length(data), 0; window=hamming, nfft=32)) ≈ expected

    # test welch_pgram configuration object
    expected = power(welch_pgram(data, length(data), 0; window=hamming, nfft=32))
    config = WelchConfig(data; n=length(data), noverlap=0, window=hamming, nfft=32)
    @test power(welch_pgram(data, config)) == expected

    # test welch_pgram!
    out = similar(expected)
    @test power(welch_pgram!(out, data, config)) == expected
    @test power(welch_pgram!(out, data, length(data), 0; window=hamming, nfft=32)) == expected
    @test_throws ArgumentError welch_pgram!(convert(Vector{Float32}, out), data, config)
    @test_throws DimensionMismatch welch_pgram!(empty!(out), data, config)
    
    # test welch_pgram! with float64 data
    out = similar(expected)
    @test power(welch_pgram!(out, convert(UnitRange{Float64}, data), config)) == expected 

    # Test fftshift
    p = periodogram(data); p_shifted = fftshift(p)
    @test power(p) == power(p_shifted)
    @test freq(p) ≈ freq(p_shifted)
    @test p_shifted == fftshift(p_shifted)

    p = periodogram(data; onesided=false)
    @test fftshift(power(p)) == power(fftshift(p))
    @test fftshift(freq(p)) == freq(fftshift(p))
end

@testset "fftshift" begin
    data = 1:100

    p = spectrogram(data); p_shifted = fftshift(p)
    @test power(p) == power(p_shifted)
    @test freq(p) ≈ freq(p_shifted)
    @test p_shifted == fftshift(p_shifted)

    p = spectrogram(data; onesided=false)
    @test fftshift(power(p), 1) == power(fftshift(p))
    @test fftshift(freq(p)) == freq(fftshift(p))

    using DSP.Periodograms: Periodogram2

    # for coverage, not very useful...
    p = Periodogram2([1 2; 3 4], 1:2, fftfreq(2)); ps = fftshift(p)
    @test ps.power == fftshift(p.power, 2)
    @test (ps.freq1, ps.freq2) == (p.freq1, fftshift(p.freq2))
end

@testset "2D" begin
    data2d = readdlm(joinpath(dirname(@__FILE__), "data", "per2dx.txt"),'\t')
    expectedsum = vec(readdlm(joinpath(dirname(@__FILE__), "data", "per2dsum.txt"),'\t'))
    expectedmean = vec(readdlm(joinpath(dirname(@__FILE__), "data", "per2dmean.txt"),'\t'))
    # 2-d periodgram (radialsum)
    # computed in octave with raPsd2d ((C) E. Ruzanski) replacing nanmean with nansum
    # P = raPsd2d(x,1)'*n^2
    @test power(periodogram(data2d,fs=1, radialsum=true)) ≈ expectedsum

    # 2-d periodgram (radialavg)
    # computed in octave with raPsd2d ((C) E. Ruzanski)
    # P = raPsd2d(x,1)'*n^2
    @test power(periodogram(data2d, fs=1, radialavg=true)) ≈ expectedmean

    # 2-d periodgram 2-d PSD
    @test power(periodogram(data2d, fs=1)) ≈ abs2.(fft(data2d))*1/prod(size(data2d))
    # 2-d periodgram 2-d PSD with padding
    pads = (size(data2d,1)+4,size(data2d,1)+7)
    data2dpad = zeros(Float64,pads...)
    data2dpad[1:size(data2d,1),1:size(data2d,2)] = data2d
    @test power(periodogram(data2d, fs=1, nfft=pads)) ≈ abs2.(fft(data2dpad))*1/prod(size(data2d))
    # 2-d periodgram radial freq
    @test freq(periodogram(data2d, fs=3.3, radialsum=true)) ≈ freq(periodogram(vec(data2d[1,:]), fs=3.3))
    # 2-d periodgram 2-d freq
    f1,f2 = freq(periodogram(data2d, fs=3.3))
    f1d = freq(periodogram(vec(data2d[1,:]), fs=3.3, onesided=false))
    @assert size(data2d,1)==size(data2d,2)
    for j=1:size(data2d,2)
        for i=1:size(data2d,1)
            @test  [f1[i],f2[j]] ≈ [f1d[i],f1d[j]]
        end
    end
    # Test fftshift
    p = periodogram(data2d); p_shifted = fftshift(p)
    @test fftshift(power(p)) == power(fftshift(p))
    @test p_shifted == fftshift(p_shifted)
    f = freq(p)
    @test (fftshift(f[1]),fftshift(f[2])) == freq(p_shifted)
end

@testset "radial" begin
    # 2-d periodgram radial test for a non-square signal sparse in fft space
    n1 = 52
    n2 = 46  # assuming n1>n2
    nf = (22,7) # the non-zero location
    F = (fftfreq(n1,1),fftfreq(n2,1))
    a = [F[1][nf[1]],F[2][nf[2]]]
    FB = Bool[[F[1][i], F[2][j]]==a || [F[1][i], F[2][j]]==-a for i = 1:n1, j = 1:n2]
    ind = findall(FB)
    x = zeros(n1,n2)*0im;
    x[ind] = [1+2im,1-2im]
    y = real(ifft(x))

    fwn = round(Int, sqrt((a[1])^2+(a[2])^2)*n2)
    pe = zeros(n2>>1 + 1)
    pe[fwn+1] = 2*abs2(x[nf...])/n1/n2
    P = periodogram(y,nfft=(n1,n2),radialsum=true)
    @test power(P) ≈ pe
    @test freq(P)[fwn+1] ≈ fwn/n2
end

# Testing STFT function and comparing results with MATLAB
@testset "stft" begin
    fs = 16000
    nfft = 512
    nwin = 400
    nhop = 160
    s = vec(readdlm(joinpath(dirname(@__FILE__), "data", "stft_x.txt"),'\t'))

    Sjl = stft(s, nwin, nwin-nhop; nfft=nfft, fs=fs, window=hanning)
    Sml_re = readdlm(joinpath(dirname(@__FILE__), "data", "stft_S_real.txt"),'\t')
    Sml_im = readdlm(joinpath(dirname(@__FILE__), "data", "stft_S_imag.txt"),'\t')
    Sml = complex.(Sml_re, Sml_im)
    @test Sjl ≈ Sml
end

@testset "fft2oneortwosided!" begin
    n = 10
    floatt = Float64
    for onesided in (true, false),
            nfft in (n, n+2, n+3),
           atype in (floatt, Complex{floatt})
        nout = nout = onesided ? (nfft >> 1)+1 : nfft
        x = zeros(atype, nfft)
        if atype <: Real
            x[1:n] = rand(atype, n)
            xrcfft = rfft(x)
        else
            x[1:n] = rand(floatt, n)+im*rand(floatt, n)
            xrcfft = fft(x)
        end
        xfft = fft(x)
        out = zeros(fftouttype(atype),nout,3)
        if !(onesided == true && atype <: Complex)
            outft = DSP.Periodograms.fft2oneortwosided!(out, xrcfft, nfft, onesided, nout)
        end
        if onesided == true && atype <: Real
            @test out[:,2] ≈ xrcfft
            @test out[:,[1,3]] ≈ xrcfft.*[0 0]
        elseif onesided == false && atype <: Real
            @test out[:,2] ≈ xfft
            @test out[:,[1,3]] ≈ xfft.*[0 0]
        elseif onesided == false && atype <: Complex
            @test out[:,2] ≈ xfft
            @test out[:,[1,3]] ≈ xfft.*[0 0]
        else
            #onesided complex
        end
    end
end

@testset "mt_pgram" begin
    # MATLAB: x = pmtm(stft_x, 4, 5000, 16000, 'unity')
    s = vec(readdlm(joinpath(dirname(@__FILE__), "data", "stft_x.txt"),'\t'))
    mtdata = vec(readdlm(joinpath(dirname(@__FILE__), "data", "mt_pgram.txt")))
    @test power(mt_pgram(s; fs=16000)) ≈ mtdata
    @test power(mt_pgram(s; fs=16000, window=dpss(length(s), 4))) ≈ mtdata

    # error tests
    @test_throws ArgumentError periodogram([1 2 3])
    @test_throws ArgumentError periodogram(rand(2,3), nfft=(3,2))
    @test_throws ArgumentError periodogram([1 2;3 4],radialsum=true, radialavg=true)

    # #124
    q = arraysplit(ones(Float64, 1000),100,10);
    @test map(mean, q) == ones(Float64, 11)

    # test that iterating ArraySplit always yields the same Vector
    q = arraysplit(-10:10,4,2)
    for x in q
        @assert isa(x, Vector)
        @test pointer(x) == pointer(q.buf)
    end

    x = vec(readdlm(joinpath(dirname(@__FILE__), "data", "pmtm_x.txt")))
    @test eltype(x) <: Float64

    # MATLAB code:
    # fft = 2^nextpow2(length(x));
    # [pxx,fx] = pmtm(x,4,nfft,1000,'unity');
    fx = vec(readdlm(joinpath(dirname(@__FILE__), "data", "pmtm_fx.txt")))
    pxx = vec(readdlm(joinpath(dirname(@__FILE__), "data", "pmtm_pxx.txt")))

    fs = 1000
    nw=4
    nfft = nextpow(2, length(x))
    result = mt_pgram(x; fs=fs, nw=nw, nfft=nfft)
    @test freq(result) ≈ fx
    @test power(result) ≈ pxx

    # Test against in-place. Full precision:
    config = MTConfig{Float64}(length(x); fs=fs, nw=nw, nfft=nfft)
    out = allocate_output(config)
    @test eltype(out) == Float64
    result2 = mt_pgram!(out, x, config)
    @test freq(result2) ≈ fx
    @test power(result2) ≈ pxx

    # in-place without config
    r = mt_pgram!(out, x; fs=fs, nw=nw, nfft=nfft)
    @test freq(r) ≈ fx
    @test power(r) ≈ pxx

    # out-of-place with config
    r = mt_pgram(x, config)
    @test freq(r) ≈ fx
    @test power(r) ≈ pxx

    # Lower precision output:
    config = MTConfig{Float32}(length(x); fs=fs, nw=nw, nfft=nfft)
    out = allocate_output(config)
    @test eltype(out) == Float32
    result3 = mt_pgram!(out, x, config)
    @test freq(result3) ≈ fx
    @test power(result3) ≈ pxx

    # with Float32 input:
    result4 = mt_pgram!(out, Float32.(x), config)
    @test freq(result4) ≈ fx
    @test power(result4) ≈ pxx

    @test_throws DimensionMismatch mt_pgram!(similar(out, length(out)+1), x, config)
    @test_throws DimensionMismatch mt_pgram!(out, vcat(x, one(eltype(x))), config)

    y = vec(readdlm(joinpath(dirname(@__FILE__), "data", "pmtm_y.txt")))
    z = x + im*y
    @test eltype(z) <: Complex{Float64}

    # MATLAB code: `[pzz,fz] = pmtm(z,4,nfft,1000, 'unity')`
    fz = vec(readdlm(joinpath(dirname(@__FILE__), "data", "pmtm_fz.txt")))
    pzz = vec(readdlm(joinpath(dirname(@__FILE__), "data", "pmtm_pzz.txt")))
    result = mt_pgram(z; fs=fs, nw=nw, nfft=nfft)
    result_mask = 0 .< freq(result) .< 500
    freqs = freq(result)[result_mask]
    @test freqs ≈ fz[2:length(freqs)+1]
    @test power(result)[result_mask] ≈ pzz[2:length(freqs)+1]

    # Test against in-place. Full precision:
    config = MTConfig{Complex{Float64}}(length(z); fs=fs, nw=nw, nfft=nfft)
    out = allocate_output(config)
    @test eltype(out) == Float64

    result2 = mt_pgram!(out, z, config)
    @test freq(result2) == freq(result)
    @test power(result2) ≈ power(result)

    # Lower precision output:
    config = MTConfig{Complex{Float32}}(length(z); fs=fs, nw=nw, nfft=nfft)
    out = allocate_output(config)
    @test eltype(out) == Float32

    result2 = mt_pgram!(out, z, config)
    @test freq(result2) == freq(result)
    @test power(result2) ≈ power(result)

    # Lower precision computation:
    result3 = mt_pgram!(out, Complex{Float32}.(z), config)
    @test freq(result3) == freq(result)
    @test power(result3) ≈ power(result)

end
