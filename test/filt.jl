!(dirname(@__FILE__) in LOAD_PATH) && push!(LOAD_PATH, dirname(@__FILE__))
using DSP, Test, Random, FilterTestHelpers
#
# filt with different filter forms
#
@testset "filt" begin
    x = randn(100)

    for n = 1:6
        for proto in (Butterworth(n), Chebyshev1(n, 1), Chebyshev2(n, 1), Elliptic(n, 0.5, 2))
            zpk = digitalfilter(Lowpass(0.2), proto)
            tf = convert(PolynomialRatio, zpk)
            if n <= 2
                bq = convert(Biquad, zpk)
            else
                @test_throws ArgumentError convert(Biquad, zpk)
            end
            sos = convert(SecondOrderSections, zpk)

            res = filt(sos, x)

            # Test with filt with tf/sos
            tfres = filt(tf, x)
            @test res ≈ tfres
            @test res ≈ filt!(similar(x), sos, x)
            @test res ≈ filt!(similar(x), tf, x)

            # For <= 2 poles, test with biquads
            if n <= 2
                @test res ≈ filt(bq, x)
                @test res ≈ filt!(similar(x), bq, x)
                f = DF2TFilter(bq)
                @test tfres ≈ [filt(f, x[1:50]); filt(f, x[51:end])]
            end

            # Test that filt with zpk converts
            @test res ≈ filt(zpk, x)
            @test res ≈ filt!(similar(x), zpk, x)

            # Test with DF2TFilter
            f = DF2TFilter(sos)
            @test res ≈ [filt(f, x[1:50]); filt(f, x[51:end])]
            f = DF2TFilter(tf)
            @test tfres ≈ [filt(f, x[1:50]); filt(f, x[51:end])]
            f = DF2TFilter(zpk)
            @test res ≈ [filt(f, x[1:50]); filt(f, x[51:end])]
        end
    end

    # Test simple scaling with DF2TFilter
    @test filt(DF2TFilter(PolynomialRatio([3.7], [4.2])), x) == x * (3.7/4.2)

    # DF2TFilter{<:PolynomialRatio} with unequal numerator/denominator orders (issue #436)
    @test filt(DF2TFilter(PolynomialRatio([0, 0, 1, 0.8], [1])), [1; zeros(9)]) == [0; 0; 1; 0.8; zeros(6)]
    @test filt([0, 0, 1, 0.8], [1], [1; zeros(9)]) == [0; 0; 1; 0.8; zeros(6)]
    @test filt(DF2TFilter(PolynomialRatio([1], [1, -0.5])), [1; zeros(9)]) ≈ 0.5.^(0:9)
    @test filt([1], [1, -0.5], [1; zeros(9)]) ≈ 0.5.^(0:9)
    @test_nowarn filt(DF2TFilter(PolynomialRatio([1:5;], [1, -0.5])), ones(10))
    @test_nowarn filt(DF2TFilter(PolynomialRatio(rand(DSP.SMALL_FILT_CUTOFF + 1), [1])), ones(10))

    # DF2TFilter{:z} state type selection for ZeroPoleGain (issue #371)
    s = rand(30) + im * rand(30)
    df = digitalfilter(Lowpass(0.25), Butterworth(4))
    f = @test_nowarn DF2TFilter(df, ComplexF64)
    @test_nowarn filt(f, s)

    # DF2TFilter{ZPG/SOS} stackoverflow error
    @test_throws MethodError DF2TFilter(ZeroPoleGain([1], [2], 3), :D, 'x', Ref(1))
end

@testset "multi-column filt $D-D" for D in 1:4
     b = [0.1, 0.1]
     a = [1.0, -0.8]
     sz = (10, ntuple(n -> n+1, Val(D))...)
     y_ref = filt(b, a, ones(sz[1]))
     x = ones(sz)
     slicedims = ntuple(n -> n+1, Val(D))
     @test all(col -> col ≈ y_ref, eachslice(filt(b, a, x); dims=slicedims))
     @test all(col -> col ≈ y_ref, eachslice(filt(PolynomialRatio(b, a), x); dims=slicedims))
     @test all(col -> col ≈ y_ref, eachslice(filt(Biquad(PolynomialRatio(b, a)), x); dims=slicedims))
     @test all(col -> col ≈ y_ref, eachslice(filt(SecondOrderSections(PolynomialRatio(b, a)), x); dims=slicedims))
     # with si given
     @test all(col -> col ≈ y_ref, eachslice(@test_deprecated(filt(b, a, x, zeros(1, sz[2:end]...))); dims=slicedims))
     @test all(col -> col ≈ y_ref, eachslice(@test_deprecated(filt(PolynomialRatio(b, a), x, zeros(1, sz[2:end]...))); dims=slicedims))
     @test all(col -> col ≈ y_ref, eachslice(@test_deprecated(filt(Biquad(PolynomialRatio(b, a)), x, zeros(2, sz[2:end]...))); dims=slicedims))
     @test all(col -> col ≈ y_ref, eachslice(@test_deprecated(filt(SecondOrderSections(PolynomialRatio(b, a)), x, zeros(2, 1, sz[2:end]...))); dims=slicedims))
     # use _small_filt_fir!
     b = [0.1, 0.1]
     a = [1.0]
     y_ref = filt(b, a, ones(sz[1]))
     @test all(col -> col ≈ y_ref, eachslice(filt(b, a, x); dims=slicedims))
     @test all(col -> col ≈ y_ref, eachslice(filt(PolynomialRatio(b, a), x); dims=slicedims))
end

@testset "multi-column DF2TFilter $D-D" for D in 1:4
    b = [0.1, 0.1]
    a = [1.0, -0.8]
    sz = (10, ntuple(n -> n+1, Val(D))...)
    y_ref = filt(b, a, ones(2*sz[1]))
    x = ones(sz)

    H = DF2TFilter(PolynomialRatio(b, a),sz[2:end])
    @test all(col -> col ≈ y_ref[1:sz[1]], eachslice(filt(H, x); dims=ntuple(n -> n+1, Val(D))))
    @test all(col -> col ≈ y_ref[sz[1]+1:end], eachslice(filt(H, x); dims=ntuple(n -> n+1, Val(D))))

    H = DF2TFilter(SecondOrderSections(PolynomialRatio(b, a)), sz[2:end])
    @test all(col -> col ≈ y_ref[1:sz[1]], eachslice(filt(H, x); dims=ntuple(n -> n+1, Val(D))))
    @test all(col -> col ≈ y_ref[sz[1]+1:end], eachslice(filt(H, x); dims=ntuple(n -> n+1, Val(D))))

    H = DF2TFilter(Biquad(PolynomialRatio(b, a)), sz[2:end])
    @test all(col -> col ≈ y_ref[1:sz[1]], eachslice(filt(H, x); dims=ntuple(n -> n+1, Val(D))))
    @test all(col -> col ≈ y_ref[sz[1]+1:end], eachslice(filt(H, x); dims=ntuple(n -> n+1, Val(D))))
end

#
# filtfilt
#

@testset "filt_stepstate" begin
    ##############
    #
    # Filter initial conditions
    # Python example 1 - http://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lfilter_zi.html
    #
    # b, a = butter(5, 0.25)
    # zi = lfilter_zi(b, a)
    #
    ##############
    zi_python = [ 0.99672078, -1.49409147,  1.28412268, -0.45244173,  0.07559489]

    b = [ 0.00327922,  0.01639608,  0.03279216,  0.03279216,  0.01639608,  0.00327922]
    a = [ 1.        , -2.47441617,  2.81100631, -1.70377224,  0.54443269, -0.07231567]

    @test ≈(zi_python, DSP.Filters.filt_stepstate(b, a)[1], atol=1e-7)

    ##############
    #
    # Filter initial conditions
    # Matlab - Random filter coefficients
    #
    # filtfilt([0.222, 0.43, 0.712], [1, 0.33, 0.22], x)
    # Then used breakpoints to extract zi
    #
    ##############

    zi_matlab = [0.6580, 0.5184]

    b = [0.222, 0.43, 0.712]
    a = [1, 0.33, 0.22]

    @test zi_matlab ≈ DSP.Filters.filt_stepstate(b, a)[1]


    ##############
    #
    # Filter initial conditions
    # Python - non 1 first coeffecient
    #
    # b = array([ 0.00327922,  0.01639608,  0.03279216,  0.03279216,  0.01639608,  0.00327922])
    # a = array([ 1.1       , -2.47441617,  2.81100631, -1.70377224,  0.54443269, -0.07231567])
    # zi = lfilter_zi(b, a)
    #
    ##############

    zi_python = [0.55996501, -0.72343165,  0.68312446, -0.2220676 ,  0.04030775]

    b = [ 0.00327922,  0.01639608,  0.03279216,  0.03279216,  0.01639608,  0.00327922]
    a = [ 1.1       , -2.47441617,  2.81100631, -1.70377224,  0.54443269, -0.07231567]

    @test ≈(zi_python, DSP.Filters.filt_stepstate(b, a)[1], atol=1e-7)
end


##############
#
# Filter in place check (with initial conditions)
#
# x = '/Users/rluke/.julia/v0.3/DSP/test/data/spectrogram_x.txt';
# filtcheck  = filter([0.4, 1], [0.9, 0.6], x, 0.4750)
# dlmwrite('filt_check.txt',[filtcheck], 'delimiter', '\t', 'precision', '%.12f')
#
##############
@testset "filt! with init. cond." begin
    matlab_filt = reference_data("filt_check.txt")

    a = [0.9, 0.6]
    b = [0.4, 1]
    z = [0.4750]
    x = vec(reference_data("spectrogram_x.txt"))
    @test_deprecated(filt!(x, b, a, x, z))

    @test matlab_filt ≈ x

    x = vec(reference_data("spectrogram_x.txt"))
    filt!(x, DF2TFilter(PolynomialRatio(b, a), z), x)

    @test matlab_filt ≈ x

    # With initial conditions: a lowpass 5-pole butterworth filter with W_n = 0.25,
    # and a stable initial filter condition matched to the initial value
    zpg = digitalfilter(Lowpass(0.25), Butterworth(5))
    si = [0.9967207836936347, -1.4940914728163142, 1.2841226760316475, -0.4524417279474106, 0.07559488540931815]
    @test filt(DF2TFilter(PolynomialRatio(zpg), si), ones(10)) ≈ ones(10) # Shouldn't affect DC offset
end

#######################################
#
# Test 1d filtfilt against matlab results
#
# x = '/Users/rluke/.julia/v0.3/DSP/test/data/spectrogram_x.txt'; x = textread(x);
# b = [ 0.00327922,  0.01639608,  0.03279216,  0.03279216,  0.01639608,  0.00327922];
# a = [ 1.        , -2.47441617,  2.81100631, -1.70377224,  0.54443269, -0.07231567];
# x2 = filtfilt(b, a, x);
# dlmwrite('filtfilt_output.txt',[x2], 'delimiter', '\t', 'precision', '%.12f')
#
#######################################
@testset "filtfilt 1D" begin
    x2_matlab = reference_data("filtfilt_output.txt")

    b = [ 0.00327922,  0.01639608,  0.03279216,  0.03279216,  0.01639608,  0.00327922]
    a = [ 1.        , -2.47441617,  2.81100631, -1.70377224,  0.54443269, -0.07231567]
    x = reference_data("spectrogram_x.txt")

    @test x2_matlab ≈ filtfilt(b, a, x)
end

# Make sure above doesn't crash for real coeffs & complex data.
@testset "filtfilt 1D Complex" begin
    x2_matlab = reference_data("filtfilt_output.txt")

    b = [ 0.00327922,  0.01639608,  0.03279216,  0.03279216,  0.01639608,  0.00327922]
    a = [ 1.        , -2.47441617,  2.81100631, -1.70377224,  0.54443269, -0.07231567]
    x = reference_data("spectrogram_x.txt")

    y = x .+ 1im .* randn(size(x, 1))

    @test x2_matlab ≈ real.(filtfilt(b, a, y))
end

#######################################
#
# Test 2d filtfilt against matlab results
#

#=x = '/Users/rluke/.julia/v0.3/DSP/test/data/spectrogram_x.txt'; x = textread(x);=#

#=x = repmat(x, 1, 3);=#
#=x(:,2) = circshift(x(:,2), 64);=#
#=x(:,3) = circshift(x(:,3), 128);=#

#=b = [ 0.00327922,  0.01639608,  0.03279216,  0.03279216,  0.01639608,  0.00327922];=#
#=a = [ 1.        , -2.47441617,  2.81100631, -1.70377224,  0.54443269, -0.07231567];=#
#=x2 = filtfilt(b, a, x);=#
#=dlmwrite('filtfilt_output_2d.txt',[x2], 'delimiter', '\t', 'precision', '%.12f')=#

#
#######################################

@testset "filtfilt 2D" begin
    b = [ 0.00327922,  0.01639608,  0.03279216,  0.03279216,  0.01639608,  0.00327922]
    a = [ 1.        , -2.47441617,  2.81100631, -1.70377224,  0.54443269, -0.07231567]

    x2_output = reference_data("filtfilt_output_2d.txt")

    x = reference_data("spectrogram_x.txt")
    x = repeat(x, outer=(1, 3))
    x[:,2] = circshift(x[:,2], 64)
    x[:,3] = circshift(x[:,3], 128)

    @test x2_output ≈ filtfilt(b, a, x)

    #######################################
    #
    # Test 2d filtfilt with filter type
    #
    #
    #######################################

    f = PolynomialRatio(b, a)

    # Use 2d data from last test
    @test x2_output ≈ filtfilt(f, x)
end

#######################################
#
# Test 2d filtfilt with SecondOrderSections
#
# Our implementation differs from MATLAB, but should match the
# PolynomialRatio provided the filter order is even. (Otherwise
# the extrapolation will differ slightly.)
#######################################
@testset "filtfilt SOS" begin
    x = reference_data("spectrogram_x.txt")

    f = digitalfilter(Lowpass(0.2), Butterworth(4))
    @test filtfilt(convert(SecondOrderSections, f), x) ≈ filtfilt(convert(PolynomialRatio, f), x)

    f = digitalfilter(Highpass(0.1), Butterworth(6))
    @test filtfilt(convert(SecondOrderSections, f), x) ≈ filtfilt(convert(PolynomialRatio, f), x)

    f = digitalfilter(Bandpass(0.1, 0.3), Butterworth(2))
    @test filtfilt(convert(SecondOrderSections, f), x) ≈ filtfilt(convert(PolynomialRatio, f), x)
    @test_nowarn filtfilt(f, [1.])  # check that pad_length won't throw oob
end

#
# fftfilt/filt
#

# make sure to include blen > SMALL_FILT_CUTOFF
@testset "fftfilt $xlen/$blen" for xlen in 2 .^ (7:18) .- 1, blen in 2 .^ (1:7) .- 1
    b = randn(blen)
    for x in (rand(xlen), rand(xlen, 2))
        out = similar(x)
        filtres = filt(b, [1.0], x)
        fftres = fftfilt(b, x)
        firres = filt(b, x)
        td_res = tdfilt(b, x)

        @test filtres ≈ fftres
        @test filtres ≈ firres
        @test filtres ≈ td_res

        fftfilt!(out, b, x)
        @test filtres ≈ out # test output of fftfilt!

        tdfilt!(out, b, x)
        @test filtres ≈ out # test output of tdfilt!
    end
end

# fir_filtfilt
@testset "fir_filtfilt" begin
    for b in (randn(10), [1:10;]), x in (randn(100), randn(100, 2), randn(10))   # also tests issue #537
        @test filtfilt(b, x) ≈ DSP.Filters.iir_filtfilt(b, [1], x)
        @test filtfilt(b, [2.0], x) ≈ DSP.Filters.iir_filtfilt(b, [2], x)
    end
end
