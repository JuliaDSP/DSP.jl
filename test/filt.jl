!(dirname(@__FILE__) in LOAD_PATH) && push!(LOAD_PATH, dirname(@__FILE__))
using DSP, Base.Test, FilterTestHelpers

#
# filt with different filter forms
#
srand(1776)
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
        @test_approx_eq res tfres
        @test_approx_eq res filt!(similar(x), sos, x)
        @test_approx_eq res filt!(similar(x), tf, x)

        # For <= 2 poles, test with biquads
        if n <= 2
            @test_approx_eq res filt(bq, x)
            @test_approx_eq res filt!(similar(x), bq, x)
            f = DF2TFilter(bq)
            @test tfres == [filt(f, x[1:50]); filt(f, x[51:end])]
        end

        # Test that filt with zpk converts
        @test res == filt(zpk, x)
        @test res == filt!(similar(x), zpk, x)

        # Test with DF2TFilter
        f = DF2TFilter(sos)
        @test res == [filt(f, x[1:50]); filt(f, x[51:end])]
        f = DF2TFilter(tf)
        @test tfres == [filt(f, x[1:50]); filt(f, x[51:end])]
        f = DF2TFilter(zpk)
        @test res == [filt(f, x[1:50]); filt(f, x[51:end])]
    end
end

# Test simple scaling with DF2TFilter
@test filt(DF2TFilter(PolynomialRatio([3.7], [4.2])), x) == x * (3.7/4.2)

#
# filtfilt
#

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

@test_approx_eq_eps zi_python DSP.Filters.filt_stepstate(b, a) 1e-7


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

@test_approx_eq zi_matlab DSP.Filters.filt_stepstate(b, a)


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

@test_approx_eq_eps zi_python DSP.Filters.filt_stepstate(b, a) 1e-7


##############
#
# Filter in place check (with initial conditions)
#
# x = '/Users/rluke/.julia/v0.3/DSP/test/data/spectrogram_x.txt';
# filtcheck  = filter([0.4, 1], [0.9, 0.6], x, 0.4750)
# dlmwrite('filt_check.txt',[filtcheck], 'delimiter', '\t', 'precision', '%.12f')
#
##############

matlab_filt  = readdlm(joinpath(dirname(@__FILE__), "data", "filt_check.txt"),'\t')

a = [0.9, 0.6]
b = [0.4, 1]
z = [0.4750]
x  = readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_x.txt"),'\t')
DSP.Filters.filt!(vec(x), b, a, vec(x), z)

@test_approx_eq matlab_filt x


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

x2_matlab = readdlm(joinpath(dirname(@__FILE__), "data", "filtfilt_output.txt"),'\t')

b = [ 0.00327922,  0.01639608,  0.03279216,  0.03279216,  0.01639608,  0.00327922]
a = [ 1.        , -2.47441617,  2.81100631, -1.70377224,  0.54443269, -0.07231567]
x  = readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_x.txt"),'\t')

@test_approx_eq x2_matlab filtfilt(b, a, x)


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

b = [ 0.00327922,  0.01639608,  0.03279216,  0.03279216,  0.01639608,  0.00327922]
a = [ 1.        , -2.47441617,  2.81100631, -1.70377224,  0.54443269, -0.07231567]

x2_output = readdlm(joinpath(dirname(@__FILE__), "data", "filtfilt_output_2d.txt"),'\t')

x  = readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_x.txt"),'\t')
x  = repmat(x, 1, 3)
x[:,2] = circshift(x[:,2], 64)
x[:,3] = circshift(x[:,3], 128)

@test_approx_eq x2_output filtfilt(b, a, x)


#######################################
#
# Test 2d filtfilt with filter type
#
#
#######################################

b = [ 0.00327922,  0.01639608,  0.03279216,  0.03279216,  0.01639608,  0.00327922]
a = [ 1.        , -2.47441617,  2.81100631, -1.70377224,  0.54443269, -0.07231567]

f = PolynomialRatio(b, a)

# Use 2d data from last test
@test_approx_eq x2_output filtfilt(f, x)


#######################################
#
# Test 2d filtfilt with SecondOrderSections
#
# Our implementation differs from MATLAB, but should match the
# PolynomialRatio provided the filter order is even. (Otherwise
# the extrapolation will differ slightly.)
#######################################

x  = readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_x.txt"),'\t')

f = DSP.digitalfilter(DSP.Lowpass(0.2), DSP.Butterworth(4))
@test_approx_eq filtfilt(convert(SecondOrderSections, f), x) filtfilt(convert(PolynomialRatio, f), x)

f = DSP.digitalfilter(DSP.Highpass(0.1), DSP.Butterworth(6))
@test_approx_eq filtfilt(convert(SecondOrderSections, f), x) filtfilt(convert(PolynomialRatio, f), x)

f = DSP.digitalfilter(DSP.Bandpass(0.1, 0.3), DSP.Butterworth(2))
@test_approx_eq filtfilt(convert(SecondOrderSections, f), x) filtfilt(convert(PolynomialRatio, f), x)

#
# fftfilt/filt
#

for xlen in 2.^(7:18).-1, blen in 2.^(1:6).-1
    b = randn(blen)
    for x in (rand(xlen), rand(xlen, 2))
        filtres = filt(b, [1.0], x)
        fftres = fftfilt(b, x)
        firres = filt(b, x)
        @test_approx_eq filtres fftres
        @test_approx_eq filtres firres
    end
end

# fir_filtfilt

b = randn(10)
for x in (randn(100), randn(100, 2))
    @test_approx_eq DSP.Filters.filtfilt(b, x) DSP.Filters.iir_filtfilt(b, [1.0], x)
    @test_approx_eq DSP.Filters.filtfilt(b, [2.0], x) DSP.Filters.iir_filtfilt(b, [2.0], x)
end
