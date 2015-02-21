require(joinpath(dirname(@__FILE__), "FilterTestHelpers.jl"))
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
            @test_throws ErrorException convert(Biquad, zpk)
        end
        sos = convert(SecondOrderSections, zpk)

        res = filt(sos, x)

        # Test with filt with tf/sos
        @test_approx_eq res filt(tf, x)
        @test_approx_eq res filt!(similar(x), sos, x)
        @test_approx_eq res filt!(similar(x), tf, x)

        # For <= 2 poles, test with biquads
        if n <= 2
            @test_approx_eq res filt(bq, x)
            @test_approx_eq res filt!(similar(x), bq, x)
        end

        # Test that filt with zpk converts
        @test res == filt(zpk, x)
        @test res == filt!(similar(x), zpk, x)
    end
end

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
# Filter initial conditions
# Matlab - Butterworth filter coefficients as SOS
#
#=
  [z, p, k] = butter(5, 0.2)
  [sos, g] = zpk2sos(z, p, k)
  filtfilt(sos, g, x)
=#
#
# Then used breakpoints to extract zi
#
##############

sos = [
   1.0000000000000000e+00   1.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00  -5.0952544949442879e-01   0.0000000000000000e+00
   1.0000000000000000e+00   2.0000000000000000e+00   1.0000000000000000e+00   1.0000000000000000e+00  -1.0965794655679613e+00   3.5544676217239030e-01
   1.0000000000000000e+00   2.0000000000000000e+00   1.0000000000000000e+00   1.0000000000000000e+00  -1.3693171946832927e+00   6.9256913538786335e-01
]
g = 1.2825810789606842e-03

sosfilter = matrix_to_sosfilter(sos, g)
zi_matlab = [0.003947378671810  14.451932524765136  11.374248987589889
             0  -4.492339385234015  -7.570022922409273];

@test_approx_eq zi_matlab DSP.Filters.filt_stepstate(sosfilter)


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

#=
  x = '/home/simon/.julia/DSP/test/data/spectrogram_x.txt'; x = textread(x);
  x2 = filtfilt(sos, g, x)
  dlmwrite('filtfilt_output_2d_sos.txt', x2, 'delimiter', '\t', 'precision', '%.12f')
=#
#
#######################################

x  = readdlm(joinpath(dirname(@__FILE__), "data", "spectrogram_x.txt"),'\t')
x2_output = readdlm(joinpath(dirname(@__FILE__), "data", "filtfilt_output_2d_sos.txt"),'\t')

@test_approx_eq x2_output filtfilt(sosfilter, x)

#
# fftfilt/firfilt
#

for xlen in 2.^(7:18).-1, blen in 2.^(1:6).-1
    b = randn(blen)
    for x in (rand(xlen), rand(xlen, 2))
        filtres = filt(b, [1.0], x)
        fftres = fftfilt(b, x)
        firres = firfilt(b, x)
        @test_approx_eq filtres fftres
        @test_approx_eq filtres firres
    end
end

# fir_filtfilt

b = randn(10)
for x in (randn(100), randn(100, 2))
    @test_approx_eq DSP.Filters.fir_filtfilt(b, x) DSP.Filters.iir_filtfilt(b, [1.0], x)
    @test_approx_eq DSP.Filters.filtfilt(b, [2.0], x) DSP.Filters.iir_filtfilt(b, [2.0], x)
end
