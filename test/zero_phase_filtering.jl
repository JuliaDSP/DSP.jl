using DSP, Base.Test


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

@test_approx_eq_eps zi_python DSP.ZeroPhaseFiltering.filt_stepstate(b, a) 1e-7


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

@test_approx_eq zi_matlab DSP.ZeroPhaseFiltering.filt_stepstate(b, a)


##############
#
# Filter check (with initial conditions)
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
DSP.ZeroPhaseFiltering.filt!(b, a, vec(x), z)

@test_approx_eq matlab_filt x


#######################################
#
# Test filtfilt against matlab results
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

@test_approx_eq x2_matlab filtfilt(b, a, vec(x))



