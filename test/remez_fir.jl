!(dirname(@__FILE__) in LOAD_PATH) && push!(LOAD_PATH, dirname(@__FILE__))
using DSP, Test, DelimitedFiles, FilterTestHelpers

@testset "remez_argument_check1" begin
    # bands not monotonically increasing
    @test_throws ArgumentError remez(151, [0, 0.25, 0.25, 0.5], [1.0, 0.0])
    @test_throws ArgumentError remez(151, [0.2, 0.1, 0.25, 0.5], [1.0, 0.0])
end

@testset "remez_argument_check2" begin
    # bands values out of range
    @test_throws ArgumentError remez(151, [0, 0.23, 0.25, 0.6], [1.0, 0.0])
end

@testset "remez_argument_check3" begin
    # bands values out of range
    @test_throws ArgumentError remez(151, [-0.01, 0.23, 0.25, 0.5], [1.0, 0.0])
end

@testset "remez_argument_check4" begin
    # length of bands not 2x length of desired
    @test_throws ArgumentError remez(151, [0, 0.23, 0.5], [1.0, 0.0])
end

@testset "remez_argument_check5" begin
    # length of bands not 2x length of weight
    @test_throws ArgumentError remez(151, [0, 0.23, 0.25, 0.5], [1.0, 0.0]; weight=[1.0, 1.0, 17.0])
end

#
# Length 151 LPF (Low Pass Filter).
#
@testset "remez_151_lpf" begin
    h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_151_lpf.txt"),'\t')
    h = remez(151, [0, 0.475, 0.5, 1.0], [1.0, 0.0]; Hz=2.0);
    @test h ≈ h_scipy
    h = remez(151, [(0, 0.475) => 1, (0.5, 1.0) => 0]; Hz=2.0);
    @test h ≈ h_scipy
end

#
# Length 152 LPF. Non-default "weight" input.
#
#    from scipy.signal import remez
#    lpf = remez(152, [0, 0.475, 0.5, 1.0], [1.0, 0.0], weight=[1,2], Hz=2.0)
#    lpf.tofile('remez_152_lpf.txt', sep='\n')
#
@testset "remez_152_lpf" begin
    h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_152_lpf.txt"),'\t')
    h = remez(152, [0, 0.475, 0.5, 1.0], [1.0, 0.0]; weight=[1,2], Hz=2.0);
    @test h ≈ h_scipy
    h = remez(152, [(0, 0.475) => (1, 1), (0.5, 1.0) => (0, 2)]; Hz=2.0);
    @test h ≈ h_scipy
end

#
# Length 51 HPF (High Pass Filter).
#
#    from scipy.signal import remez
#    hpf = remez(51, [0, 0.75, 0.8, 1.0], [0.0, 1.0], Hz=2.0)
#    hpf.tofile('remez_51_hpf.txt', sep='\n')
#
@testset "remez_51_hpf" begin
    h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_51_hpf.txt"),'\t')
    h = remez(51, [0, 0.75, 0.8, 1.0], [0.0, 1.0]; Hz=2.0);
    @test h ≈ h_scipy
    h = remez(51, [(0, 0.75) => 0, (0.8, 1.0) => 1]; Hz=2.0);
    @test h ≈ h_scipy
end

#
# Length 180 BPF (Band Pass Filter).
#
#    from scipy.signal import remez
#    bpf = remez(180, [0, 0.375, 0.4, 0.5, 0.525, 1.0], [0.0, 1.0, 0.0], Hz=2.0, maxiter=30)
#    bpf.tofile('remez_180_bpf.txt', sep='\n')
#
@testset "remez_180_bpf" begin
h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_180_bpf.txt"),'\t')
    h = remez(180, [0, 0.375, 0.4, 0.5, 0.525, 1.0], [0.0, 1.0, 0.0]; Hz=2.0, maxiter=30);
    @test h ≈ h_scipy
    h = remez(180, [(0, 0.375) => 0, (0.4, 0.5) => 1, (0.525, 1.0) => 0]; Hz=2.0, maxiter=30);
    @test h ≈ h_scipy
end

@testset "remez_warn_no_converge_after_maxiter_iterations" begin
    @test_logs (:warn, r"filter is not converged") remez(180, [0, 0.375, 0.4, 0.5, 0.525, 1.0], [0.0, 1.0, 0.0]; Hz=2.0)
    @test_logs (:warn, r"filter is not converged") remez(180, [(0, 0.375) => 0, (0.4, 0.5) => 1, (0.525, 1.0) => 0]; Hz=2.0)
end

@testset "remez_error_no_converge_transition_band_too_wide" begin
    @test_throws ErrorException remez(151, [0, 0.1, 0.4, 0.5], [1.0, 0.0])
    @test_throws ErrorException remez(151, [(0, 0.1) => 1, (0.4, 0.5) => 0])
end

#
#  Odd-symmetric filters - hilbert and differentiators type.
#  Even length - much better approximation since it is not constrained to 0 at
#  the nyquist frequency
#
# Length 20 hilbert
#
#    from scipy.signal import remez
#    h = remez(20, [0.1, 0.95], [1], type="hilbert", Hz=2.0)
#    h.tofile('remez_20_hilbert.txt', sep='\n')
#
@testset "remez_20_hilbert" begin
    h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_20_hilbert.txt"),'\t')
    h = remez(20, [0.1, 0.95], [1]; filter_type=filter_type_hilbert, Hz=2.0);
    @test h ≈ h_scipy
    h = remez(20, [(0.1, 0.95) => 1]; neg=true, Hz=2.0);
    @test h ≈ h_scipy
end

#
# Length 21 hilbert
#
#    from scipy.signal import remez
#    h = remez(21, [0.1, 0.95], [1], type="hilbert", Hz=2.0)
#    h.tofile('remez_21_hilbert.txt', sep='\n')
#
@testset "remez_21_hilbert" begin
    h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_21_hilbert.txt"),'\t')
    h = remez(21, [0.1, 0.95], [1]; filter_type=filter_type_hilbert, Hz=2.0);
    @test h ≈ h_scipy
    h = remez(21, [(0.1, 0.95) => 1]; neg=true, Hz=2.0);
    @test h ≈ h_scipy
end

#
# Length 200 differentiator
#
#    from scipy.signal import remez
#    h = remez(200,[0.01, 0.99],[1],type="differentiator" Hz=2.0)
#    h.tofile('remez_200_differentiator.txt', sep='\n')
#
@testset "remez_200_differentiator" begin
    h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_200_differentiator.txt"),'\t')
    h = remez(200, [0.01, 0.99], [1]; filter_type=filter_type_differentiator, Hz=2.0);
    @test h ≈ h_scipy
    h = remez(200, [(0.01, 0.99) => (f -> f/2, f -> 1/f)]; neg=true, Hz=2.0);
    @test h ≈ h_scipy
end

#
# Length 201 differentiator
#
#    from scipy.signal import remez
#    h = remez(201,[0.05, 0.95],[1],type="differentiator" Hz=2.0)
#    h.tofile('remez_201_differentiator.txt', sep='\n')
#
@testset "remez_201_differentiator" begin
    h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_201_differentiator.txt"),'\t')
    h = remez(201, [0.05, 0.95], [1]; filter_type=filter_type_differentiator, Hz=2.0);
    @test h ≈ h_scipy
    h = remez(201, [(0.05, 0.95) => (f -> f/2, f -> 1/f)]; neg=true, Hz=2.0);
    @test h ≈ h_scipy
end


#
# Inverse sinc filter - custom response function
#
@testset "inverse_sinc_response_function" begin
    L = 64

    Fs = 4800*L
    f = range(0, stop=0.5, length=10000)

    P = (π * f * Fs / 4800) ./ sinpi.(f * Fs / 4800)
    Pdb = 20 * log10.(abs.(P))
    Pdb[1] = 0.0

    g_vec = remez(201, [
            (    0.0, 2880.0) => (f -> (f==0) ? 1.0 : abs.((π*f/4800) ./ sinpi.(f/4800)), 1.0),
            (10000.0,  Fs/2) => (0.0, 100.0)
        ]; Hz=Fs)
    g = PolynomialRatio(g_vec, [1.0])
    Gdb = 20*log10.(abs.(freqresp(g, 2π*f)))

    passband_indices = (f*Fs) .< 2880.0
    # Test that maximum passband error is less than 1/4 dB.
    @test maximum(abs.(Pdb[passband_indices] - Gdb[passband_indices])) < 0.25
end

