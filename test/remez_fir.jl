!(dirname(@__FILE__) in LOAD_PATH) && push!(LOAD_PATH, dirname(@__FILE__))
using DSP, Compat, Compat.Test, FilterTestHelpers

@testset "remez_argument_check1" begin
    # bands not monotonically increasing
    @test_throws ArgumentError remez(151, [0, 0.25, 0.23, 0.5], [1.0, 0.0])
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
    h = remez(151, [0, 0.475, 0.5, 1.0], [1.0, 0.0]; Hz=2.0);
    h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_151_lpf.txt"),'\t')
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
    h = remez(152, [0, 0.475, 0.5, 1.0], [1.0, 0.0]; weight=[1,2], Hz=2.0);
    h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_152_lpf.txt"),'\t')
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
    h = remez(51, [0, 0.75, 0.8, 1.0], [0.0, 1.0]; Hz=2.0);
    h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_51_hpf.txt"),'\t')
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
    h = remez(180, [0, 0.375, 0.4, 0.5, 0.525, 1.0], [0.0, 1.0, 0.0]; Hz=2.0, maxiter=30);
    h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_180_bpf.txt"),'\t')
    @test h ≈ h_scipy
end

@testset "remez_warn_no_converge_after_maxiter_iterations" begin
    @static if VERSION ≥ v"0.7.0-DEV.2988"
        @test_logs (:warn, r"filter is not converged") remez(180, [0, 0.375, 0.4, 0.5, 0.525, 1.0], [0.0, 1.0, 0.0]; Hz=2.0)
    else
        warn_check(msg) = contains(msg, "filter is not converged")
        @test_warn warn_check remez(180, [0, 0.375, 0.4, 0.5, 0.525, 1.0], [0.0, 1.0, 0.0]; Hz=2.0)
    end
end

@testset "remez_error_no_converge_transition_band_too_wide" begin
    @test_throws ErrorException remez(151, [0, 0.1, 0.4, 0.5], [1.0, 0.0])
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
    h = remez(20, [0.1, 0.95],[1]; filter_type=filter_type_hilbert, Hz=2.0);
    h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_20_hilbert.txt"),'\t')
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
    h = remez(21, [0.1, 0.95],[1]; filter_type=filter_type_hilbert, Hz=2.0);
    h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_21_hilbert.txt"),'\t')
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
    h = remez(200, [0.01, 0.99],[1]; filter_type=filter_type_differentiator, Hz=2.0);
    h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_200_differentiator.txt"),'\t')
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
    h = remez(201, [0.05, 0.95],[1]; filter_type=filter_type_differentiator, Hz=2.0);
    h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_201_differentiator.txt"),'\t')
    @test h ≈ h_scipy
end
