!(@__DIR__() in LOAD_PATH) && push!(LOAD_PATH, @__DIR__)
using DSP, Test, DelimitedFiles, FilterTestHelpers

@testset "remez_argument_check1" begin
    # bands not monotonically increasing
    @test_throws ArgumentError remez(151, [0, 0.25, 0.25, 0.5], [1.0, 0.0])
    @test_throws ArgumentError remez(151, [0.2, 0.1, 0.25, 0.5], [1.0, 0.0])
end

@testset "remez API equivalence and symmetry" begin
    for n in (20, 21), density in (8, 32)
        bands = [.05, .2, .3, .45]
        desired, weight = [1., 0.], [2., 3.]
        definitions = [(.05, .2) => (1., 2.), (.3, .45) => (0., 3.)]
        h = remez(n, bands, desired; weight, grid_density=density)
        @test h == remez(n, definitions; grid_density=density)
        @test h == reverse(h)
        @test eltype(h) === Float64
        @test length(h) == n
        # Changing frequency units must leave the design unchanged.
        @test h ≈ remez(n, 1000bands, desired; weight, Hz=1000, grid_density=density)
        if iseven(n)
            @test abs(sum(h .* (-1) .^ (0:n-1))) < 1e-14
        end

        h = remez(n, [.05, .45], [1.]; filter_type=filter_type_hilbert,
                  grid_density=density)
        @test h == remez(n, [(.05, .45) => 1.]; neg=true, grid_density=density)
        @test h == -reverse(h)
        @test abs(sum(h)) < 1e-14
        if isodd(n)
            @test h[cld(n, 2)] == 0
            @test abs(sum(h .* (-1) .^ (0:n-1))) < 1e-14
        end

        # Both branches of the differentiator's relative weighting rule.
        h = remez(n, bands, desired; weight, filter_type=filter_type_differentiator,
                  grid_density=density)
        functional = [(.05, .2) => (f -> f, f -> 2 / f),
                      (.3, .45) => (f -> 0f, f -> 3.)]
        @test h == remez(n, functional; neg=true, grid_density=density)
        @test bands == [.05, .2, .3, .45]
        @test desired == [1., 0.]
        @test weight == [2., 3.]
    end
end

@testset "remez heterogeneous responses" begin
    definitions = ((.02, .2) => (f -> 1 + f, f -> 1 + 2f), (.3, .48) => (0., 3.))
    h = remez(32, definitions)
    @test h == remez(32, Any[definitions...])
    @test h == remez(32, [(.02, .2) => (f -> 1 + f, f -> 1 + 2f),
                         (.3, .48) => (Returns(0.), Returns(3.))])
    @test all(isfinite, h)
end

@testset "remez exact grid endpoints" begin
    # For 15 taps and density 4 the step is 1/64. The upper edge replaces
    # the last range point; a band narrower than a step still contributes one.
    seen = Float64[]
    response = f -> (push!(seen, f); 1.)
    bands = [(0., .1) => response, (.2, .201) => response, (.3, .3) => response]
    grid, desired, weight = DSP.Filters._remez_grid(15, bands, 1., 4, false)
    expected = [0., 1/64, 2/64, 3/64, 4/64, 5/64, .1, .201, .3]
    @test seen == expected
    @test grid == cospi.(2expected)
    @test desired == ones(length(expected))
    @test weight == ones(length(expected))
end

@testset "remez iteration limit coefficients" begin
    # Snapshot of the original implementation, not just a warning assertion.
    half = [-0.0057324947447139635, -0.012308085906187437, 0.026588205554288386,
            0.026835890017840345, -0.07953783263891252, -0.03818487882449576,
            0.3086821218293382, 0.5424953165331308]
    expected = [half; reverse(half[1:end-1])]
    h = @test_logs (:warn, r"filter is not converged") remez(15,
        [(0., .2) => 1., (.3, .5) => 0.]; maxiter=1)
    @test h ≈ expected
    h = @test_logs (:warn, r"filter is not converged") remez(15,
        [0., .2, .3, .5], [1., 0.]; maxiter=1)
    @test h ≈ expected
end

@testset "remez weighted equiripple response" begin
    h = remez(31, [(0., .2) => (1., 1.), (.25, .5) => (0., 3.)])
    errors = Float64[]
    for (band, desired, weight) in ((range(0., .2; length=2001), 1., 1.),
                                    (range(.25, .5; length=2501), 0., 3.))
        # Evaluate the zero-phase response directly, independently of the
        # barycentric interpolation used by the implementation.
        response = [sum(h[k] * cospi(2f * (k - 16)) for k in eachindex(h)) for f in band]
        error = weight .* (response .- desired)
        for i in eachindex(error)
            left = i == firstindex(error) || abs(error[i]) >= abs(error[i-1])
            right = i == lastindex(error) || abs(error[i]) >= abs(error[i+1])
            left && right && push!(errors, error[i])
        end
    end
    @test length(errors) == 17
    @test all(errors[i] * errors[i+1] < 0 for i in 1:length(errors)-1)
    # Allow dense-grid discretization error relative to the continuous extrema.
    @test maximum(abs, errors) / minimum(abs, errors) < 1.02
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
    h_scipy = read_reference_data("remez_151_lpf.txt")
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
    h_scipy = read_reference_data("remez_152_lpf.txt")
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
    h_scipy = read_reference_data("remez_51_hpf.txt")
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
h_scipy = read_reference_data("remez_180_bpf.txt")
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
    h_scipy = read_reference_data("remez_20_hilbert.txt")
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
    h_scipy = read_reference_data("remez_21_hilbert.txt")
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
    h_scipy = read_reference_data("remez_200_differentiator.txt")
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
    h_scipy = read_reference_data("remez_201_differentiator.txt")
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
