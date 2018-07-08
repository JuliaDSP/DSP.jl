!(dirname(@__FILE__) in LOAD_PATH) && push!(LOAD_PATH, dirname(@__FILE__))
using DSP, Compat, Compat.Test, FilterTestHelpers

#
# Length 151 LPF (Low Pass Filter).
#
@testset "remez_151_lpf" begin
    h = remez_jl2(151, [0 0.475 0.5 1.0], [1.0 0.0]; Hz=2.0);
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
    h = remez_jl2(152, [0 0.475 0.5 1.0], [1.0 0.0]; weight=[1,2], Hz=2.0);
    h_scipy = readdlm(joinpath(dirname(@__FILE__), "data", "remez_152_lpf.txt"),'\t')
    @test h ≈ h_scipy
end

"""


Add non-default "weight" input
Add bandpass, highpass

Odd-symmetric filters - hilbert and differentiators type.
Even length - much better approximation since it is not constrained to 0 at 
the nyquist frequency
h1=remez(20,[0.1, 0.95],[1],type='hilbert', Hz=2.0)
h2=remez(21,[0.1, 0.95],[1],type='hilbert', Hz=2.0)

h=remez(200,[0.01, 0.99],[1],type='differentiator', Hz=2.0)
h=remez(201,[0.05, 0.95],[1],type='differentiator', Hz=2.0)
"""