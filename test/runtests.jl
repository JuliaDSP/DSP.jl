using DSP, FFTW, Test

using DSP: allocate_output
using Random: seed!

module_tests = [
    :DSPBase => ["dsp.jl", "diric.jl"],
    :Utils => ["util.jl"],
    :Windows => ["windows.jl"],
    :Filters => [
        "filter_conversion.jl", "filter_design.jl", "filter_response.jl", "filt.jl",
        "filt_order.jl", "filt_stream.jl", "resample.jl", "remez_fir.jl"
    ],
    :Periodograms => ["periodograms.jl", "multitaper.jl"],
    :LPC => ["lpc.jl"],
    :Estimation => ["estimation.jl"],
    :Unwrap => ["unwrap.jl"]
]

seed!(1776)

@testset verbose=true "DSP.jl" begin
    @testset verbose=true "$(first(module_test))" for module_test in module_tests
        @testset "$testfile" for testfile in last(module_test)
            include(testfile)
        end
    end
end
