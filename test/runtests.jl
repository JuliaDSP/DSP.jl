using DSP, FFTW, Test

using DSP: allocate_output
using Random: seed!

testfiles = ["dsp.jl", "util.jl", "windows.jl", "filter_conversion.jl",
    "diric.jl",
    "filter_design.jl", "filter_response.jl", "filt.jl", "filt_stream.jl",
    "periodograms.jl", "multitaper.jl", "resample.jl", "lpc.jl", "estimation.jl", "unwrap.jl",
    "remez_fir.jl" ]

seed!(1776)

for testfile in testfiles
    @testset "$testfile" begin
        time = @elapsed eval(:(include($testfile)))
        @show time
    end
end
