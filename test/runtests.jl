using Compat, DSP, AbstractFFTs, FFTW, Compat.Test
if VERSION >= v"0.7.0-beta.171"
    using Random: seed!
else
    seed!(x) = srand(x)
end

testfiles = [ "dsp.jl", "util.jl", "windows.jl", "filter_conversion.jl",
    "filter_design.jl", "filter_response.jl", "filt.jl", "filt_stream.jl",
    "periodograms.jl", "resample.jl", "lpc.jl", "estimation.jl", "unwrap.jl",
    "remez_fir.jl" ]

seed!(1776)

for testfile in testfiles
    eval(:(@testset $testfile begin include($testfile) end))
end
