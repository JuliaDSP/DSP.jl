using  DSP, FFTW, Test
import DSP: Frequencies, fftfreq, rfftfreq

using Random: seed!

testfiles = [ "dsp.jl", "util.jl", "windows.jl", "filter_conversion.jl",
    "filter_design.jl", "filter_response.jl", "filt.jl", "filt_stream.jl",
    "periodograms.jl", "resample.jl", "lpc.jl", "estimation.jl", "unwrap.jl",
    "remez_fir.jl" ]

seed!(1776)

for testfile in testfiles
    eval(:(@testset $testfile begin include($testfile) end))
end
