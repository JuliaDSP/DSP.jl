using Compat, DSP, AbstractFFTs, FFTW, Compat.Test

testfiles = [ "dsp.jl", "util.jl", "windows.jl", "filter_conversion.jl",
    "filter_design.jl", "filter_response.jl", "filt.jl", "filt_stream.jl",
    "periodograms.jl", "resample.jl", "lpc.jl"]

for testfile in testfiles
    eval(@__MODULE__, :(@testset $testfile begin include($testfile) end))
end
