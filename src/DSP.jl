module DSP
using Base
using Base.FFTW

include("aux.jl")
include("base.jl")
include("windows.jl")
include("periodogram.jl")
include("filter_design.jl")
end
