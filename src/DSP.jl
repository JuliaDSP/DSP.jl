__precompile__()

module DSP

# We want to be very sure we don't pull in Base names unless we're very sure we want them
# This macro will be called in each submodule herein to do the appropriate imports
macro importffts()
    quote
        using AbstractFFTs
        using FFTW
    end
end

@importffts
using Compat.LinearAlgebra: mul!, rmul!



export conv, conv2, deconv, filt, filt!, xcorr

using Compat: copyto!
import Compat
include("dspbase.jl")

include("util.jl")
include("unwrap.jl")
include("windows.jl")
include("periodograms.jl")
include("Filters/Filters.jl")
include("lpc.jl")
include("estimation.jl")

using Reexport
@reexport using .Util, .Windows, .Periodograms, .Filters, .LPC, .Unwrap, .Estimation

include("deprecated.jl")
end
