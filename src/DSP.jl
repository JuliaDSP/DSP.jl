module DSP

using LinearAlgebra: mul!

export filt, filt!

# This function has methods added in `periodograms` but is not exported,
# so we define it here so one can do `DSP.allocate_output` instead of
# `DSP.Periodograms.allocate_output`.
function allocate_output end

include("convolutions.jl")

include("dspbase.jl")
include("util.jl")

include("unwrap.jl")
include("windows.jl")
include("periodograms.jl")
include("Filters/Filters.jl")
include("lpc.jl")
include("estimation.jl")
include("diric.jl")

using Reexport
@reexport using .Util, .Windows, .Periodograms, .Filters, .LPC, .Unwrap, .Estimation, .Convolutions

include("deprecated.jl")
end
