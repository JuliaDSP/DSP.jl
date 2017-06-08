__precompile__()

module DSP

include("util.jl")
include("windows.jl")
include("periodograms.jl")
include("Filters/Filters.jl")
include("lpc.jl")

using Reexport
@reexport using  .Util, .Windows, .Periodograms, .Filters, .LPC

include("deprecated.jl")
end
