module DSP

include("util.jl")
include("windows.jl")
include("periodograms.jl")
include("Filters/Filters.jl")

using Reexport
@reexport using  .Util, .Windows, .Periodograms, .Filters

include("deprecated.jl")
end
