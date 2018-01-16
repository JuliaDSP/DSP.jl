__precompile__()

module DSP

# We want to be very sure we don't pull in Base names unless we're very sure we want them
# This macro will be called in each submodule herein to do the appropriate imports
macro importffts()
    quote
        using AbstractFFTs
        if VERSION >= v"0.7.0-DEV.602"
            using FFTW
        end
    end
end

@importffts

if VERSION >= v"0.7.0-DEV.602"
    if VERSION < v"0.7.0-DEV.986" # JuliaLang/julia#22763
        import Base: conv, conv2, deconv, filt, filt!, xcorr
    else
        export conv, conv2, deconv, filt, filt!, xcorr
    end
    using Compat: copyto!
    include("dspbase.jl")
end

include("util.jl")
include("windows.jl")
include("periodograms.jl")
include("Filters/Filters.jl")
include("lpc.jl")

using Reexport
@reexport using .Util, .Windows, .Periodograms, .Filters, .LPC

include("deprecated.jl")
end
