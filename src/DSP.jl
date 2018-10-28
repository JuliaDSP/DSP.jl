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

if VERSION < v"0.7.0-DEV.3204"
    # N.B. this is not a generally valid replacement, but suffices for DSP, and
    # doesn't interfere with other packages as it is unexported
    mul!(Y, A, B) = A_mul_B!(Y, A, B)
    mul!(Y::AbstractArray, A::AbstractArray, B::Number) = scale!(Y, A, B)
    mul!(Y::AbstractArray, A::Number, B::AbstractArray) = scale!(Y, A, B)
else
    using Compat.LinearAlgebra: mul!
end

if VERSION < v"0.7.0-DEV.3665"
    # As above, this is not a generally valid replacement, but suffices for DSP,
    # and doesn't interfere with other packages as it is unexported
    rmul!(A::AbstractArray, B::Number) = scale!(A, B)
    rmul!(A::AbstractArray, B::Diagonal) = scale!(A, diag(B))
else
    using Compat.LinearAlgebra: rmul!
end

if VERSION >= v"0.7.0-DEV.602"
    if VERSION < v"0.7.0-DEV.986" # JuliaLang/julia#22763
        import Base: conv, conv2, deconv, filt, filt!, xcorr
    else
        export conv, conv2, deconv, filt, filt!, xcorr
    end
    using Compat: copyto!
    import Compat
    include("dspbase.jl")
end

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
