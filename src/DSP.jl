module DSP

using FFTW
using LinearAlgebra: Transpose, mul!, rmul!
using IterTools: subsets

export conv, conv!, deconv, filt, filt!, xcorr

# This function has methods added in `periodograms` but is not exported,
# so we define it here so one can do `DSP.allocate_output` instead of
# `DSP.Periodograms.allocate_output`.
function allocate_output end

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
@reexport using .Util, .Windows, .Periodograms, .Filters, .LPC, .Unwrap, .Estimation

include("deprecated.jl")

# Compile the convolution and filtering paths a test suite or first use pays
# for otherwise: `conv` specializes per element type, dimensionality and
# algorithm (the overlap-save edge kernels alone are ~24 method instances).
using PrecompileTools: @setup_workload, @compile_workload
@setup_workload begin
    @compile_workload begin
        for T in (Float64, Float32, ComplexF64)
            u1 = ones(T, 16); v1 = ones(T, 5)
            u2 = ones(T, 8, 8); v2 = ones(T, 3, 3)
            u3 = ones(T, 6, 6, 6); v3 = ones(T, 2, 2, 2)
            for algorithm in (:direct, :fft_simple, :fft_overlapsave)
                conv(u1, v1; algorithm)
                conv(u2, v2; algorithm)
                conv(u3, v3; algorithm)
            end
        end
        filt([0.5, 0.5], [1.0, 0.1], ones(16))
    end
end
end
