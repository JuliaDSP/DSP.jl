module FFTAImpl

using FFTA
using ..FFTBackend
import ..FFTBackend: AbstractFFTBackend, fft, ifft, ifft!, bfft, rfft, irfft, brfft,
                     plan_fft, plan_ifft, plan_bfft, plan_rfft, plan_irfft, plan_brfft,
                     plan_fft!, plan_bfft!

export FFTABackend

"""
    FFTABackend <: AbstractFFTBackend

FFT backend using FFTA.jl. This is the default backend and works with arbitrary numeric types
including BigFloat. For better performance with Float32/Float64, use FFTWBackend by loading FFTW.jl.
"""
struct FFTABackend <: AbstractFFTBackend end

# Transform implementations
fft(::FFTABackend, x) = FFTA.fft(x)
bfft(::FFTABackend, x) = FFTA.bfft(x)

function ifft(::FFTABackend, x)
    result = FFTA.bfft(x)
    return result ./ length(x)
end

function ifft!(::FFTABackend, x)
    # FFTA doesn't have in-place operations, so we compute and copy back
    result = FFTA.bfft(x) ./ length(x)
    copyto!(x, result)
    return x
end

function rfft(::FFTABackend, x::AbstractVector{<:Real})
    X = FFTA.fft(complex.(x))
    return X[1:length(x)÷2+1]
end

function rfft(::FFTABackend, x::AbstractArray{<:Real})
    # For arrays, transform along first dimension
    X = FFTA.fft(complex.(x))
    return X[1:size(x, 1)÷2+1, axes(X)[2:end]...]
end

function irfft(::FFTABackend, X::AbstractVector{<:Complex}, d::Int)
    if iseven(d)
        Xfull = vcat(X, conj.(X[end-1:-1:2]))
    else
        Xfull = vcat(X, conj.(X[end:-1:2]))
    end
    return real.(FFTA.bfft(Xfull)) ./ d
end

function brfft(::FFTABackend, X::AbstractVector{<:Complex}, d::Int)
    if iseven(d)
        Xfull = vcat(X, conj.(X[end-1:-1:2]))
    else
        Xfull = vcat(X, conj.(X[end:-1:2]))
    end
    return real.(FFTA.bfft(Xfull))
end

# Plan wrapper (FFTA doesn't have true planning, so we simulate it)
struct FFTAPlan{B,F,S}
    backend::B
    f::F
    sz::S
end

(p::FFTAPlan)(x) = p.f(p.backend, x)
Base.size(p::FFTAPlan) = p.sz

# LinearAlgebra.mul! support for plans
import LinearAlgebra: mul!
function mul!(Y, p::FFTAPlan, X)
    result = p(X)
    copyto!(Y, result)
    return Y
end

plan_fft(b::FFTABackend, x; kw...) = FFTAPlan(b, fft, size(x))
plan_ifft(b::FFTABackend, x; kw...) = FFTAPlan(b, ifft, size(x))
plan_bfft(b::FFTABackend, x; kw...) = FFTAPlan(b, bfft, size(x))
plan_rfft(b::FFTABackend, x; kw...) = FFTAPlan(b, rfft, size(x))
plan_irfft(b::FFTABackend, x, d; kw...) = FFTAPlan(b, (b, X) -> irfft(b, X, d), size(x))
plan_brfft(b::FFTABackend, x, d; kw...) = FFTAPlan(b, (b, X) -> brfft(b, X, d), size(x))

# In-place variants (not truly in-place for FFTA, but API compatible)
plan_fft!(b::FFTABackend, x; kw...) = FFTAPlan(b, (b, x) -> copyto!(x, fft(b, x)), size(x))
plan_bfft!(b::FFTABackend, x; kw...) = FFTAPlan(b, (b, x) -> copyto!(x, bfft(b, x)), size(x))

# Plan inversion
struct InversePlan{P}
    p::P
end

function Base.inv(p::FFTAPlan)
    if p.f === fft
        return InversePlan(FFTAPlan(p.backend, ifft, p.sz))
    elseif p.f === bfft
        # inv of bfft plan would need to track normalization
        return InversePlan(FFTAPlan(p.backend, (b, x) -> fft(b, x) ./ prod(p.sz), p.sz))
    else
        error("Inverse not defined for this plan type")
    end
end

# Set FFTA as default on module load
function __init__()
    FFTBackend.set_fft_backend!(FFTABackend())
end

end # module
