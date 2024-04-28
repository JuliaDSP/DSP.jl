module Estimation

using LinearAlgebra: eigen, svd
using FFTW
using Statistics: mean

export esprit, jacobsen, quinn

"""
    esprit(x::AbstractArray, M::Integer, p::Integer, Fs::Real=1.0)

ESPRIT [^Roy1986] algorithm for frequency estimation.
Estimation of Signal Parameters via Rotational Invariance Techniques

Given a signal `x` of length N that is the sum of `p` complex exponentials (cisoids) with unknown constant amplitudes and phases,
estimate and return an array of the `p` frequencies.

```math
x(t)=\\sum_{i=1}^p A_ie^{j2\\pi f_i t}
```
where ``A_i`` and ``f_i`` are the amplitude and frequency of the ``i``th cisoid.

# Arguments
- `x::AbstractArray`: complex length N signal array
- `M::Integer`: size of correlation matrix, must be <= N.
      The signal subspace is computed from the SVD of an M x (N-M+1) signal matrix
      formed from N-M+1 length-M shifts of the signal x in its columns.
      For best performance for 1 cisoid, use M = (N+1)/3 (according to van der Veen and Leus).
      For faster execution (due to smaller SVD), use small M or small N-M
- `p::Integer`: number of complex exponentials/frequencies to estimate.
- `Fs::Float64`: sampling frequency, in Hz.

# Returns
length `p` real array of frequencies in units of Hz.

[^Roy1986]: R Roy, A Paulraj and T Kailath, ESPRIT -
    A subspace approach to estimation of parameters of cisoids in noise,
    IEEE Trans. Acoustics, Speech, Signal Process., 34, 1340-1342 (1986).
    <http://ieeexplore.ieee.org/abstract/document/1164935/>.

# Examples
Suppose a signal ``x(t)=2e^{j2\\pi(2500)t} + 5e^{j2\\pi(400)t}`` is corrupted by noise. The goal is
to retrieve the frequencies in the signal. Let the sampling frequency be ``8kHz``.

```@example
julia> Fs = 8000;                 # sampling frequency (Hz)
julia> true_freq = [2500 400];    # true frequencies (Hz)
julia> t = collect(1:Fs)/Fs;      # time vector
julia> x = 2exp.(1im*2π*true_freq[1]*t) + 5exp.(1im*2π*true_freq[2]*t); # original signal x
```
Add gaussian noise in complex form to the signal ``x`` to mimic noise corruption.
```@example
julia> noise = randn(Fs, 2)*[1; 1im];    # complex random noise 
julia> x += noise;                       # add noise to signal x
```
Run the ESPRIT algorithm to retrieve approximate frequencies.
```@example
julia> M = 5;    # window length of 5
julia> p = 2;    # 2 frequencies to estimate
julia> esprit(x, M, p, Fs)
2-element Vector{Float64}:
2502.2704154274957
399.8420984461712
```

"""
function esprit(x::AbstractArray, M::Integer, p::Integer, Fs::Real=1.0)
    count(!isone, size(x)) <= 1 || throw(ArgumentError("`x` must be a 1D signal"))
    N = length(x)
    X = x[ (1:M) .+ (0:N-M)' ]
    U, _ = svd(X)
    D, _ = eigen( U[1:end-1, 1:p] \ U[2:end, 1:p] )

    angle.(D)*Fs/2π
end

"""
    jacobsen(x::AbstractVector, Fs::Real = 1.0)

Estimate the largest frequency in the complex signal `x` using Jacobsen's
algorithm [^Jacobsen2007]. Argument `Fs` is the sampling frequency.  All
frequencies are expressed in Hz.

If the signal `x` is real, the estimated frequency is guaranteed to be
positive, but it may be highly inaccurate (especially for frequencies close to
zero or to `Fs/2`).

If the sampling frequency `Fs` is not provided, then it is assumed that `Fs =
1.0`.

[^Jacobsen2007]: E Jacobsen and P Kootsookos, "Fast, Accurate Frequency Estimators", Chapter
    10 in "Streamlining Digital Signal Processing", edited by R. Lyons, 2007, IEEE Press.
"""
function jacobsen(x::AbstractVector, Fs::Real = 1.0)
    N = length(x)
    X = fft(x)
    k = argmax(abs.(X)) # index of DFT peak
    fpeak = fftfreq(N, Fs)[k]  # peak frequency
    if (k == N)     # k+1 is OOB -- X[k+1] = X[1]
        Xkm1 = X[N-1]
        Xkp1 = X[1]
    elseif (k == 1) # k-1 is OOB -- X[k-1] = X[N]
        Xkp1 = X[2]
        Xkm1 = X[N]
    else
        Xkp1 = X[k+1]
        Xkm1 = X[k-1]
    end
    δ = -real( (Xkp1 - Xkm1) / (2X[k] - Xkm1 - Xkp1) )
    estimate = fpeak + δ*Fs/N
    # if signal is real, return positive frequency
    if eltype(x) <: Real
        return abs(estimate)
    end
    return estimate
end

"""
    quinn(x::Vector, f0::Real, Fs::Real = 1.0 ; tol = 1e-6, maxiters = 20)

    quinn(x::Vector, Fs::Real = 1.0 ; kwargs...)

    quinn(x::Vector ; kwargs...)

Algorithms by Quinn and Quinn & Fernandes for frequency estimation. Given a
signal `x` and an initial guess `f0`, estimate and return the frequency of the
largest sinusoid in `x`. `Fs` is the sampling frequency. All frequencies are
expressed in Hz.

If the initial guess `f0` is not provided, then a guess is calculated using
Jacobsen's estimator. If the sampling frequency `Fs` is not provided, then it
is assumed that `Fs = 1.0`.

The following keyword arguments control the algorithm's behavior:

- `tol`: the algorithm stops when the absolut value of the difference between
   two consecutive estimates is less than `tol`. Defaults to `1e-6`.
- `maxiters`: the maximum number of iterations to run. Defaults to `20`.

Returns a tuple `(estimate, reachedmaxiters)`, where `estimate` is the
estimated frequency, and `reachedmaxiters` is `true` if the algorithm finished
after running for `maxiters` iterations (this may indicate that the algorithm
did not converge).

If the signal `x` is real, then the algorithm used is [^Quinn1991]. If the signal is
complex, the algorithm is [^Quinn2009].

[^Quinn1991]: B Quinn and J Fernandes, "A fast efficient technique for the
    estimation of frequency", Biometrika, Vol. 78 (1991).

[^Quinn2009]: B Quinn, "Recent advances in rapid frequency estimation", Digital
    Signal Processing, Vol. 19 (2009), Elsevier.
"""
quinn(x ; kwargs...) = quinn(x, jacobsen(x, 1.0), 1.0 ; kwargs...)

quinn(x, Fs ; kwargs...) = quinn(x, jacobsen(x, Fs), Fs ; kwargs...)

function quinn(x::Vector{<:Real}, f0::Real, Fs::Real ; tol = 1e-6, maxiters = 20)
    fₙ = Fs/2
    N = length(x)

    # Run a quick estimate of largest sinusoid in x
    ω̂ = π*f0/fₙ

    # remove DC
    x .= x .- mean(x)

    # Initialize algorithm
    α = 2cos(ω̂)

    # iteration
    ξ = zeros(eltype(x), N)
    β = zero(eltype(x))
    iter = 0
    @inbounds while iter < maxiters
        iter += 1
        ξ[1] = x[1]
        ξ[2] = x[2] + α*ξ[1]
        for t in 3:N
            ξ[t] = x[t] + α*ξ[t-1] - ξ[t-2]
        end
        β = ξ[2]/ξ[1]
        for t = 3:N
            β += (ξ[t]+ξ[t-2])*ξ[t-1]
        end
        β = β/(ξ[1:end-1]'*ξ[1:end-1])
        abs(α - β) < tol && break
        α = 2β-α
    end

    fₙ*acos(0.5*β)/π, iter == maxiters
end

function quinn(x::Vector{<:Complex}, f0::Real, Fs::Real ; tol = 1e-6, maxiters = 20)
    fₙ = Fs/2
    N = length(x)

    ω̂ = π*f0/fₙ

    # Remove any DC term in x
    x .= x .- mean(x)

    # iteration
    ξ = zeros(eltype(x), N)
    iter = 0
    @inbounds while iter < maxiters
        iter += 1
        # step 2
        ξ[1] = x[1]
        for t in 2:N
            ξ[t] = x[t] + exp(complex(0,ω̂))*ξ[t-1]
        end
        # step 3
        S = let s = 0.0
                for t=2:N
                    s += x[t]*conj(ξ[t-1])
                end
                s
            end
        num = imag(S*cis(-ω̂))
        den = sum(abs2.(ξ[1:end-1]))
        ω̂ += 2*num/den

        # stop condition
        (abs(2*num/den) < tol) && break
    end

    fₙ*ω̂/π, iter == maxiters
end

end # end module definition
