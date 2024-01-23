module Estimation

using LinearAlgebra: eigen, svd
using FFTW
using Statistics: mean

export esprit, jacobsen, quinn

"""
    esprit(x::AbstractArray, M::Integer, p::Integer, Fs::Real=1.0)

ESPRIT [^Roy1986] algorithm for frequency estimation.
Estimation of Signal Parameters via Rotational Invariance Techniques

Given length N signal "x" that is the sum of p sinusoids of unknown frequencies,
estimate and return an array of the p frequencies.

# Arguments
- `x::AbstractArray`: complex length N signal array
- `M::Integer`: size of correlation matrix, must be <= N.
      The signal subspace is computed from the SVD of an M x (N-M+1) signal matrix
      formed from N-M+1 length-M shifts of the signal x in its columns.
      For best performance for 1 sinusoid, use M = (N+1)/3 (according to van der Veen and Leus).
      For faster execution (due to smaller SVD), use small M or small N-M
- `p::Integer`: number of sinusoids to estimate.
- `Fs::Float64`: sampling frequency, in Hz.

# Returns
length p real array of frequencies in units of Hz.

[^Roy1986]: R Roy, A Paulraj and T Kailath, ESPRIT - A subspace approach to estimation of parameters of cisoids in noise, IEEE Trans. Acoustics, Speech, Signal Process., 34, 1340-1342 (1986). [url](http://ieeexplore.ieee.org/abstract/document/1164935/).
"""
function esprit(x::AbstractArray, M::Integer, p::Integer, Fs::Real=1.0)
    count(v->v != 1, size(x)) <= 1 || error("`x` must be a 1D signal")
    N = length(x)
    X = x[ (1:M) .+ (0:N-M)' ]
    U,s,V = svd(X)
    D,_ = eigen( U[1:end-1,1:p] \ U[2:end,1:p] )

    angle.(D)*Fs/2π
end

"""
    jacobsen(x::Vector, Fs::Real = 1.0)

Estimate the largest frequency (in Hz) in the signal `x` using Jacobsen's
algorithm [^Jacobsen2007].

[^Jacobsen2007]: E Jacobsen and P Kootsookos, "Fast, Accurate Frequency Estimators", Chapter
10 in "Streamlining Digital Signal Processing", edited by R. Lyons, 2007, IEEE Press.
"""
function jacobsen(x::Vector{<:Real}, Fs::Real = 1.0)
    N = length(x)
    X = rfft(x)
    k = argmax(abs.(X))  # index of DFT peak
    if (k+1 <= N) && (k-1 >= 1)
        δ = -real((X[k+1] - X[k-1]) / (2X[k] - X[k-1] - X[k+1]) )
    else
        δ = 0.0
    end

    (k + δ)*Fs/N
end

function jacobsen(x::Vector{<:Complex}, Fs::Real = 1.0)
    N = length(x)
    X = fftshift(fft(x))
    k = argmax(abs.(X))  # index of DFT peak
    f = fftshift(fftfreq(N, Fs))[k]  # frequency at index k
    if (k+1 <= N) && (k-1 >= 1)
        δ = -real((X[k+1] - X[k-1]) / (2X[k] - X[k-1] - X[k+1]) )
    else
        δ = 0.0
    end

    f + δ*Fs/N
end

"""
    quinn(x::Vector, f0::Real = 0.0, Fs::Real = 1.0 ; tol = 1e-6, maxiters = 20)

    quinn(x::Vector, Fs::Real = 1.0 ; kwargs...)

    quinn(x::Vector ; kwargs...)

Algorithms by Quinn and Quinn & Fernandes for frequency estimation. Given a
signal `x` and an initial guess `f0`, estimate and return the frequency of the
largest sinusoid in `x`. `Fs` is the sampling frequency. All frequencies are
expressed in Hz.

If the initial guess `f0` is not provided or if it is equal to zero, then a guess
is calculated by `quinn` using Jacobsen's estimator. The sampling frequency `Fs`
defaults to `1.0`.

The following keyword arguments control the algorithm's behavior:

- `tol`: () the algorithm stops when the absolut value of the
  difference between two consecutive estimates is less than `tol`. Defaults to
  `1e-6`.
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
quinn(x ; kwargs...) = quinn(x, 0.0, 1.0 ; kwargs...)

quinn(x, Fs ; kwargs...) = quinn(x, 0.0, Fs ; kwargs...)

function quinn(x::Vector{<:Real}, f0::Real, Fs::Real ; tol = 1e-6, maxiters = 20)
    fₙ = Fs/2
    T = length(x)

    # Run a quick estimate of largest sinusoid in x
    if f0 == 0.0
        f_est = jacobsen(x, Fs)
        ω̂ = π*f_est/fₙ  # estimate ranges from 0 to π
    else
        ω̂ = π*f0/fₙ
    end

    # remove DC
    x .= x .- mean(x)

    # Initialize algorithm
    α = 2cos(ω̂)

    # iteration
    ξ = zeros(eltype(x), T)
    β = zero(eltype(x))
    iter = 0
    @inbounds while iter < maxiters
        iter += 1
        ξ[1] = x[1]
        ξ[2] = x[2] + α*ξ[1]
        for t in 3:T
            ξ[t] = x[t] + α*ξ[t-1] - ξ[t-2]
        end
        β = ξ[2]/ξ[1]
        for t = 3:T
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
    T = length(x)

    # Run a quick estimate of largest sinusoid in x
    if f0 == 0.0
        f_est = jacobsen(x, Fs)
        ω̂ = π*f_est/fₙ
    else
        ω̂ = π*f0/fₙ
    end

    # Remove any DC term in x
    x .= x .- complex(mean(real.(x)), mean(imag.(x)))

    # iteration
    ξ = zeros(eltype(x), T)
    iter = 0
    @inbounds while iter < maxiters
        iter += 1
        # step 2
        ξ[1] = x[1]
        for t in 2:T
            ξ[t] = x[t] + exp(complex(0,ω̂))*ξ[t-1]
        end
        # step 3
        S = let s = 0.0
                for t=2:T
                    s += x[t]*conj(ξ[t-1])
                end
                s
            end
        num = imag(S*exp(complex(0,-ω̂)))
        den = sum(abs2.(ξ[1:end-1]))
        ω̂ += 2*num/den

        # stop condition
        (abs(2*num/den) < tol) && break
    end

    fₙ*ω̂/π, iter == maxiters
end

end # end module definition
