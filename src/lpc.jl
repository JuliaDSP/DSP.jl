module LPC


using ..DSP: xcorr

using LinearAlgebra: dot, BlasComplex, BLAS

export lpc, arburg, levinson, LPCBurg, LPCLevinson

# Dispatch types for lpc()
struct LPCBurg end
struct LPCLevinson end

"""
    lpc(x::AbstractVector, p::Integer, [LPCBurg()])

Given input signal `x` and prediction order `p`, returns IIR coefficients `a`
and average reconstruction error `prediction_err`. Note that this method does
NOT return the leading ``1`` present in the true autocorrelative estimate; it
omits it as it is implicit in every LPC estimate, and must be manually
reintroduced if the returned vector should be treated as a polynomial.

The algorithm used is determined by the last optional parameter, and can be
either `LPCBurg` ([`arburg`](@ref)) or `LPCLevinson` ([`levinson`](@ref)).
"""
function lpc end

function lpc(x::AbstractVector{<:Number}, p::Integer, ::LPCBurg)
    a, prediction_err = arburg(x, p)
    popfirst!(a)
    a, prediction_err
end

"""
    arburg(x::AbstractVector, p::Integer)

LPC (Linear Predictive Coding) estimation, using the Burg method.
This function implements the mathematics published in [^Lagrange], and
the recursion relation as noted in [^Vos], in turn referenced from [^Andersen].

[^Lagrange]: Enhanced Partial Tracking Using Linear Prediction.
    [DAFX 2003 article, Lagrange et al]
    (http://www.sylvain-marchand.info/Publications/dafx03.pdf)

[^Vos]: [A Fast Implementation of Burg’s Method]
    (https://www.opus-codec.org/docs/vos_fastburg.pdf).
    © 2013 Koen Vos, licensed under [CC BY 3.0]
    (https://creativecommons.org/licenses/by/3.0/)

[^Andersen]: N. Andersen. Comments on the performance of maximum entropy algorithms.
    Proceedings of the IEEE 66.11: 1581-1582, 1978.
"""
function arburg(x::AbstractVector{T}, p::Integer) where T<:Number
    # Initialize prediction error with the variance of the signal
    prediction_err = sum(abs2, x) / length(x)
    R = typeof(prediction_err)
    F = promote_type(R, Base.promote_union(T))

    ef = collect(F, x)                  # forward error
    eb = copy(ef)                       # backwards error
    a = zeros(F, p + 1); a[1] = 1       # prediction coefficients
    rev_buf = similar(a, p)             # buffer to store a in reverse
    reflection_coeffs = similar(a, p)   # reflection coefficients

    den = 2sum(abs2, ef)
    ratio = one(R)

    @views for m in 1:p
        cf = pop!(ef)
        cb = popfirst!(eb)
        den = ratio * den - abs2(cf) - abs2(cb)

        k = -2 * dot(eb, ef) / den
        reflection_coeffs[m] = k

        copyto!(rev_buf, CartesianIndices((1:m,)), a, CartesianIndices((m:-1:1,)))
        @. a[2:m+1] += k * conj(rev_buf[1:m])

        # update prediction errors
        for i in eachindex(ef, eb)
            ef_i, eb_i = ef[i], eb[i]
            ef[i] += k * eb_i
            eb[i] += conj(k) * ef_i
        end

        ratio = one(R) - abs2(k)
        prediction_err *= ratio
    end

    return conj!(a), prediction_err, reflection_coeffs
end

function lpc(x::AbstractVector{<:Number}, p::Integer, ::LPCLevinson)
    R_xx = xcorr(x; scaling=:biased)[length(x):end]
    a, prediction_err = levinson(R_xx, p)
    a, prediction_err
end

"""
    levinson(x::AbstractVector, p::Integer)

Implements Levinson recursion, as described in [^Levinson], to find
the solution `a` of the linear equation
```math
\\mathbf{T} (-\\vec{a})
=
\\begin{bmatrix}
    x_2 \\\\
    \\vdots \\\\
    x_{p+1}
\\end{bmatrix}
```
in the case where ``\\mathbf{T}`` is Hermitian and Toeplitz, with first column `x[1:p]`.
This function can be used for LPC (Linear Predictive Coding) estimation,
by providing `LPCLevinson()` as an argument to `lpc`.

[^Levinson]: The Wiener (RMS) Error Criterion in Filter Design and Prediction.
    N. Levinson, Studies in Applied Mathematics 25(1946), 261-278.\\
    <https://doi.org/10.1002/sapm1946251261>
"""
function levinson(R_xx::AbstractVector{U}, p::Integer) where U<:Number
    # for m = 1
    k = -R_xx[2] / R_xx[1]
    F = promote_type(Base.promote_union(U), typeof(k))
    prediction_err = real(R_xx[1] * (one(F) - abs2(k)))
    R = typeof(prediction_err)
    T = promote_type(F, R)

    a = zeros(T, p)
    reflection_coeffs = zeros(T, p)
    a[1] = reflection_coeffs[1] = k
    rev_a = similar(a, p - 1)     # buffer to store a in reverse

    @views for m = 2:p
        copyto!(rev_a, CartesianIndices((1:m-1,)), a, CartesianIndices((m-1:-1:1,)))
        k = -(R_xx[m+1] + dotu(R_xx[2:m], rev_a[1:m-1])) / prediction_err
        @. a[1:m-1] += k * conj(rev_a[1:m-1])
        a[m] = reflection_coeffs[m] = k
        prediction_err *= (one(R) - abs2(k))
    end

    # Return autocorrelation coefficients, error estimate, and reflection coefficients
    a, prediction_err, reflection_coeffs
end

# for convenience, define dotu as dot for real vectors
dotu(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}) = dot(x, y)
dotu(x::AbstractVector{T}, y::AbstractVector{T}) where T<:BlasComplex = BLAS.dotu(x, y)
function dotu(x::AbstractVector{T}, y::AbstractVector{V}) where {T,V}
    dotprod = zero(promote_type(T, V))
    for i in eachindex(x, y)
        dotprod += x[i] * y[i]
    end
    dotprod
end

# Default users to using Burg estimation as it is in general more stable
lpc(x, p) = lpc(x, p, LPCBurg())

end # module LPC
