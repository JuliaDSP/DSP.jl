module LPC


using ..DSP: xcorr

using LinearAlgebra: dot

export lpc, arburg, levinson, LPCBurg, LPCLevinson

# Dispatch types for lpc()
struct LPCBurg end
struct LPCLevinson end

"""
    lpc(x::AbstractVector, p::Integer, [LPCBurg()])

Given input signal `x` and prediction order `p`, returns IIR coefficients `a`
and average reconstruction error `prediction_err`. Note that this method does
NOT return the leading `1` present in the true autocorrelative estimate; it
omits it as it is implicit in every LPC estimate, and must be manually
reintroduced if the returned vector should be treated as a polynomial.

The algorithm used is determined by the last optional parameter, and can be
either `LPCBurg` or `LPCLevinson`.
"""
function lpc end

function lpc(x::AbstractVector{<:Number}, p::Integer, ::LPCBurg)
    a, prediction_err = arburg(x, p)
    popfirst!(a)
    a, prediction_err
end

"""
    arburg(x::AbstractVector, p::Integer)

LPC (Linear-Predictive-Code) estimation, using the Burg method. This function
implements the mathematics published in [1].

[1] - Enhanced Partial Tracking Using Linear Prediction
(DAFX 2003 article, Lagrange et al)
http://www.sylvain-marchand.info/Publications/dafx03.pdf
"""
function arburg(x::AbstractVector{<:Number}, p::Integer)
    # Initialize prediction error with the variance of the signal
    prediction_err = sum(abs2, x) / length(x)
    T = typeof(prediction_err)

    ef = collect(T, x)                  # forward error
    eb = copy(ef)                       # backwards error
    a = zeros(T, p + 1); a[1] = 1       # prediction coefficients
    rev_buf = similar(a, p)             # buffer to store a in reverse
    reflection_coeffs = zeros(T, p)     # reflection coefficients

    @views for m in 1:p
        pop!(ef)
        popfirst!(eb)
        k = -2 * dot(ef, eb) / (sum(abs2, eb) + sum(abs2, ef))
        for i in eachindex(ef, eb)
            ef_i, eb_i = ef[i], eb[i]
            ef[i] += k * eb_i
            eb[i] += k * ef_i
        end
        rev_buf[1:m] = a[m:-1:1]
        @. a[2:m+1] += k * rev_buf[1:m]
        prediction_err *= (1 - k^2)
        reflection_coeffs[m] = k
    end

    return a, prediction_err, reflection_coeffs
end

function lpc(x::AbstractVector{<:Number}, p::Integer, ::LPCLevinson)
    R_xx = xcorr(x,x)[length(x):end]
    a, prediction_err = levinson(R_xx, p)
    a, prediction_err
end

"""
    levinson(x::AbstractVector, p::Integer, LPCLevinson())

LPC (Linear-Predictive-Code) estimation, using the Levinson method. This
function implements the mathematics described in [1].

[1] - The Wiener (RMS) Error Criterion in Filter Design and Prediction
(N. Levinson, Studies in Applied Mathematics 25(1946), 261-278,
https://doi.org/10.1002/sapm1946251261)
"""
function levinson(R_xx::AbstractVector{<:Number}, p::Integer)
    # for m = 1
    a_1 = -R_xx[2] / R_xx[1]
    prediction_err = R_xx[1] * (one(a_1) - a_1^2)
    T = typeof(prediction_err)

    a = zeros(T, p)
    reflection_coeffs = zeros(T, p)
    a[1] = reflection_coeffs[1] = a_1
    rev_buf = similar(a, p - 1)     # buffer to store a in reverse

    @views for m = 2:p
        k = (-(R_xx[m+1] + dot(a[1:m-1], R_xx[m:-1:2])) / prediction_err)
        rev_buf[1:m-1] = a[m-1:-1:1]
        @. a[1:m-1] += k * rev_buf[1:m-1]
        a[m] = reflection_coeffs[m] = k
        prediction_err *= (1 - k^2)
    end

    # Return autocorrelation coefficients, error estimate, and reflection coefficients
    a, prediction_err, reflection_coeffs
end


# Default users to using Burg estimation as it is in general more stable
lpc(x, p) = lpc(x, p, LPCBurg())

end # module LPC
