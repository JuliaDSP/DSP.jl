module LPC


using ..DSP: xcorr

using LinearAlgebra: dot

export lpc, LPCBurg, LPCLevinson

# Dispatch types for lpc()
struct LPCBurg end
struct LPCLevinson end

"""
    lpc(x::AbstractVector, p::Int, [LPCBurg()])

Given input signal `x` and prediction order `p`, returns IIR coefficients `a`
and average reconstruction error `prediction_err`. Note that this method does
NOT return the leading `1` present in the true autocorrelative estimate; it
omits it as it is implicit in every LPC estimate, and must be manually
reintroduced if the returned vector should be treated as a polynomial.

The algorithm used is determined by the last optional parameter, and can be
either `LPCBurg` or `LPCLevinson`.
"""
function lpc end

"""
    lpc(x::AbstractVector, p::Int, LPCBurg())

LPC (Linear-Predictive-Code) estimation, using the Burg method. This function
implements the mathematics published in [1].

[1] - Enhanced Partial Tracking Using Linear Prediction
(DAFX 2003 article, Lagrange et al)
http://www.sylvain-marchand.info/Publications/dafx03.pdf
"""
function lpc(x::AbstractVector{<:Number}, p::Int, ::LPCBurg)
    # Initialize prediction error with the variance of the signal
    prediction_err = sum(abs2, x) / length(x)
    T = typeof(prediction_err)

    ef = collect(T, x)                  # forward error
    eb = copy(ef)                       # backwards error
    a = zeros(T, p + 1); a[1] = 1       # prediction coefficients
    rev_buf = similar(a, p)             # buffer to store a in reverse

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
    end

    popfirst!(a)
    return a, prediction_err
end

"""
    lpc(x::AbstractVector, p::Int, LPCLevinson())

LPC (Linear-Predictive-Code) estimation, using the Levinson method. This
function implements the mathematics described in [1].

[1] - The Wiener (RMS) Error Criterion in Filter Design and Prediction
(N. Levinson, Studies in Applied Mathematics 25(1946), 261-278,
https://doi.org/10.1002/sapm1946251261)
"""
function lpc(x::AbstractVector{<:Number}, p::Int, ::LPCLevinson)
    R_xx = xcorr(x,x)[length(x):end]
    a = zeros(p)

    # for m = 1
    a[1] = -R_xx[2] / R_xx[1]
    prediction_err = R_xx[1] * (1 - a[1]^2)

    # for m = 2,3,4,..p
    @views for m = 2:p
        a[m] = (-(R_xx[m+1] + dot(a[1:m-1], R_xx[m:-1:2])) / prediction_err)
        @. a[1:m-1] += a[m] * a[m-1:-1:1]
        prediction_err *= (1 - a[m]^2)
    end

    # Return autocorrelation coefficients and error estimate
    a, prediction_err
end


# Default users to using Burg estimation as it is in general more stable
lpc(x, p) = lpc(x, p, LPCBurg())

end # module LPC
