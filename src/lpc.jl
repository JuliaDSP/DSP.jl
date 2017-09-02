module LPC
export lpc, LPCBurg, LPCLevinson
export arburg

# Dispatch types for lpc()
mutable struct LPCBurg; end
mutable struct LPCLevinson; end

"""
`ar_burg(x::AbstractVector, p::Int)`

LPC (Linear-Predictive-Code) estimation, using the Burg method.

Given input signal `x` and prediction order `p`, returns IIR coefficients `a`
and average reconstruction error `prediction_err`. Note that this method does
NOT return the leading `1` present in the true autocorrelative estimate; it
omits it as it is implicit in every LPC estimate, and must be manually
reintroduced if the returned vector should be treated as a polynomial. This
function implements the mathematics published in [1].

[1] - Enhanced Partial Tracking Using Linear Prediction
(DAFX 2003 article, Lagrange et al)
http://www.sylvain-marchand.info/Publications/dafx03.pdf
"""
function arburg{T <: Number}(x::AbstractVector{T}, p::Int)
    ef = x                      # forward error
    eb = x                      # backwards error
    a = [1; zeros(T, p)]        # prediction coefficients
    # Initialize prediction error wtih the variance of the signal
    prediction_err :: T = sum(abs2, x) ./ length(x)
    reflection_coeffs = zeros(T, p)

    for m in 1:p
        efp = ef[1:end-1]
        ebp = eb[2:end]
        k = -2 * dot(ebp,efp) / (sum(abs2, ebp) + sum(abs2,efp))
        ef = efp + k .* ebp
        eb = ebp + k .* efp
        a[1:m+1] = [a[1:m]; 0] + k .* [0; a[m:-1:1]]
        prediction_err *= (1 - k*k)
        reflection_coeffs[m] = k
    end

    return a, prediction_err, reflection_coeffs
end

function lpc{T <: Number}(x::AbstractVector{T}, p::Int, ::LPCBurg)
    a, prediction_err = arburg(x, p)
    return a[2:end], prediction_err
end

"""
lpc_levinson(x::AbstractVector, p::Int)

LPC (Linear-Predictive-Code) estimation, using the Levinson method.

Given input signal `x` and prediction order `p`, returns IIR coefficients `a`
and average reconstruction error `prediction_err`. Note that this method does
NOT return the leading `1` present in the true autocorrelative estimate; it
omits it as it is implicit in every LPC estimate, and must be manually
reintroduced if the returned vector should be treated as a polynomial.  This
function implements the mathematics described in [1].

[1] - The Wiener (RMS) Error Criterion in Filter Design and Prediction
(Studies in Applied Mathematics 1947 article, N. Levison)
"""
function lpc{T <: Number}(x::AbstractVector{T}, p::Int, ::LPCLevinson)
    R_xx = xcorr(x,x)[length(x):end]
    a = zeros(T,p,p)
    prediction_err = zeros(1,p)

    # for m = 1
    a[1,1] = -R_xx[2]/R_xx[1]
    prediction_err[1] = R_xx[1]*(1-a[1,1]^2)

    # for m = 2,3,4,..p
    for m = 2:p
        a[m,m] = (-(R_xx[m+1] + dot(a[m-1,1:m-1],R_xx[m:-1:2]))/prediction_err[m-1])
        a[m,1:m-1] = a[m-1,1:m-1] + a[m,m] * a[m-1,m-1:-1:1]
        prediction_err[m] = prediction_err[m-1]*(1-a[m,m]^2)
    end

    # Return autocorrelation coefficients and error estimate
    a[p,:], prediction_err[p]
end


# Default users to using Burg estimation as it is in general more stable
lpc(x, p) = lpc(x, p, LPCBurg())

end # module LPC
