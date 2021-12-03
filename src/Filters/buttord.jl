
#
# Butterworth prototype filter transformations
#

_toprototype(Wp::Real, Ws::Real, ftype::Type{Lowpass}) = Ws / Wp
_toprototype(Wp::Real, Ws::Real, ftype::Type{Highpass}) = Wp / Ws

function _toprototype(Wp::Vector{Real}, Ws::Vector{Real}, ftype::Type{Bandpass}) 
    # bandpass filter must have two corner frequencies we're computing with
    length(Wp) == 2 || error("2 passband frequencies were expected.")
    length(Ws) == 2 || error("2 stopband frequencies were expected.")
    Wn = (Ws.^2 .- Wp[1]*Wp[2]) ./ (Ws .* (Wp[1]-Wp[2]))
    Wn
end

_fromprototype(Wp::Real, Wscale::Real, ftype::Type{Lowpass}) = Wp * Wscale
_fromprototype(Wp::Real, Wscale::Real, ftype::Type{Highpass}) = Wp / Wscale

function _fromprototype(Wp::Vector{Real}, Wscale::Real, ftype::Type{Bandpass})
    length(Wp) == 2 || error("2 passband frequencies were expected.")
    Wsc = [-Wscale, Wscale]
    Wn = -Wsc .* (Wp[2]-Wp[1])./2 .+ sqrt.( Wsc.^2/4*(Wp[2]-Wp[1]).^2 + Wp[1]*Wp[2])
    sort!(Wn)
end

_order_estimate(Rp::Real, Rs::Real, warp::Real) = ceil((log10(^(10, 0.1*Rs) - 1) - log10(^(10, 0.1*Rp) - 1)) / (2*log10(warp)))
_natfreq_estimate(warp::Real, Rs::Real, order::Integer) = warp / (^(10, 0.1*Rs) - 1)^(1/(2*order))

function buttord(Wp::Array{Real}, Ws::Array{Real}, Rs::Real, Rp::Real, ftype::Union{Type{Bandpass}, Type{Bandstop}})
    
    # pre-warp both components
    wp = tan.(π/2 .* Wp)
    ws = tan.(π/2 .* Ws)
    wa = _toprototype(wp, ws, ftype)

    N = _order_estimate(Rp, Rs, wa)

    wscale = _natfreq_estimate(wa, Rs, N)
    ωn = (2/π).*atan.(_fromprototype(wp, wscale, ftype))
    N, ωn
end

function buttord(Wp::Real, Ws::Real, Rp::Real, Rs::Real, ftype::Union{Type{Lowpass}, Type{Highpass}})
    """
        buttord(Wp, Ws, Rp, Rs, filtertype)

    Butterworth filter order and -3 dB frequency approximation. The lowpass and highpass filters
    are supported, currently working on the bandpass and bandstop prototype transformations. `Wp`
    and `Ws` are the passband and stopband frequencies, whereas Rp and Rs are the passband and
    stopband ripple attenuations in dB.
    """

    if ((ftype == Lowpass) && (Wp > Ws))
        error("Passband frequency must be less than stopband frequency for Lowpass filters.")
    elseif ((ftype == Highpass) && (Ws > Wp))
        error("Passband frequency must be greater than stopband frequency for Highpass filters.")
    end

    # need to pre-warp since we want to use formulae for analog case.
    wp = tan(π/2 * Wp)
    ws = tan(π/2 * Ws)
    wa = _toprototype(wp, ws, ftype)

    # rounding up fractional order. Using differences of logs instead of division. 
    N = _order_estimate(Rp, Rs, wa)
    
    # specifications for the stopband ripple are met precisely.
    wscale = _natfreq_estimate(wa, Rs, Integer(N))
    
    # convert back to the original analog filter and bilinear xform.
    ωn = (2/π)*atan(_fromprototype(wp, wscale, ftype))
    N, ωn
end