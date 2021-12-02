
#
# Prototype filter transformations
#

_toprototype(Wp::Real, Ws::Real, ftype::Lowpass) = Ws/Wp
_toprototype(Wp::Real, Ws::Real, ftype::Highpass) = Wp/Ws

function _toprototype(Wp::Vector{Real}, Ws::Vector{Real}, ftype::Bandpass) 
    # bandpass filter must have two corner frequencies we're computing with
    @assert length(Ws) == 2
    @assert length(Wp) == 2
    Wn = (Ws.^2 .- Wp[1]*Wp[2]) ./ (Ws .* (Wp[1]-Wp[2]))
    Wn
end

_fromprototype(Wp::Real, Wscale::Real, ftype::Lowpass) = Wp*Wscale
_fromprototype(Wp::Real, Wscale::Real, ftype::Highpass) = Wp/Wscale

function _fromprototype(Wp::Vector{Real}, Wscale::Real, ftype::Bandpass)
    @assert length(Wp) == 2
    Wsc = [-Wscale, Wscale]
    Wn = -Wsc .* (Wp[2]-Wp[1])./2 .+ sqrt.( Wsc.^2/4*(Wp[2]-Wp[1]).^2 + Wp[1]*Wp[2])
    sort!(Wn)
end

function buttord(Wp::Real, Ws::Real, Rp::Real, Rs::Real, ftype::Union{Lowpass{:z}, Highpass{:z}})
    """
        buttord(Wp, Ws, Rp, Rs, filtertype)

    Butterworth filter order and -3 dB frequency approximation. The lowpass and highpass filters
    are supported, currently working on the bandpass and bandstop prototype transformations. `Wp`
    and `Ws` are the passband and stopband frequencies, whereas Rp and Rs are the passband and
    stopband ripple attenuations in dB.
    """
    # need to pre-warp since we want to use formulae for analog case.
    wp = tan(π/2 * Wp)
    ws = tan(π/2 * Ws)
    wa = _toprototype(ws, wp, ftype)

    # rounding up fractional order. Using differences of logs instead of division. 
    N = ceil((log10(^(10, 0.1*Rs) - 1) - log10(^(10, 0.1*Rp) - 1)) / (2*log10(wa)))
    
    # specifications for the stopband ripple are met precisely.
    wscale = wa / (^(10, 0.1*Rs) - 1)^(1/(2*N))
    
    # convert back to the original analog filter and bilinear xform.
    Wn = (2/π)*atan(_fromprototype(wp, wscale, ftype))
    order, Wn
end