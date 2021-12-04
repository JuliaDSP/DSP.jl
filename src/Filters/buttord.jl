
#
# Butterworth prototype filter transformations
#


function bsfcost(Wx::Real, uselowband::Bool, Wp::AbstractArray{<:Real}, Ws::AbstractArray{<:Real}, Rp::Real, Rs::Real)
    """
        bsfcost(Wx, uselowband, Wp, Ws, Rs, Rp)

    Bandstop filter order estimation. The primary variables are `Wx` and `uselowband`,
    which indicate the test passband edge and if to use the lower edge. If false, the
    test frequency is used in high-band, ([low, high] ordering in Wp.) This function
    returns a non-integer BSF order estimate.
    """

    # override one of the passband edges with the test frequency.
    Wpc = uselowband ? [Wx, Wp[2]] : [Wp[1], Wx]
    
    # get the new warp frequency.
    warp = minimum(abs.((Ws .* (Wpc[1]-Wpc[2])) ./ (Ws.^2 .- (Wpc[1]*Wpc[2]))))

    # use the new frequency to determine the new filter order.
    order_estimate(Rp, Rs, warp)
end

toprototype(Wp::Real, Ws::Real, ftype::Type{Lowpass}) = Ws / Wp
toprototype(Wp::Real, Ws::Real, ftype::Type{Highpass}) = Wp / Ws

function toprototype(Wp::AbstractArray{<:Real}, Ws::AbstractArray{<:Real}, ftype::Type{Bandpass}) 
    # bandpass filter must have two corner frequencies we're computing with
    length(Wp) == 2 || error("2 passband frequencies were expected.")
    length(Ws) == 2 || error("2 stopband frequencies were expected.")
    Wa = (Ws.^2 .- Wp[1]*Wp[2]) ./ (Ws .* (Wp[1]-Wp[2]))
    minimum(abs.(Wa))
end

function toprototype!(Wp::AbstractArray{<:Real}, Ws::AbstractArray{<:Real}, Rp::Real, Rs::Real, ftype::Type{Bandstop})
    # NOTE: the optimization function will adjust the corner frequencies in Wp, modifying the original input.
    length(Wp) == 2 || error("2 passband frequencies were expected.")
    length(Ws) == 2 || error("2 stopband frequencies were expected.")

    C₁(w, px) = bsfcost(w, true, Wp, Ws, Rp, Rs)
    C₂(w, px) = bsfcost(w, false, Wp, Ws, Rp, Rs)
    tol = eps(typeof(Wp[1]))^(2/3)

    # optimize order on bound [passband low < w < stopband-tolerance].
    p1 = minimizer(optimize(C₁, Wp[1], Ws[1]-tol))
    Wp[1] = p1 # modifying band edge for next optimization run.
    
    p2 = minimizer(optimize(C₂, Ws[2]+tol, Wp[2]))
    Wp[2] = p2

    Wa = (Ws .* (Wp[1]-Wp[2])) ./ (Ws.^2 .- (Wp[1]-Wp[2]))
    minimum(abs.(Wa))
end

fromprototype(Wp::Real, Wscale::Real, ftype::Type{Lowpass}) = Wp * Wscale
fromprototype(Wp::Real, Wscale::Real, ftype::Type{Highpass}) = Wp / Wscale

function fromprototype(Wp::AbstractArray{<:Real}, Wscale::Real, ftype::Type{Bandstop})
    length(Wp) == 2 || error("2 passband frequencies were  expected.")
    Wa = zeros(2)
    diff = Wp[2]-Wp[1]
    prod = Wp[2]*Wp[1]
    Wa[1] = (diff + sqrt(diff^2 + 4*(Wscale^2)*prod)) / (2*Wscale)
    Wa[2] = (diff - sqrt(diff^2 + 4*(Wscale^2)*prod)) / (2*Wscale)
    sort!(abs.(Wa))
    Wa
end

function fromprototype(Wp::AbstractArray{<:Real}, Wscale::Real, ftype::Type{Bandpass})
    length(Wp) == 2 || error("2 passband frequencies were expected.")
    Wsc = [-Wscale, Wscale]
    Wa = -Wsc .* (Wp[2]-Wp[1])./2 .+ sqrt.( Wsc.^2/4*(Wp[2]-Wp[1]).^2 + Wp[1]*Wp[2])
    sort!(abs.(Wa))
    Wa
end

order_estimate(Rp::Real, Rs::Real, warp::Real) = (log10(^(10, 0.1*Rs) - 1) - log10(^(10, 0.1*Rp) - 1)) / (2*log10(warp))
natfreq_estimate(warp::Real, Rs::Real, order::Integer) = warp / (^(10, 0.1*Rs) - 1)^(1/(2*order))

function buttord(Wp::AbstractArray{<:Real}, Ws::AbstractArray{<:Real}, Rp::Real, Rs::Real, ftype::Union{Type{Bandpass}, Type{Bandstop}})
    """
        buttord(Wp, Ws, Rp, Rs, filttype)

    Butterworth order estimate for bandpass and bandstop filter types. 
    `Wp` and `Ws` are 2-element pass/stopband frequency edges, with filttype specifying
    the `Bandpass` or `Bandstop` filter types.
    """
    
    # pre-warp both components
    wp = tan.(π/2 .* Wp)
    ws = tan.(π/2 .* Ws)
    if (type == Bandstop)
        # optimizer modifies passband frequencies. this is intended behavior.
        wa = toprototype!(wp, ws, Rp, Rs, ftype)
    else
        wa = toprototype(wp, ws, ftype)
    end

    # get the integer order estimate.
    N = ceil(order_estimate(Rp, Rs, wa))

    wscale = natfreq_estimate(wa, Rs, N)
    ωn = (2/π).*atan.(fromprototype(wp, wscale, ftype))
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
    wa = toprototype(wp, ws, ftype)

    # rounding up fractional order. Using differences of logs instead of division. 
    N = ceil(order_estimate(Rp, Rs, wa))
    
    # specifications for the stopband ripple are met precisely.
    wscale = natfreq_estimate(wa, Rs, Integer(N))
    
    # convert back to the original analog filter and bilinear xform.
    ωn = (2/π)*atan(fromprototype(wp, wscale, ftype))
    N, ωn
end