
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
    Wa = (Ws.^2 .- Wp[1]*Wp[2]) ./ (Ws .* (Wp[1]-Wp[2]))
    minimum(abs.(Wa))
end

function toprototype!(Wp::AbstractArray{<:Real}, Ws::AbstractArray{<:Real}, Rp::Real, Rs::Real, ftype::Type{Bandstop})
    # NOTE: the optimization function will adjust the corner frequencies in Wp, modifying the original input.
    Δ = eps(typeof(Wp[1]))^(2/3)

    # optimize order on bound [passband low < w < stopband-tolerance].
    C₁(w) = bsfcost(w, true, Wp, Ws, Rp, Rs)
    p1 = minimizer(optimize(C₁, Wp[1], Ws[1]-Δ))
    Wp[1] = p1 # modifying band edge for next optimization run.
    
    # declaring the 2nd cost function here to make sure the overwritten 
    # passband edge vector Wp is used in the new cost function.
    C₂(w) = bsfcost(w, false, Wp, Ws, Rp, Rs)
    p2 = minimizer(optimize(C₂, Ws[2]+Δ, Wp[2]))
    Wp[2] = p2

    Wa = (Ws .* (Wp[1]-Wp[2])) ./ (Ws.^2 .- (Wp[1]*Wp[2]))
    minimum(abs.(Wa))
end

fromprototype(Wp::Real, Wscale::Real, ftype::Type{Lowpass}) = Wp * Wscale
fromprototype(Wp::Real, Wscale::Real, ftype::Type{Highpass}) = Wp / Wscale

function fromprototype(Wp::AbstractArray{<:Real}, Wscale::Real, ftype::Type{Bandstop})
    Wa = zeros(2)
    diff = Wp[2]-Wp[1]
    prod = Wp[2]*Wp[1]
    Wa[1] = (diff + sqrt(diff^2 + 4*(Wscale^2)*prod)) / (2*Wscale)
    Wa[2] = (diff - sqrt(diff^2 + 4*(Wscale^2)*prod)) / (2*Wscale)
    sort!(abs.(Wa))
end

function fromprototype(Wp::AbstractArray{<:Real}, Wscale::Real, ftype::Type{Bandpass})
    Wsc = [-Wscale, Wscale]
    Wa = -Wsc .* (Wp[2]-Wp[1])./2 .+ sqrt.( Wsc.^2/4*(Wp[2]-Wp[1]).^2 .+ Wp[1]*Wp[2])
    sort!(abs.(Wa))
end

order_estimate(Rp::Real, Rs::Real, warp::Real) = (log10(^(10, 0.1*Rs) - 1) - log10(^(10, 0.1*Rp) - 1)) / (2*log10(warp))
natfreq_estimate(warp::Real, Rs::Real, order::Integer) = warp / (^(10, 0.1*Rs) - 1)^(1/(2*order))

function buttord(Wp::AbstractArray{<:Real}, Ws::AbstractArray{<:Real}, Rp::Real, Rs::Real)
    """
        buttord(Wp, Ws, Rp, Rs)

    Butterworth order estimate for bandpass and bandstop filter types. 
    `Wp` and `Ws` are 2-element pass/stopband frequency edges, with filttype specifying
    the `Bandpass` or `Bandstop` filter types.
    """

    length(Wp) == 2 || error("2 passband frequencies were expected.")
    length(Ws) == 2 || error("2 stopband frequencies were expected.")
    
    # make sure the band edges are in increasing order.
    sort!(Wp)
    sort!(Ws)

    # infer filter type based on ordering of edges.
    ftype = (Wp[1] < Ws[1]) ? Bandstop : Bandpass

    # pre-warp both components
    wp = tan.(π/2 .* Wp)
    ws = tan.(π/2 .* Ws)
    if (ftype == Bandstop)
        # optimizer modifies passband frequencies: this is intended behavior.
        wa = toprototype!(wp, ws, Rp, Rs, ftype)
    else
        wa = toprototype(wp, ws, ftype)
    end

    # get the integer order estimate.
    N = ceil(order_estimate(Rp, Rs, wa))

    wscale = natfreq_estimate(wa, Rs, Integer(N))
    ωn = (2/π).*atan.(fromprototype(wp, wscale, ftype))
    N, ωn
end

function buttord(Wp::Real, Ws::Real, Rp::Real, Rs::Real)
    """
        buttord(Wp, Ws, Rp, Rs)

    LPF/HPF Butterworth filter order and -3 dB frequency approximation. `Wp`
    and `Ws` are the passband and stopband frequencies, whereas Rp and Rs 
    are the passband and stopband ripple attenuations in dB. 
    If the passband is greater than stopband, the filter type is inferred to 
    be for estimating the order of a highpass filter.
    """

    # infer which filter type based on the frequency ordering.
    ftype = (Wp < Ws) ? Lowpass : Highpass

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