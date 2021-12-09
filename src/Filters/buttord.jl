
#====================================================================

Butterworth filter order estimation.
Jordan R. Smith, 2021.

Filter transformations translated from Scipy to Julia.

SCIPY license:
Copyright (c) 2001, 2002 Enthought, Inc.
All rights reserved.

Copyright (c) 2003-2017 SciPy Developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of Enthought nor the names of the SciPy Developers
     may be used to endorse or promote products derived from this software
     without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.

--- end of scipy license

Design equations based on [1], Chapter 10.

[1] Proakis, J. G., & Manolakis, D. G. (1996). Digital Signal Processing, Fourth Edition. 
Prentice Hall, New Jersey.
[2] Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., 
Cournapeau, D., ... & Van Mulbregt, P. (2020). SciPy 1.0: fundamental algorithms for scientific 
computing in Python. Nature methods, 17(3), 261-272.

====================================================================#

"""
    N = bsfcost(Wx, uselowband, Wp, Ws, Rs, Rp)

Bandstop filter order estimation. The primary variables are `Wx` and `uselowband`,
which indicate the test passband edge and if to use the lower edge. If false, the
test frequency is used in high-band, ([low, high] ordering in Wp.) This function
returns a non-integer BSF order estimate.
"""
function bsfcost(Wx::Real, uselowband::Bool, Wp::Tuple{Real,Real}, Ws::Tuple{Real,Real}, Rp::Real, Rs::Real)

    # override one of the passband edges with the test frequency.
    Wpc = uselowband ? tuple(Wx, Wp[2]) : tuple(Wp[1], Wx)
    
    # get the new warp frequency.
    warp = minimum(abs.((Ws .* (Wpc[1]-Wpc[2])) ./ (Ws.^2 .- (Wpc[1]*Wpc[2]))))

    # use the new frequency to determine the new filter order.
    order_estimate(Rp, Rs, warp)
end

toprototype(Wp::Real, Ws::Real, ftype::Type{Lowpass}) = Ws / Wp
toprototype(Wp::Real, Ws::Real, ftype::Type{Highpass}) = Wp / Ws

function toprototype(Wp::Tuple{Real,Real}, Ws::Tuple{Real,Real}, ftype::Type{Bandpass}) 
    # bandpass filter must have two corner frequencies we're computing with
    Wa = (Ws.^2 .- Wp[1]*Wp[2]) ./ (Ws .* (Wp[1]-Wp[2]))
    minimum(abs.(Wa))
end

function toprototype(Wp::Tuple{Real,Real}, Ws::Tuple{Real,Real}, Rp::Real, Rs::Real, ftype::Type{Bandstop})
    # NOTE: the optimization function will adjust the corner frequencies in Wp, returning a new passband tuple.
    Δ = eps(typeof(Wp[1]))^(2/3)

    # optimize order on bound [passband low < w < stopband-tolerance].
    C₁(w) = bsfcost(w, true, Wp, Ws, Rp, Rs)
    p1 = minimizer(optimize(C₁, Wp[1], Ws[1]-Δ))
    
    # declaring the 2nd cost function here to use the new passband tuple.
    C₂(w) = bsfcost(w, false, tuple(p1, Wp[2]), Ws, Rp, Rs)
    p2 = minimizer(optimize(C₂, Ws[2]+Δ, Wp[2]))
    Wadj = tuple(p1, p2)

    Wa = (Ws .* (Wadj[1]-Wadj[2])) ./ (Ws.^2 .- (Wadj[1]*Wadj[2]))
    minimum(abs.(Wa)), Wadj
end

fromprototype(Wp::Real, Wscale::Real, ftype::Type{Lowpass}) = Wp * Wscale
fromprototype(Wp::Real, Wscale::Real, ftype::Type{Highpass}) = Wp / Wscale

function fromprototype(Wp::Tuple{Real,Real}, Wscale::Real, ftype::Type{Bandstop})
    Wa = zeros(2)
    diff = Wp[2]-Wp[1]
    prod = Wp[2]*Wp[1]
    Wa[1] = (diff + sqrt(diff^2 + 4*(Wscale^2)*prod)) / (2*Wscale)
    Wa[2] = (diff - sqrt(diff^2 + 4*(Wscale^2)*prod)) / (2*Wscale)
    sort!(abs.(Wa))
end

function fromprototype(Wp::Tuple{Real,Real}, Wscale::Real, ftype::Type{Bandpass})
    Wsc = [-Wscale, Wscale]
    Wa = -Wsc .* (Wp[2]-Wp[1])./2 .+ sqrt.( Wsc.^2/4*(Wp[2]-Wp[1]).^2 .+ Wp[1]*Wp[2])
    sort!(abs.(Wa))
end

order_estimate(Rp::Real, Rs::Real, warp::Real) = (log(db2pow(Rs) - 1) - log(db2pow(Rp) - 1)) / (2*log(warp))
natfreq_estimate(warp::Real, Rs::Real, order::Integer) = warp / (db2pow(Rs) - 1)^(1/(2*order))


"""
    (N, ωn) = buttord(Wp::Tuple{Real,Real}, Ws::Tuple{Real,Real}, Rp::Real, Rs::Real; domain=:z)

Butterworth order estimate for bandpass and bandstop filter types. 
`Wp` and `Ws` are 2-element pass and stopband frequency edges, with 
no more than `Rp` dB passband ripple and at least `Rs` dB stopband 
attenuation. Based on the ordering of passband and bandstop edges, 
the Bandstop or Bandpass filter type is inferred. `N` is an integer 
indicating the lowest estimated filter order, with `ωn` specifying
the cutoff or "-3 dB" frequencies. If a domain of `:s` is specified,
the passband and stopband frequencies are interpretted as radians/second,
giving an order and natural frequencies for an analog filter. The default
domain is `:z`, interpretting the input frequencies as normalized from 0 to 1,
where 1 corresponds to π radians/sample.
"""
function buttord(Wp::Tuple{Real,Real}, Ws::Tuple{Real,Real}, Rp::Real, Rs::Real; domain::Symbol=:z)

    # make sure the band edges are in increasing order.
    Wps = (Wp[1] > Wp[2]) ? tuple(Wp[2], Wp[1]) : tuple(Wp[1], Wp[2])
    Wss = (Ws[1] > Ws[2]) ? tuple(Ws[2], Ws[1]) : tuple(Ws[1], Ws[2])

    # infer filter type based on ordering of edges.
    ftype = (Wps[1] < Wss[1]) ? Bandstop : Bandpass

    # pre-warp both components, (if Z-domain specified)
    if (domain == :z)
        wp = tan.(π/2 .* Wps)
        ws = tan.(π/2 .* Wss)
    else
        wp = Wps
        ws = Wss
    end

    if (ftype == Bandstop)
        # optimizer modifies passband frequencies.
        (wa, wpadj) = toprototype(wp, ws, Rp, Rs, ftype)
    else
        # no optimization in BPF case, use original warped passband edges.
        wa = toprototype(wp, ws, ftype)
        wpadj = wp
    end

    # get the integer order estimate.
    N = ceil(Int, order_estimate(Rp, Rs, wa))

    wscale = natfreq_estimate(wa, Rs, N)
    if (domain == :z)
        ωn = (2/π).*atan.(fromprototype(wpadj, wscale, ftype))
    else
        ωn = fromprototype(wpadj, wscale, ftype)
    end
    N, ωn
end

"""
    (N, ωn) = buttord(Wp::Real, Ws::Real, Rp::Real, Rs::Real; domain=:z)

LPF/HPF Butterworth filter order and -3 dB frequency approximation. `Wp`
and `Ws` are the passband and stopband frequencies, whereas Rp and Rs 
are the passband and stopband ripple attenuations in dB. 
If the passband is greater than stopband, the filter type is inferred to 
be for estimating the order of a highpass filter. `N` specifies the lowest
possible integer filter order, whereas `ωn` is the cutoff or "-3 dB" frequency.
If a domain of `:s` is specified, the passband and stopband edges are interpretted
as radians/second, giving an order and natural frequency result for an analog filter.
The default domain is `:z`, interpretting the input frequencies as normalized from 0 to 1,
where 1 corresponds to π radians/sample.
"""
function buttord(Wp::Real, Ws::Real, Rp::Real, Rs::Real; domain::Symbol=:z)
    # infer which filter type based on the frequency ordering.
    ftype = (Wp < Ws) ? Lowpass : Highpass

    if (domain == :z)
        # need to pre-warp since we want to use formulae for analog case.
        wp = tan(π/2 * Wp)
        ws = tan(π/2 * Ws)
    else
        # already analog
        wp = Wp
        ws = Ws
    end
    wa = toprototype(wp, ws, ftype)

    # rounding up fractional order. Using differences of logs instead of division. 
    N = ceil(Int, order_estimate(Rp, Rs, wa))
    
    # specifications for the stopband ripple are met precisely.
    wscale = natfreq_estimate(wa, Rs, N)
    
    # convert back to the original analog filter
    if (domain == :z)
        # bilinear xform to digital
        ωn = (2/π)*atan(fromprototype(wp, wscale, ftype))
    else
        # s-domain, no atan call.
        ωn = fromprototype(wp, wscale, ftype)
    end
    N, ωn
end