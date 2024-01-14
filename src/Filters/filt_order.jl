#====================================================================

Filter order estimation routines.
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

Design equations based on [1], Chapter 10. Brent's method with Parabolic interpolation
uses code from [3], modifications from Optim.jl's Brent method in [4].

[1] Proakis, J. G., & Manolakis, D. G. (1996). Digital Signal Processing, Fourth Edition.
Prentice Hall, New Jersey.

[2] Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T.,
Cournapeau, D., ... & Van Mulbregt, P. (2020). SciPy 1.0: fundamental algorithms for scientific
computing in Python. Nature methods, 17(3), 261-272.

[3] W. H. Press, et al. Numerical recipes. Third Edition. Cambridge university press, 2007. Chapter 10.3

[4] P. K. Mogensen and A. N. Riseth, Optim: A mathematical optimization package for Julia. Journal of Open
Source Software, 2018. https://github.com/JuliaNLSolvers/Optim.jl/blob/master/src/univariate/solvers/brent.jl

====================================================================#

sort_W((a, b)::Tuple{Real,Real}) = (a > b) ? (b, a) : (a, b)

toprototype(Wp::Real, Ws::Real, ::Type{Lowpass}) = Ws / Wp
toprototype(Wp::Real, Ws::Real, ::Type{Highpass}) = Wp / Ws
function toprototype(Wp::Tuple{Real,Real}, Ws::Tuple{Real,Real}, ::Type{Bandpass})
    # bandpass filter must have two corner frequencies we're computing with
    Wa = muladd.(-Wp[1], Wp[2], Ws .^ 2) ./ (Ws .* (Wp[1] - Wp[2]))
    min(abs.(Wa)...)
end
toprototype(Wp::Tuple{Real,Real}, Ws::Tuple{Real,Real}, Rp::Real, Rs::Real, ::Type{Bandstop}) = butterworth_bsfmin(Wp, Ws, Rp, Rs)
fromprototype(Wp::Real, Wscale::Real, ::Type{Lowpass}) = Wp * Wscale
fromprototype(Wp::Real, Wscale::Real, ::Type{Highpass}) = Wp / Wscale

function fromprototype(Wp::Tuple{Real,Real}, Wscale::Real, ::Type{Bandstop})
    diff = Wp[2] - Wp[1]
    prod = Wp[2] * Wp[1]
    k = sqrt(muladd(4 * (Wscale^2), prod, diff^2))
    den = 2 * Wscale
    Wa = (diff + k) / den, (diff - k) / den
    sort_W(abs.(Wa))
end

function fromprototype(Wp::Tuple{Real,Real}, Wscale::Real, ::Type{Bandpass})
    Wsc = (-Wscale, Wscale)
    diff = Wp[2] - Wp[1]
    prod = Wp[2] * Wp[1]
    Wa = muladd.(.-Wsc, diff ./ 2, sqrt(muladd(Wscale^2 / 4, diff^2, prod)))
    sort_W(abs.(Wa))
end

butterworth_order_estimate(Rp::Real, Rs::Real, warp::Real) = (log(db2pow(Rs) - 1) - log(db2pow(Rp) - 1)) / (2 * log(warp))
butterworth_natfreq_estimate(warp::Real, Rs::Real, order::Integer) = warp / (db2pow(Rs) - 1)^(1 / (2 * order))

function elliptic_order_estimate(Rp::Real, Rs::Real, Wa::Real)
    # Elliptic integer order estmate. Requires Complete Elliptic integral of first kind.
    ϵ = √(db2pow(Rp) - 1) # define selectivity/discrimination parameters.
    k₁ = ϵ / √(db2pow(Rs) - 1)
    k = Wa^-1
    (k^2 < 1) || throw(DomainError(k^2, "Selectivity parameter specifies too narrow of a transition width.")) # check if in-bounds (k₁ << k <~ 1)
    (1 - k₁^2 < 1) || throw(DomainError(1 - k₁^2, "Discrimination parameter specifies too deep of a stopband."))
    K = (ellipk(k^2), ellipk(1 - k^2)) # define the complementary moduli
    K₁ = (ellipk(k₁^2), ellipk(1 - k₁^2))

    # other order approach would be using (k₁'/k₁) / (k′/k).
    (K[1] * K₁[2]) / (K[2] * K₁[1])
end

function chebyshev_order_estimate(Rp::Real, Rs::Real, Wa::Real)
    # Chebyshev1/2 order estimate. Requires inverse hyperbolic cosine.
    ϵs, ϵp = db2pow(Rs) - 1, db2pow(Rp) - 1
    acosh(√(ϵs / ϵp)) / acosh(Wa) # Eq. (10.3.63) in [1]
end

function brent(f, x1::T, x2::T) where {T<:AbstractFloat}
    a, b = x1, x2
    fa, fb = f(a), f(b)
    ϵ, rtol = eps(T), √(eps(T))
    g = 1 / 2 * (3 - √(5)) # golden ratio
    k = zero(T)
    k1 = zero(T)

    # first interval calculation.
    m = muladd(g, (b - a), a)
    m1, m2 = m, m
    fm = f(m)
    fm1, fm2 = fm, fm
    iter = 0

    while (iter < 1_000)
        p, q = zero(T), zero(T)
        xt = rtol * abs(m) + ϵ
        c = (b + a) / 2
        if (abs(m - c) ≤ 2 * xt - (b - a) / 2)
            break
        end
        iter += 1
        if (abs(k1) > xt) # trial parabolic fit
            r = (m - m1) * (fm - fm2)
            q = (m - m2) * (fm - fm1)
            p = muladd((m - m2), q, -(m - m1) * r)
            q = 2 * (q - r)
            if (q > 0)
                p = -p
            else
                q = -q
            end
        end
        if abs(p) < abs(q * k1 / 2) && (p < q * (b - m)) && (p < q * (m - a)) # parabolic step
            k1 = k
            k = p / q
            if ((m + k) < 2 * xt) || ((b - (m + k)) < 2 * xt)
                k = (m < c) ? xt : -xt
            end
        else
            k1 = (m < c) ? b - m : a - m
            k = g * k1
        end

        if (abs(k) > xt)
            xn = m + k
        else
            xn = m + ((k > 0) ? xt : -xt)
        end
        # re-assign variables for next iteration, ("Housekeeping" in [3])
        fn = f(xn)
        if (fn < fm)
            if (xn < m)
                b = m
            else
                a = m
            end
            m2 = m1
            fm2 = fm1
            m1 = m
            fm1 = fm
            m = xn
            fm = fn
        else
            if (xn < m)
                a = xn
            else
                b = xn
            end
            if (fn ≤ fm1) || (m1 == m)
                m2 = m1
                fm2 = fm1
                m1 = xn
                fm1 = fn
            elseif (fn ≤ fm2 || m2 == m || m2 == m1)
                m2 = xn
                fm2 = fn
            end
        end
    end
    m
end

#
# Bandstop cost functions and passband minimization
#
for filt in (:butterworth, :elliptic, :chebyshev)
    @eval begin
        function $(Symbol(filt, :_bsfcost))(Wx::Real, uselowband::Bool, Wp::Tuple{Real,Real}, Ws::Tuple{Real,Real}, Rp::Real, Rs::Real)
            # override one of the passband edges with the test frequency.
            Wpc = uselowband ? (Wx, Wp[2]) : (Wp[1], Wx)

            # get the new warp frequency.
            warp = min(abs.((Ws .* (Wpc[1] - Wpc[2])) ./ muladd.(-Wpc[1], Wpc[2], Ws .^ 2))...)

            # use the new frequency to determine the filter order.
            $(Symbol(filt, :_order_estimate))(Rp, Rs, warp)
        end

        function $(Symbol(filt, :_bsfmin))(Wp::Tuple{Real,Real}, Ws::Tuple{Real,Real}, Rp::Real, Rs::Real)
            # NOTE: the optimization function will adjust the corner frequencies in Wp, returning a new passband tuple.
            Δ = eps(typeof(Wp[1]))^(2 / 3)
            C₁(w) = $(Symbol(filt, :_bsfcost))(w, true, Wp, Ws, Rp, Rs)
            p1 = brent(C₁, Wp[1], Ws[1] - Δ)

            C₂(w) = $(Symbol(filt, :_bsfcost))(w, false, (p1, Wp[2]), Ws, Rp, Rs)
            p2 = brent(C₂, Ws[2] + Δ, Wp[2])
            Wadj = (p1, p2)

            Wa = (Ws .* (p1 - p2)) ./ muladd.(-p1, p2, Ws .^ 2)
            min(abs.(Wa)...), Wadj
        end
    end
end

"""
    (N, ωn) = buttord(Wp::Tuple{Real,Real}, Ws::Tuple{Real,Real}, Rp::Real, Rs::Real; domain=:z)

Butterworth order estimate for bandpass and bandstop filter types.  `Wp` and `Ws` are 2-element pass
and stopband frequency edges, with no more than `Rp` dB passband ripple and at least `Rs` dB stopband
attenuation. Based on the ordering of passband and bandstop edges,  the Bandstop or Bandpass filter
type is inferred. `N` is an integer indicating the lowest estimated filter order, with `ωn` specifying
the cutoff or "-3 dB" frequencies. If a domain of `:s` is specified, the passband and stopband frequencies
are interpreted as radians/second, giving an order and natural frequencies for an analog filter. The default
domain is `:z`, interpretting the input frequencies as normalized from 0 to 1, where 1 corresponds to π radians/sample.
"""
function buttord(Wp::Tuple{Real,Real}, Ws::Tuple{Real,Real}, Rp::Real, Rs::Real; domain::Symbol=:z)

    # make sure the band edges are in increasing order.
    Wps = sort_W(Wp)
    Wss = sort_W(Ws)

    # infer filter type based on ordering of edges.
    (Wps[1] < Wss[1]) != (Wps[2] > Wss[2]) && throw(ArgumentError("Pass and stopband edges must be ordered for Bandpass/Bandstop filters."))
    ftype = (Wps[1] < Wss[1]) ? Bandstop : Bandpass

    # pre-warp both components, (if Z-domain specified)
    if (domain == :z)
        Ωp = tan.(π / 2 .* Wps)
        Ωs = tan.(π / 2 .* Wss)
    else
        Ωp = Wps
        Ωs = Wss
    end

    if (ftype == Bandstop)
        # optimizer modifies passband frequencies.
        (wa, wpadj) = toprototype(Ωp, Ωs, Rp, Rs, ftype)
    else
        # no optimization in BPF case, use original warped passband edges.
        wa = toprototype(Ωp, Ωs, ftype)
        wpadj = Ωp
    end

    # get the integer order estimate.
    N = ceil(Int, butterworth_order_estimate(Rp, Rs, wa))

    wscale = butterworth_natfreq_estimate(wa, Rs, N)
    if (domain == :z)
        ωn = (2 / π) .* atan.(fromprototype(wpadj, wscale, ftype))
    else
        ωn = fromprototype(wpadj, wscale, ftype)
    end
    N, ωn
end

"""
    (N, ωn) = buttord(Wp::Real, Ws::Real, Rp::Real, Rs::Real; domain=:z)

LPF/HPF Butterworth filter order and -3 dB frequency approximation. `Wp` and `Ws` are
the passband and stopband frequencies, whereas Rp and Rs  are the passband and stopband
ripple attenuations in dB. If the passband is greater than stopband, the filter type
is inferred to be for estimating the order of a highpass filter. `N` specifies the lowest
possible integer filter order, whereas `ωn` is the cutoff or "-3 dB" frequency. If a domain
of `:s` is specified, the passband and stopband edges are interpreted as radians/second,
giving an order and natural frequency result for an analog filter. The default domain
is `:z`, interpretting the input frequencies as normalized from 0 to 1, where 1 corresponds
to π radians/sample.
"""
function buttord(Wp::Real, Ws::Real, Rp::Real, Rs::Real; domain::Symbol=:z)
    # infer which filter type based on the frequency ordering.
    ftype = (Wp < Ws) ? Lowpass : Highpass

    if (domain == :z)
        # need to pre-warp since we want to use formulae for analog case.
        Ωp = tan(π / 2 * Wp)
        Ωs = tan(π / 2 * Ws)
    else
        # already analog
        Ωp = Wp
        Ωs = Ws
    end
    wa = toprototype(Ωp, Ωs, ftype)

    # rounding up fractional order. Using differences of logs instead of division.
    N = ceil(Int, butterworth_order_estimate(Rp, Rs, wa))

    # specifications for the stopband ripple are met precisely.
    wscale = butterworth_natfreq_estimate(wa, Rs, N)

    # convert back to the original analog filter
    if (domain == :z)
        # bilinear xform
        ωn = (2 / π) * atan(fromprototype(Ωp, wscale, ftype))
    else
        # s-domain, no atan call.
        ωn = fromprototype(Ωp, wscale, ftype)
    end
    N, ωn
end

#
# Elliptic/Chebyshev1 Estimation
#
for (fcn, est, filt) in ((:ellipord, :elliptic, "Elliptic (Cauer)"),
    (:cheb1ord, :chebyshev, "Chebyshev Type I"))
    @eval begin
        """
            (N, ωn) = $($fcn)(Wp::Real, Ws::Real, Rp::Real, Rs::Real; domain::Symbol=:z)

        Integer and natural frequency order estimation for $($filt) Filters. `Wp`
        and `Ws` indicate the passband and stopband frequency edges, and `Rp` and `Rs` indicate
        the maximum loss in the passband and the minimum attenuation in the stopband, (in dB.)
        Based on the ordering of passband and stopband edges, the Lowpass or Highpass filter type is inferred.
        `N` indicates the smallest integer filter order that achieves the desired specifications,
        and `ωn` contains the natural frequency of the filter, (in this case, simply the passband edge.)
        If a domain of `:s` is specified, the passband and stopband edges are interpreted
        as radians/second, giving an order and natural frequency result for an analog filter.
        The default domain is `:z`, interpretting the input frequencies as normalized from 0 to 1,
        where 1 corresponds to π radians/sample.
        """
        function $fcn(Wp::Real, Ws::Real, Rp::Real, Rs::Real; domain::Symbol=:z)
            ftype = (Wp < Ws) ? Lowpass : Highpass
            if (domain == :z)
                Ωp = tan(π / 2 * Wp)
                Ωs = tan(π / 2 * Ws)
            else
                Ωp = Wp
                Ωs = Ws
            end
            # Lowpass/Highpass prototype xform is same as Butterworth.
            wa = toprototype(Ωp, Ωs, ftype)
            N = ceil(Int, $(Symbol(est, :_order_estimate))(Rp, Rs, wa))
            if (domain == :z)
                ωn = (2 / π) * atan(Ωp)
            else
                ωn = Ωp
            end
            N, ωn
        end

        """
            (N, ωn) = $($fcn)(Wp::Tuple{Real, Real}, Ws::Tuple{Real, Real}, Rp::Real, Rs::Real; domain::Symbol=:z)

        Integer and natural frequency order estimation for $($filt) Filters. `Wp` and `Ws` are 2-element
        frequency edges for Bandpass/Bandstop cases, with `Rp` and `Rs` representing the ripple maximum loss
        in the passband and minimum ripple attenuation in the stopband in dB. Based on the ordering of passband
        and bandstop edges, the Bandstop or Bandpass filter type is inferred. `N` is an integer indicating the
        lowest estimated filter order, with `ωn` specifying the cutoff or "-3 dB" frequencies. If a domain of
        `:s` is specified, the passband and stopband frequencies are interpreted as radians/second, giving an
        order and natural frequencies for an analog filter. The default domain is `:z`, interpreting the input
        frequencies as normalized from 0 to 1, where 1 corresponds to π radians/sample.
        """
        function $fcn(Wp::Tuple{Real,Real}, Ws::Tuple{Real,Real}, Rp::Real, Rs::Real; domain::Symbol=:z)
            Wps = sort_W(Wp)
            Wss = sort_W(Ws)
            (Wps[1] < Wss[1]) != (Wps[2] > Wss[2]) && throw(ArgumentError("Pass and stopband edges must be ordered for Bandpass/Bandstop filters."))
            ftype = (Wps[1] < Wss[1]) ? Bandstop : Bandpass
            # pre-warp to analog if z-domain.
            (Ωp, Ωs) = (domain == :z) ? (tan.(π / 2 .* Wps), tan.(π / 2 .* Wss)) : (Wps, Wss)
            if (ftype == Bandpass)
                Wa = muladd.(-Ωp[1], Ωp[2], Ωs .^ 2) ./ (Ωs .* (Ωp[1] - Ωp[2]))
                Ωpadj = Ωp
            else
                (Wa, Ωpadj) = $(Symbol(est, :_bsfmin))(Ωp, Ωs, Rp, Rs) # check scipy.
            end
            N = ceil(Int, $(Symbol(est, :_order_estimate))(Rp, Rs, min(abs.(Wa)...)))
            ωn = (domain == :z) ? Wps : Ωpadj
            N, ωn
        end
    end
end


"""
    (N, ωn) = cheb2ord(Wp::Real, Ws::Real, Rp::Real, Rs::Real; domain::Symbol=:z)

Integer and natural frequency order estimation for Chebyshev Type II (inverse) Filters. `Wp` and `Ws`
are the frequency edges for Bandpass/Bandstop cases, with `Rp` and `Rs` representing the ripple maximum
loss in the passband and minimum ripple attenuation in the stopband in dB. Based on the ordering of
the passband and stopband edges, the Lowpass or Highpass filter type is inferred. `N` is the integer
indicating the lowest filter order, with `ωn` specifying the "-3 dB" cutoff frequency. If the domain
is specified as `:s`, the passband and stopband frequencies are interpreted as radians/second, giving
the order and natural frequencies for an analog filter. The default domain is `:z`, interpretting the
input frequencies as normalized from 0 to 1, where 1 corresponds to π radians/sample.
"""
function cheb2ord(Wp::Real, Ws::Real, Rp::Real, Rs::Real; domain::Symbol=:z)
    ftype = (Wp < Ws) ? Lowpass : Highpass
    (Ωp, Ωs) = (domain == :z) ? (tan(π / 2 * Wp), tan(π / 2 * Ws)) : (Wp, Ws)
    wa = toprototype(Ωp, Ωs, ftype)
    N = ceil(Int, chebyshev_order_estimate(Rp, Rs, wa))

    # new frequency for stopband spec.
    wnew = 1 / cosh(1 / N * acosh(√(db2pow(Rs) - 1) / √(db2pow(Rp) - 1)))
    wa = (ftype == Lowpass) ? Ωp / wnew : Ωp * wnew
    ωn = (domain == :z) ? (2 / π) * atan(wa) : wa
    N, ωn
end


"""
    (N, ωn) = cheb2ord(Wp::Tuple{Real, Real}, Ws::Tuple{Real, Real}, Rp::Real, Rs::Real; domain::Symbol=:z)

Integer and natural frequency order estimation for Chebyshev Type II (inverse) Filters. `Wp` and `Ws` are
2-element frequency edges for Bandpass/Bandstop cases, with `Rp` and `Rs` representing  the ripple maximum
loss in the passband and minimum ripple attenuation in the stopband in dB. Based on the ordering of passband
and bandstop edges, the Bandstop or Bandpass filter type is inferred. `N` is an integer indicating the
lowest estimated filter order, with `ωn` specifying the cutoff or "-3 dB" frequencies. If a domain of
`:s` is specified, the passband and stopband frequencies are interpreted as radians/second, giving an
order and natural frequencies for an analog filter. The default domain is `:z`, interpretting the input
frequencies as normalized from 0 to 1, where 1 corresponds to π radians/sample.
"""
function cheb2ord(Wp::Tuple{Real,Real}, Ws::Tuple{Real,Real}, Rp::Real, Rs::Real; domain::Symbol=:z)
    Wps = sort_W(Wp)
    Wss = sort_W(Ws)
    (Wps[1] < Wss[1]) != (Wps[2] > Wss[2]) && throw(ArgumentError("Pass and stopband edges must be ordered for Bandpass/Bandstop filters."))
    ftype = (Wps[1] < Wss[1]) ? Bandstop : Bandpass
    (Ωp, Ωs) = (domain == :z) ? (tan.(π / 2 .* Wps), tan.(π / 2 .* Wss)) : (Wps, Wss)
    if (ftype == Bandpass)
        prod = Ωp[1] * Ωp[2]
        diff = Ωp[1] - Ωp[2]
        Wa = muladd.(Ωs, Ωs, -prod) ./ (Ωs .* diff)
    else
        (Wa, Ωpadj) = chebyshev_bsfmin(Ωp, Ωs, Rp, Rs)
        prod = Ωpadj[1] * Ωpadj[2]
        diff = Ωpadj[1] - Ωpadj[2]
    end
    N = ceil(Int, chebyshev_order_estimate(Rp, Rs, min(abs.(Wa)...)))

    # new frequency for stopband spec.
    wnew = 1 / cosh(1 / N * acosh(√(db2pow(Rs) - 1) / √(db2pow(Rp) - 1)))

    # lowpass prototype to analog filter re-map.
    if (ftype == Bandpass)
        Wna1 = diff / (2 * wnew) + √(diff^2 / (4 * wnew^2) + prod)
    else
        Wna1 = (diff * wnew) / 2 + √(muladd(diff^2, wnew^2 / 4, prod))
    end
    Wna2 = prod / Wna1
    ωn = (domain == :z) ? ((2 / π) * atan(Wna1), (2 / π) * atan(Wna2)) : (Wna1, Wna2)
    N, ωn
end

"""
    N = remezord(Wp, Ws, Rp, Rs)

Order estimation for lowpass digital filter cases based
on the equations and coefficients in [^Rabiner]. The original
equation returned the minimum filterlength, whereas this implementation
returns the order (N=L-1). `Wp` and `Ws` are the normalized passband and
stopband frequencies, with `Rp` indicating the passband ripple and `Rs` is the
stopband attenuation (linear.)

NOTE: The value of N is an initial approximation. If under-estimated, the order
should be increased until the design specifications are met.

[^Rabiner]: Herrmann, O., Lawrence R. Rabiner, and D. S. K. Chan.
    "Practical design rules for optimum finite impulse response lowpass digital filters."
    Bell System Technical Journal 52.6 (1973): 769-799.
"""
function remezord(Wp::Real, Ws::Real, Rp::Real, Rs::Real)
    (0 > Wp > 0.5) || (0 > Ws > 0.5) && throw(ArgumentError("Pass and stopband edges must be greater than DC and less than Nyquist."))
    L1, L2 = log10(Rp), log10(Rs)
    df = abs(Ws - Wp) # works in HPF case if passband/stopband edges are flipped
    A = muladd(5.309e-3, L1^2, muladd(7.114e-2, L1, -0.4761))
    B = muladd(2.66e-3, L1^2, muladd(0.5941, L1, 0.4278))
    Kf = muladd(0.51244, (L1 - L2), 11.01217)
    D = muladd(A, L2, -B)
    return ceil(Int, (muladd(-Kf, df^2, D) / df))
end