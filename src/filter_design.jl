# Filter design functions for Julia
# Created and (c) Simon Kornblith <simon@simonster.com>
#
# I know very little about filter design. While these functions seem to work,
# their numerical stability may be poor. You should really make sure that
# these functions do what you want before relying upon them.
#
# Insipred by scipy.signal's filter_design.py
#
# Copyright (c) 2001, 2002 Enthought, Inc.
# All rights reserved.
#
# Copyright (c) 2003-2012 SciPy Developers.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#   a. Redistributions of source code must retain the above copyright notice,
#      this list of conditions and the following disclaimer.
#   b. Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#   c. Neither the name of Enthought nor the names of the SciPy Developers
#      may be used to endorse or promote products derived from this software
#      without specific prior written permission.
#
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

module FilterDesign
using Polynomial

export ZPKFilter, TFFilter, BiquadFilter, SOSFilter, Butterworth, Lowpass, Highpass, Bandpass,
       Bandstop, analogfilter, digitalfilter, Filter
import Base: convert, filt

#
# Utility functions
#

# Pad a vector with f on the left to a specified length
function lfill{T}(v::Vector{T}, n::Integer, f::Number=zero(T))
    d = n-length(v)
    if d > 0
        nv = similar(v, n)
        nv[1:d] = f
        nv[d+1:end] = v
        nv
    else
        v
    end
end

# Make a new array of zeros of the same type as another vector
similarzeros{T}(v::Array{T}, args...) = zeros(T, args...)

# Get coefficients of a polynomial
coeffs{T}(p::Poly{T}) = p.a[1+p.nzfirst:end]

#
# Filter types
#

abstract Filter

# Filter in zero-pole-gain form
immutable ZPKFilter{Z,P,K} <: Filter
    z::Vector{Z}
    p::Vector{P}
    k::K
end

# Filter in transfer function (numerator and denominator) form
immutable TFFilter{T} <: Filter
    b::Poly{T}
    a::Poly{T}

    function TFFilter(b::Poly{T}, a::Poly{T})
        new(b/a[1], a/a[1])
    end
end
TFFilter{T}(b::Poly{T}, a::Poly{T}) = TFFilter{T}(b, a)
TFFilter{T}(b::Vector{T}, a::Vector{T}) = TFFilter{T}(Poly(b), Poly(a))
function TFFilter{T,S}(b::Vector{T}, a::Vector{S})
    V = promote_type(T, S)
    TFFilter(convert(Vector{V}, b), convert(Vector{V}, a))
end

function convert(::Type{TFFilter}, f::ZPKFilter)
    b = f.k*poly(f.z)
    a = poly(f.p)
    TFFilter(Poly(real(b.a)), Poly(real(a.a)))
end

function convert{T}(::Type{ZPKFilter}, f::TFFilter{T})
    k = real(f.b[1])
    b = f.b / k
    z = convert(Vector{Complex{T}}, roots(b))
    p = convert(Vector{Complex{T}}, roots(f.a))
    ZPKFilter(z, p, k)
end

#
# Biquad filter in transfer function form
# A separate immutable to improve efficiency of filtering using SOSFilters
#

immutable BiquadFilter{T} <: Filter
    b0::T
    b1::T
    b2::T
    a1::T
    a2::T
end
BiquadFilter{T}(b0::T, b1::T, b2::T, a0::T, a1::T, a2::T, g::Real=1) =
    BiquadFilter(g*b0/a0, g*b1/a0, g*b2/a0, a1/a0, a2/a0)

convert(::Type{ZPKFilter}, f::BiquadFilter) = convert(ZPKFilter, convert(TFFilter, f))

function convert{T}(::Type{TFFilter}, f::BiquadFilter{T})
    if f.b2 == zero(T) && f.a2 == zero(T)
        if f.b1 == zero(T) && f.a1 == zero(T)
            b = [f.b0]
            a = [one(T)]
        else
            b = [f.b0, f.b1]
            a = [one(T), f.a1]
        end
    else
        b = [f.b0, f.b1, f.b2]
        a = [one(T), f.a1, f.a2]
    end

    TFFilter(b, a)
end

convert(::Type{BiquadFilter}, f::ZPKFilter) = convert(BiquadFilter, convert(TFFilter, f))

function convert{T}(::Type{BiquadFilter}, f::TFFilter{T})
    xs = max(length(f.b), length(f.a))
    b = lfill(coeffs(f.b), xs)
    a = lfill(coeffs(f.a), xs)

    if xs == 3
        BiquadFilter(b[1], b[2], b[3], a[2], a[3])
    elseif xs == 2
        BiquadFilter(b[1], b[2], zero(T), a[2], zero(T))
    elseif xs == 1
        BiquadFilter(b[1], zero(T), zero(T), zero(T), zero(T))
    elseif xs == 0
        error("cannot convert an empty TFFilter to BiquadFilter")
    else
        error("cannot convert a filter of length > 3 to BiquadFilter")
    end
end

#
# Filtering as second-order sections
#

immutable SOSFilter{T,G} <: Filter
    biquads::Vector{BiquadFilter{T}}
    g::G
end

realtype(x::DataType) = x
realtype{T}(::Type{Complex{T}}) = T
complextype(T::DataType) = Complex{T}
complextype{T}(::Type{Complex{T}}) = Complex{T}

function convert{T}(::Type{ZPKFilter}, f::SOSFilter{T})
    t = complextype(T)
    z = t[]
    p = t[]
    k = f.g
    for biquad in f.biquads
        biquadzpk = convert(ZPKFilter, biquad)
        append!(z, biquadzpk.z)
        append!(p, biquadzpk.p)
        k *= biquadzpk.k
    end
    ZPKFilter(z, p, k)
end

convert(to::Union(Type{TFFilter}, Type{BiquadFilter}), f::SOSFilter) =
    convert(to, convert(ZPKFilter, f))

# Lexicographic less than function
# This is only necessary in Julia 0.2
function lexlt(a::Complex, b::Complex)
    if real(a) == real(b)
        isless(imag(a), imag(b))
    else
        isless(real(a), real(b))
    end
end
lexlt(a::Real, b::Real) = isless(a, b)

# Convert a filter to second-order sections
# The returned sections are in ZPK form
function convert{Z,P}(::Type{SOSFilter}, f::ZPKFilter{Z,P})
    z = f.z
    p = f.p
    n = max(length(z), length(p))

    # Sort poles lexicographically so that matched poles are adjacent
    p = sort(p, lt=lexlt)

    # Sort poles according to distance to unit circle (farthest first)
    p = sort!(p, by=x->abs(abs(x) - 1), rev=true)

    # Move real poles to the end
    complexp = P[]
    realp = P[]
    for x in p
        push!(ifelse(imag(x) == zero(P), realp, complexp), x)
    end
    append!(complexp, realp)
    p = complexp

    # Group each pole with a zero
    zorder = zeros(Int, length(z))
    for i = 1:length(z)
        closest_idx = 1
        closest_val = Inf
        for j = 1:length(p)
            zorder[j] == 0 || continue
            val = abs(z[j] - p[i])
            if val < closest_val
                closest_idx = j
                closest_val = val
            end
        end
        zorder[closest_idx] = i
    end

    # Build second-order sections
    T = promote_type(realtype(Z), realtype(P))
    biquads = Array(BiquadFilter{T}, div(n, 2)+isodd(n))
    for i = 1:div(n, 2)
        biquads[i] = convert(BiquadFilter, ZPKFilter(z[2i-1 .<= zorder .<= 2i], p[2i-1:2i], one(T)))
    end
    if isodd(n)
        biquads[end] = convert(BiquadFilter, ZPKFilter([z[zorder .== length(p)]], [p[end]], one(T)))
    end

    SOSFilter(biquads, f.k)
end

convert(::Type{SOSFilter}, f::Filter) = convert(SOSFilter, convert(ZPKFilter, f))

filt(f::Filter, x) = filt(convert(TFFilter, f), x)
filt(f::TFFilter, x) = filt(coeffs(f.b), coeffs(f.a), x)

function filt(f::SOSFilter, x::AbstractVector)
    y = copy(x)
    biquads = f.biquads
    si = zeros(2, length(biquads))

    @inbounds begin
        for i = 1:size(x, 1), fi = 1:length(biquads)
            biquad = biquads[fi]
            yp = y[i]
            y[i] = si[1, fi] + biquad.b0*yp
            si[1, fi] = si[2, fi] + biquad.b1*yp - biquad.a1*y[i]
            si[2, fi] = biquad.b2*yp - biquad.a2*y[i]
        end
    end
    scale!(y, f.g)
end

abstract FilterType

immutable Lowpass{T} <: FilterType
    w::T
end

immutable Highpass{T} <: FilterType
    w::T
end

immutable Bandpass{T} <: FilterType
    w1::T
    w2::T
end

immutable Bandstop{T} <: FilterType
    w1::T
    w2::T
end

function Butterworth(N::Integer)
    poles = zeros(Complex128, N)
    for i = 1:div(N, 2)
        w = (2*i-1)/2N
        pole = complex(-sinpi(w), cospi(w))
        poles[i*2-1] = pole
        poles[i*2] = conj(pole)
    end
    if isodd(N)
        poles[end] = -1.0+0.0im
    end
    ZPKFilter(Float64[], poles, 1)
end

# Create a lowpass filter from a lowpass filter prototype
function transform_prototype(ftype::Lowpass, proto::TFFilter)
    b = lfill(coeffs(proto.b), length(proto.a))
    a = lfill(coeffs(proto.a), length(proto.b))
    c = ftype.w .^ [length(a)-1:-1:0]
    TFFilter(b ./ c, a ./ c)
end

# Create a highpass filter from a lowpass filter prototype
function transform_prototype(ftype::Highpass, proto::TFFilter)
    b = lfill(coeffs(proto.b), length(proto.a))
    a = lfill(coeffs(proto.a), length(proto.b))
    c = ftype.w .^ [0:length(a)-1]
    TFFilter(flipud(b) .* c, flipud(a) .* c)
end

# Create a bandpass filter from a lowpass filter prototype
# Thus is a direct port of Scipy's lp2bp
function transform_prototype(ftype::Bandpass, proto::TFFilter)
    bw = ftype.w2 - ftype.w1
    wo = sqrt(ftype.w1 * ftype.w2)
    tf = convert(TFFilter, proto)
    b = tf.b.a
    a = tf.a.a
    D = length(a) - 1
    N = length(b) - 1
    M = max(N, D)
    Np = N + M
    Dp = D + M
    bprime = similarzeros(b, Np+1)
    aprime = similarzeros(a, Dp+1)
    wosq = wo^2
    for j = 0:Np
        val = 0.0
        for i = 0:N
            for k = 0:i
                if M - i + 2 * k == j
                    val += binomial(i, k) * b[N - i + 1] * wosq^(i - k) / bw^i
                end
            end
        end

        bprime[Np - j + 1] = val
    end
    for j = 0:Dp
        val = 0.0
        for i = 0:D
            for k in 0:i+1
                if M - i + 2 * k == j
                    val += binomial(i, k) * a[D - i + 1] * wosq ^(i - k) / bw^i
                end
            end
        end
        aprime[Dp - j + 1] = val
    end
    TFFilter(bprime, aprime)
end

# Create a bandstop filter from a lowpass filter prototype
# Thus is a direct port of Scipy's lp2bs
function transform_prototype(ftype::Bandstop, proto::TFFilter)
    bw = ftype.w2 - ftype.w1
    wo = sqrt(ftype.w1 * ftype.w2)
    tf = convert(TFFilter, proto)
    b = tf.b.a
    a = tf.a.a
    D = length(a) - 1
    N = length(b) - 1
    M = max(N, D)
    Np = 2 * M
    Dp = 2 * M
    bprime = similarzeros(b, Np+1)
    aprime = similarzeros(a, Dp+1)
    wosq = wo^2
    for j = 0:Np
        val = 0.0
        for i = 0:N
            for k = 0:M-i
                if i + 2 * k == j
                    val += binomial(M - i, k) * b[N - i + 1] * wosq^(M - i - k) * bw^i
                end
            end
        end
        bprime[Np - j + 1] = val
    end
    for j = 0:Dp
        val = 0.0
        for i = 0:D
            for k in 0:M-i
                if i + 2 * k == j
                    val += binomial(M - i, k) * a[D - i + 1] * wosq^(M - i - k) * bw^i
                end
            end
        end
        aprime[Dp - j + 1] = val
    end
    TFFilter(bprime, aprime)
end

transform_prototype(ftype::FilterType, proto::Filter) =
    transform_prototype(ftype, convert(TFFilter, proto))

analogfilter(ftype::FilterType, proto::Filter) = transform_prototype(ftype, proto)

# Do bilinear transform
bilinear(f::Filter, fs::Real) = bilinear(convert(ZPKFilter, f), fs)
bilinear(f::ZPKFilter, fs::Real) =
    ZPKFilter(lfill((2 .+ f.z / fs)./(2 .- f.z / fs), length(f.p), -1),
    (2 .+ f.p / fs)./(2 .- f.p / fs), real(f.k * prod(2 * fs .- f.z) ./ prod(2 * fs .- f.p)))

# Pre-warp filter frequencies for digital filtering
prewarp(ftype::Union(Lowpass, Highpass)) = (typeof(ftype))(4*tan(pi*ftype.w/2))
prewarp(ftype::Union(Bandpass, Bandstop)) = (typeof(ftype))(4*tan(pi*ftype.w1/2), 4*tan(pi*ftype.w2/2))

# Digital filter design using ZPKFilter->TFFilter->ZPKFilter conversion on all poles
digitalfilter(ftype::FilterType, proto::Filter, as::Type{ZPKFilter}=ZPKFilter) =
    bilinear(transform_prototype(prewarp(ftype), proto), 2)

# Digital filter design using second-order sections
function digitalfilter(ftype::FilterType, proto::SOSFilter)
    ftype = prewarp(ftype)
    g = proto.g
    biquads = vcat([begin
                        s = convert(SOSFilter, bilinear(transform_prototype(ftype, f), 2))
                        g *= s.g
                        s.biquads
                    end for f in proto.biquads]...)
    SOSFilter(biquads, g)
end
end
