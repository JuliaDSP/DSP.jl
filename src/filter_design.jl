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

export Butterworth, Lowpass, Highpass, Bandpass, Bandstop, analogfilter, digitalfilter
import Base: convert, filt

# Pad a vector with f on the left to a specified length
lfill{T}(v::Vector{T}, n::Integer, f::Number) = length(v) < n ? [fill(f, n - length(v)), v] : v
lfill{T}(v::Vector{T}, n::Integer) = lfill(v, n, 0)

# Make a new array of zeros of the same type as another vector
similarzeros{T}(v::Array{T}, args...) = zeros(T, args...)

# Get coefficients of a polynomial
coeffs{T}(p::Poly{T}) = p.a[1+p.nzfirst:end]

abstract Filter

immutable ZPKFilter{Z,P,K} <: Filter
    z::Vector{Z}
    p::Vector{P}
    k::K
end

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

function convert(::Type{ZPKFilter}, f::TFFilter)
    k = real(f.b[1])
    b = f.b / k
    z = roots(b)
    p = roots(f.a)
    ZPKFilter(z, p, k)
end

function convert(::Type{TFFilter}, f::ZPKFilter)
    b = f.k*poly(f.z)
    a = poly(f.p)
    TFFilter(Poly(real(b.a)), Poly(real(a.a)))
end

filt(f::Filter, x) = filt(convert(TFFilter, f), x)
filt(f::TFFilter, x) = filt(coeffs(f.b), coeffs(f.a), x)

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

transform_prototype(ftype::Union(Lowpass, Highpass, Bandpass, Bandstop), proto::Filter) =
    transform_prototype(ftype, convert(TFFilter, proto))

# Do bilinear transform
bilinear(f::Filter, fs::Real) = bilinear(convert(ZPKFilter, f), fs)
bilinear(f::ZPKFilter, fs::Real) =
    ZPKFilter(lfill((2 + f.z / fs)./(2 - f.z / fs), length(f.p), -1),
    (2 + f.p / fs)./(2 - f.p / fs), real(f.k * prod(2 * fs - f.z) ./ prod(2 * fs - f.p)))

analogfilter(ftype::FilterType, proto::Filter) = transform_prototype(ftype, proto)
digitalfilter(ftype::Union(Lowpass, Highpass), proto::Filter) =
    bilinear(transform_prototype((typeof(ftype))(4*tan(pi*ftype.w/2)), proto), 2)
digitalfilter(ftype::Union(Bandpass, Bandstop), proto::Filter) =
    bilinear(transform_prototype((typeof(ftype))(4*tan(pi*ftype.w1/2), 4*tan(pi*ftype.w2/2)), proto), 2)
end
