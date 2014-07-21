# Filter prototypes, transformations, and transforms

abstract FilterType

#
# Prototypes
#

function Butterworth(T::Type, n::Integer)
    poles = zeros(Complex{T}, n)
    for i = 1:div(n, 2)
        w = convert(T, 2i-1)/2n
        pole = complex(-sinpi(w), cospi(w))
        poles[2i-1] = pole
        poles[2i] = conj(pole)
    end
    if isodd(n)
        poles[end] = -1
    end
    ZPKFilter(T[], poles, 1)
end
Butterworth(n::Integer) = Butterworth(Float64, n)

function chebyshev_poles(T::Type, n::Integer, ε::Real)
    p = zeros(Complex{T}, n)
    μ = asinh(convert(T, 1)/ε)/n
    b = -sinh(μ)
    c = cosh(μ)
    for i = 1:div(n, 2)
        w = convert(T, 2i-1)/2n
        pole = complex(b*sinpi(w), c*cospi(w))
        p[2i-1] = pole
        p[2i] = conj(pole)
    end
    if isodd(n)
        w = convert(T, 2*div(n, 2)+1)/2n
        pole = b*sinpi(w)
        p[end] = pole
    end
    p
end

function Chebyshev1(T::Type, n::Integer, ripple::Real)
    ε = sqrt(10^(convert(T, ripple)/10)-1)
    p = chebyshev_poles(T, n, ε)
    k = one(T)
    for i = 1:div(n, 2)
        k *= abs2(p[2i])
    end
    if iseven(n)
        k /= sqrt(1+abs2(ε))
    else
        k *= real(-p[end])
    end
    ZPKFilter(Float64[], p, k)
end
Chebyshev1(n::Integer, ripple::Real) = Chebyshev1(Float64, n, ripple)

function Chebyshev2(T::Type, n::Integer, ripple::Real)
    ε = 1/sqrt(10^(convert(T, ripple)/10)-1)
    p = chebyshev_poles(T, n, ε)
    for i = 1:length(p)
        p[i] = inv(p[i])
    end

    z = zeros(Complex{T}, n-isodd(n))
    k = one(T)
    for i = 1:div(n, 2)
        w = convert(T, 2i-1)/2n
        ze = Complex(zero(T), -inv(cospi(w)))
        z[2i-1] = ze
        z[2i] = conj(ze)
        k *= abs2(p[2i])/abs2(ze)
    end
    isodd(n) && (k *= -real(p[end]))

    ZPKFilter(z, p, k)
end
Chebyshev2(n::Integer, ripple::Real) = Chebyshev2(Float64, n, ripple)

#
# Prototype transformations
#

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

# Create a lowpass filter from a lowpass filter prototype
function transform_prototype(ftype::Lowpass, proto::TFFilter)
    TFFilter(Poly([proto.b[i]/ftype.w^(i) for i = 0:length(proto.b)-1]),
             Poly([proto.a[i]/ftype.w^(i) for i = 0:length(proto.a)-1]))
end

# Create a highpass filter from a lowpass filter prototype
function transform_prototype(ftype::Highpass, proto::TFFilter)
    n = max(length(proto.b), length(proto.a))
    TFFilter(Poly([proto.b[n-i-1]/ftype.w^(i) for i = 0:n-1]),
             Poly([proto.a[n-i-1]/ftype.w^(i) for i = 0:n-1]))
end

# Create a bandpass filter from a lowpass filter prototype
# Thus is a direct port of Scipy's lp2bp
function transform_prototype(ftype::Bandpass, proto::TFFilter)
    bw = ftype.w2 - ftype.w1
    wo = sqrt(ftype.w1 * ftype.w2)
    b = proto.b
    a = proto.a
    D = length(a) - 1
    N = length(b) - 1
    M = max(N, D)
    Np = N + M
    Dp = D + M
    bprime = zeros(eltype(b), Np+1)
    aprime = zeros(eltype(a), Dp+1)
    wosq = wo^2
    for j = 0:Np
        val = 0.0
        for i = 0:N
            for k = 0:i
                if M - i + 2 * k == j
                    val += binomial(i, k) * b[i] * wosq^(i - k) / bw^i
                end
            end
        end

        bprime[j+1] = val
    end
    for j = 0:Dp
        val = 0.0
        for i = 0:D
            for k in 0:i+1
                if M - i + 2 * k == j
                    val += binomial(i, k) * a[i] * wosq ^(i - k) / bw^i
                end
            end
        end
        aprime[j+1] = val
    end
    TFFilter(Poly(bprime), Poly(aprime))
end

# Create a bandstop filter from a lowpass filter prototype
# Thus is a direct port of Scipy's lp2bs
function transform_prototype(ftype::Bandstop, proto::TFFilter)
    bw = ftype.w2 - ftype.w1
    wo = sqrt(ftype.w1 * ftype.w2)
    b = proto.b
    a = proto.a
    D = length(a) - 1
    N = length(b) - 1
    M = max(N, D)
    Np = 2 * M
    Dp = 2 * M
    bprime = zeros(eltype(b), Np+1)
    aprime = zeros(eltype(a), Dp+1)
    wosq = wo^2
    for j = 0:Np
        val = 0.0
        for i = 0:N
            for k = 0:M-i
                if i + 2 * k == j
                    val += binomial(M - i, k) * b[i] * wosq^(M - i - k) * bw^i
                end
            end
        end
        bprime[j+1] = val
    end
    for j = 0:Dp
        val = 0.0
        for i = 0:D
            for k in 0:M-i
                if i + 2 * k == j
                    val += binomial(M - i, k) * a[i] * wosq^(M - i - k) * bw^i
                end
            end
        end
        aprime[j+1] = val
    end
    TFFilter(Poly(bprime), Poly(aprime))
end

transform_prototype(ftype::FilterType, proto::Filter) =
    transform_prototype(ftype, convert(TFFilter, proto))

analogfilter(ftype::FilterType, proto::Filter) = transform_prototype(ftype, proto)

# Do bilinear transform
bilinear(f::Filter, fs::Real) = bilinear(convert(ZPKFilter, f), fs)
function bilinear{Z,P,K}(f::ZPKFilter{Z,P,K}, fs::Real)
    ztype = typeof(0 + zero(Z)/fs)
    z = fill(convert(ztype, -1), max(length(f.p), length(f.z)))

    ptype = typeof(0 + zero(P)/fs)
    p = Array(typeof(zero(P)/fs), length(f.p))

    num = one(one(fs) - one(Z))
    for i = 1:length(f.z)
        z[i] = (2 + f.z[i] / fs)/(2 - f.z[i] / fs)
        num *= (2 * fs - f.z[i])
    end

    den = one(one(fs) - one(P))
    for i = 1:length(f.p)
        p[i] = (2 + f.p[i] / fs)/(2 - f.p[i]/fs)
        den *= (2 * fs - f.p[i])
    end

    ZPKFilter(z, p, f.k * real(num/den))
end

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
