# Filter prototypes, transformations, and transforms

abstract FilterType

#
# Butterworth prototype
#

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
