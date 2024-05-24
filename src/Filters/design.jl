# Filter prototypes, transformations, and transforms

using ..Windows

abstract type FilterType end

#
# Butterworth prototype
#

function Butterworth(::Type{T}, n::Integer) where {T<:Real}
    n > 0 || throw(DomainError(n, "n must be positive"))

    poles = Vector{Complex{T}}(undef, n)
    for i = 1:n÷2
        w = convert(T, 2i - 1) / 2n
        sinpi_w, cospi_w = @inline sincospi(w)
        pole = complex(-sinpi_w, cospi_w)
        poles[2i-1] = pole
        poles[2i] = conj(pole)
    end
    if isodd(n)
        poles[end] = -1
    end
    ZeroPoleGain{:s}(Complex{T}[], poles, one(T))
end

"""
    Butterworth(n::Integer)

``n`` pole Butterworth filter.
"""
Butterworth(n::Integer) = Butterworth(Float64, n)

#
# Chebyshev type I and II prototypes
#

function chebyshev_poles(::Type{T}, n::Integer, ε::Real) where {T<:Real}
    p = Vector{Complex{T}}(undef, n)
    μ = asinh(convert(T, 1)/ε)/n
    b = -sinh(μ)
    c = cosh(μ)
    for i = 1:n÷2
        w = convert(T, 2i - 1) / 2n
        sinpi_w, cospi_w = @inline sincospi(w)
        pole = complex(b * sinpi_w, c * cospi_w)
        p[2i-1] = pole
        p[2i] = conj(pole)
    end
    if isodd(n)
        w = convert(T, 2 * (n ÷ 2) + 1) / 2n
        pole = b * sinpi(w)
        p[end] = pole
    end
    p
end

function Chebyshev1(::Type{T}, n::Integer, ripple::Real) where {T<:Real}
    n > 0 || throw(DomainError(n, "n must be positive"))
    ripple >= 0 || throw(DomainError(ripple, "ripple must be non-negative"))

    ε = sqrt(10^(convert(T, ripple)/10)-1)
    p = chebyshev_poles(T, n, ε)
    k = one(T)
    for i = 1:n÷2
        k *= abs2(p[2i])
    end
    if iseven(n)
        k /= sqrt(1+abs2(ε))
    else
        k *= real(-p[end])
    end
    ZeroPoleGain{:s}(Complex{T}[], p, k)
end

"""
    Chebyshev1(n::Integer, ripple::Real)

`n` pole Chebyshev type I filter with `ripple` dB ripple in
the passband.
"""
Chebyshev1(n::Integer, ripple::Real) = Chebyshev1(Float64, n, ripple)

function Chebyshev2(::Type{T}, n::Integer, ripple::Real) where {T<:Real}
    n > 0 || throw(DomainError(n, "n must be positive"))
    ripple >= 0 || throw(DomainError(ripple, "ripple must be non-negative"))

    ε = 1/sqrt(10^(convert(T, ripple)/10)-1)
    p = chebyshev_poles(T, n, ε)
    map!(inv, p, p)

    z = Vector{Complex{T}}(undef, n - isodd(n))
    k = one(T)
    for i = 1:n÷2
        w = convert(T, 2i - 1) / 2n
        ze = @inline complex(zero(T), -inv(cospi(w)))
        z[2i-1] = ze
        z[2i] = conj(ze)
        k *= abs2(p[2i]) / abs2(ze)
    end
    isodd(n) && (k *= -real(p[end]))

    ZeroPoleGain{:s}(z, p, k)
end

"""
    Chebyshev2(n::Integer, ripple::Real)

`n` pole Chebyshev type II filter with `ripple` dB ripple in
the stopband.
"""
Chebyshev2(n::Integer, ripple::Real) = Chebyshev2(Float64, n, ripple)

#
# Elliptic prototype
#
# See Orfanidis, S. J. (2007). Lecture notes on elliptic filter design.
# Retrieved from http://www.ece.rutgers.edu/~orfanidi/ece521/notes.pdf

# Compute Landen sequence for evaluation of elliptic functions
function landen(k::Real)
    niter = 7
    kn = Vector{typeof(k)}(undef, niter)
    # Eq. (50)
    for i = 1:niter
        kn[i] = k = abs2(k/(1+sqrt(1-abs2(k))))
    end
    kn
end

# cde computes cd(u*K(k), k)
# sne computes sn(u*K(k), k)
# Both accept the Landen sequence as generated by landen above
function _ellip(init::Number, landen::Vector{<:Real})
    winv = inv(init)
    # Eq. (55)
    for x in Iterators.reverse(landen)
        winv = 1 / (1 + x) * (winv + x / winv)
    end
    w = inv(winv)
end
@inline cde(u::Number, landen::Vector{<:Real}) = @inline _ellip(cospi(u / 2), landen)
@inline sne(u::Number, landen::Vector{<:Real}) = @inline _ellip(sinpi(u / 2), landen)

# sne inverse
function asne(w::Number, k::Real)
    oldw = NaN
    while w != oldw
        oldw = w
        kold = k
        # Eq. (50)
        k = abs2(k/(1+sqrt(1-abs2(k))))
        # Eq. (56)
        w = 2w / ((1 + k) * (1 + sqrt(muladd(-abs2(kold), w^2, 1))))
    end
    2*asin(w)/π
end

function Elliptic(::Type{T}, n::Integer, rp::Real, rs::Real) where {T<:Real}
    n > 0 || throw(DomainError(n, "n must be positive"))
    rp > 0 || throw(DomainError(rp, "rp must be positive"))
    rp < rs || throw(DomainError(rp, "rp must be less than rs"))

    # Eq. (2)
    εp = sqrt(10^(convert(T, rp)/10)-1)
    εs = sqrt(10^(convert(T, rs)/10)-1)

    # Eq. (3)
    k1 = εp/εs
    k1 >= 1 && throw(ArgumentError("filter order is too high for parameters"))

    # Eq. (20)
    k1′² = 1 - abs2(k1)
    k1′ = sqrt(k1′²)
    k1′_landen = landen(k1′)

    # Eq. (47)
    k′ = one(T)
    for i = 1:n÷2
        k′ *= sne(convert(T, 2i-1)/n, k1′_landen)
    end
    k′ = k1′²^(convert(T, n)/2)*k′^4

    k = sqrt(1 - abs2(k′))
    k_landen = landen(k)

    # Eq. (65)
    v0 = -im/convert(T, n)*asne(im/εp, k1)

    z = Vector{Complex{T}}(undef, 2*(n÷2))
    p = Vector{Complex{T}}(undef, n)
    gain = one(T)
    for i = 1:n÷2
        # Eq. (43)
        w = convert(T, 2i-1)/n

        # Eq. (62)
        ze = complex(zero(T), -inv(k*cde(w, k_landen)))
        z[2i-1] = ze
        z[2i] = conj(ze)

        # Eq. (64)
        pole = im*cde(w - im*v0, k_landen)
        p[2i-1] = conj(pole)
        p[2i] = pole

        gain *= abs2(pole)/abs2(ze)
    end

    if isodd(n)
        pole = im*sne(im*v0, k_landen)
        p[end] = pole
        gain *= abs(pole)
    else
        gain *= 10^(-convert(T, rp)/20)
    end

    ZeroPoleGain{:s}(z, p, gain)
end

"""
    Elliptic(n::Integer, rp::Real, rs::Real)

`n` pole elliptic (Cauer) filter with `rp` dB ripple in the
passband and `rs` dB attentuation in the stopband.
"""
Elliptic(n::Integer, rp::Real, rs::Real) = Elliptic(Float64, n, rp, rs)

#
# Prototype transformation types
#

# returns frequency in half-cycles per sample ∈ (0, 1)
function normalize_freq(w::Real, fs::Real)
    w <= 0 && throw(DomainError(w, "frequencies must be positive"))
    f = 2 * w / fs
    f >= 1 && throw(DomainError(f, "frequencies must be less than the Nyquist frequency $(fs/2)"))
    f
end
function normalize_complex_freq(w::Real, fs::Real)
    f = 2 * w / fs
    f >= 2 && throw(DomainError(f, "frequencies must be less than the Nyquist frequency $(fs)"))
    f
end

"""
    Lowpass(Wn::Real)

Low pass filter with cutoff frequency `Wn`.
"""
struct Lowpass{T<:Real} <: FilterType
    w::T
    Lowpass{T}(w::Real) where {T<:Real} = new{T}(w)
    Lowpass(w::Real) = Lowpass{typeof(w / 1)}(w)
end

"""
    Highpass(Wn::Real)

High pass filter with cutoff frequency `Wn`.
"""
struct Highpass{T<:Real} <: FilterType
    w::T
    Highpass{T}(w::Real) where {T<:Real} = new{T}(w)
    Highpass(w::Real) = Highpass{typeof(w / 1)}(w)
end

"""
    Bandpass(Wn1::Real, Wn2::Real)

Band pass filter with pass band frequencies (`Wn1`, `Wn2`).
"""
struct Bandpass{T<:Real} <: FilterType
    w1::T
    w2::T
    function Bandpass{T}(w1::Real, w2::Real) where {T<:Real}
        w1 < w2 || throw(ArgumentError("w1 must be less than w2"))
        new{T}(w1, w2)
    end
    Bandpass(w1::T, w2::V) where {T<:Real,V<:Real} =
        Bandpass{typeof(one(promote_type(T, V)) / 1)}(w1, w2)
end

"""
    ComplexBandpass(Wn1, Wn2)

Complex band pass filter with pass band frequencies (`Wn1`, `Wn2`).
"""
struct ComplexBandpass{T<:Real} <: FilterType
    w1::T
    w2::T
    function ComplexBandpass{T}(w1::Real, w2::Real) where {T<:Real}
        w1 < w2 || throw(ArgumentError("w1 must be less than w2"))
        new{T}(w1, w2)
    end
    ComplexBandpass(w1::T, w2::V) where {T,V} =
        ComplexBandpass{typeof(one(promote_type(T, V)) / 1)}(w1, w2)
end

"""
    Bandstop(Wn1::Real, Wn2::Real)

Band stop filter with stop band frequencies (`Wn1`, `Wn2`).
"""
struct Bandstop{T<:Real} <: FilterType
    w1::T
    w2::T
    function Bandstop{T}(w1::Real, w2::Real) where {T<:Real}
        w1 < w2 || throw(ArgumentError("w1 must be less than w2"))
        new{T}(w1, w2)
    end
    Bandstop(w1::T, w2::V) where {T,V} =
        Bandstop{typeof(one(promote_type(T, V)) / 1)}(w1, w2)
end

#
# Prototype transformation implementations
#
# The formulas implemented here come from the documentation for the
# corresponding functionality in Octave, available at
# https://staff.ti.bfh.ch/sha1/Octave/index/f/sftrans.html
# The Octave implementation was not consulted in creating this code.

# Create a lowpass filter from a lowpass filter prototype
transform_prototype(ftype::Lowpass, proto::ZeroPoleGain{:s, Z, P, K}) where {Z, P, K} =
    ZeroPoleGain{:s, Z, P, K}(ftype.w * proto.z, ftype.w * proto.p,
              proto.k * ftype.w^(length(proto.p)-length(proto.z)))

# Create a highpass filter from a lowpass filter prototype
function transform_prototype(ftype::Highpass, proto::ZeroPoleGain{:s})
    z = proto.z
    p = proto.p
    k = proto.k
    nz = length(z)
    np = length(p)
    TR = Base.promote_eltype(z, p)
    newz = zeros(TR, max(nz, np))
    newp = zeros(TR, max(nz, np))

    num = one(eltype(z))
    for i = 1:nz
        num *= -z[i]
        newz[i] = ftype.w / z[i]
    end

    den = one(eltype(p))
    for i = 1:np
        den *= -p[i]
        newp[i] = ftype.w / p[i]
    end

    abs(real(num) - 1) < np*eps(real(num)) && (num = 1)
    abs(real(den) - 1) < np*eps(real(den)) && (den = 1)
    ZeroPoleGain{:s}(newz, newp, oftype(k, k * real(num)/real(den)))
end

# Create a bandpass filter from a lowpass filter prototype
function transform_prototype(ftype::Bandpass, proto::ZeroPoleGain{:s})
    z = proto.z
    p = proto.p
    k = proto.k
    nz = length(z)
    np = length(p)
    ncommon = min(nz, np)
    TR = Base.promote_eltype(z, p)
    newz = zeros(TR, 2*nz+np-ncommon)
    newp = zeros(TR, 2*np+nz-ncommon)
    for (oldc, newc) in ((p, newp), (z, newz))
        for i in eachindex(oldc)
            b = oldc[i] * ((ftype.w2 - ftype.w1)/2)
            pm = sqrt(muladd(-ftype.w2, ftype.w1, b^2))
            newc[2i-1] = b + pm
            newc[2i] = b - pm
        end
    end
    ZeroPoleGain{:s}(newz, newp, oftype(k, k * (ftype.w2 - ftype.w1) ^ (np - nz)))
end

# Create a bandstop filter from a lowpass filter prototype
function transform_prototype(ftype::Bandstop, proto::ZeroPoleGain{:s})
    z = proto.z
    p = proto.p
    k = proto.k
    nz = length(z)
    np = length(p)
    npairs = max(nz, np)
    TR = Base.promote_eltype(z, p)
    newz = Vector{TR}(undef, 2*npairs)
    newp = Vector{TR}(undef, 2*npairs)

    num = one(eltype(z))
    for i = 1:nz
        num *= -z[i]
        b = (ftype.w2 - ftype.w1)/2z[i]
        pm = sqrt(muladd(-ftype.w2, ftype.w1, b^2))
        newz[2i-1] = b - pm
        newz[2i] = b + pm
    end

    den = one(eltype(p))
    for i = 1:np
        den *= -p[i]
        b = (ftype.w2 - ftype.w1)/2p[i]
        pm = sqrt(muladd(-ftype.w2, ftype.w1, b^2))
        newp[2i-1] = b - pm
        newp[2i] = b + pm
    end

    # Any emaining poles/zeros are real and not cancelled
    npm = sqrt(-complex(ftype.w2 * ftype.w1))
    for (n, newc) in ((np, newp), (nz, newz))
        for i = n+1:npairs
            newc[2i-1] = -npm
            newc[2i] = npm
        end
    end

    abs(real(num) - 1) < np*eps(real(num)) && (num = 1)
    abs(real(den) - 1) < np*eps(real(den)) && (den = 1)
    ZeroPoleGain{:s}(newz, newp, oftype(k, k * real(num)/real(den)))
end

transform_prototype(ftype, proto::FilterCoefficients{:s}) =
    transform_prototype(ftype, convert(ZeroPoleGain, proto))

"""
    analogfilter(responsetype::FilterType, designmethod::FilterCoefficients)

Construct an analog filter. See below for possible response and
filter types.
"""
analogfilter(ftype::FilterType, proto::FilterCoefficients) =
    transform_prototype(ftype, proto)

# Bilinear transform
"""
    bilinear(f::FilterCoefficients{:s}, fs::Real)

Calculate the digital filter (z-domain) ZPK representation of an analog filter defined
in s-domain using bilinear transform with sampling frequency `fs`. The s-domain
representation is first converted to a ZPK representation in s-domain and then
transformed to z-domain using bilinear transform.
"""
bilinear(f::FilterCoefficients{:s}, fs::Real) = bilinear(convert(ZeroPoleGain, f), fs)

"""
    bilinear(f::ZeroPoleGain{:s,Z,P,K}, fs::Real) where {Z,P,K}

Calculate the digital filter (z-domain) ZPK representation of an analog filter defined
as a ZPK representation in s-domain using bilinear transform with sampling frequency
`fs`.

Input s-domain representation must be a `ZeroPoleGain{:s, Z, P, K}` object:
```math
H(s) = f.k\\frac{(s - \\verb!f.z[1]!) \\ldots (s - \\verb!f.z[m]!)}{(s - \\verb!f.p[1]!) \\ldots (s - \\verb!f.p[n]!)}
```
Output z-domain representation is a `ZeroPoleGain{:z, Z, P, K}` object:
```math
H(z) = K\\frac{(z - \\verb!Z[1]!) \\ldots (z - \\verb!Z[m]!)}{(z - \\verb!P[1]!) \\ldots (z - \\verb!P[n]!)}
```
where `Z, P, K` are calculated as:
```math
Z[i] = \\frac{(2 + \\verb!f.z[i]!/\\verb!fs!)}{(2 - \\verb!f.z[i]!/\\verb!fs!)} \\quad \\text{for } i = 1, \\ldots, m
```
```math
P[i] = \\frac{(2 + \\verb!f.p[i]!/\\verb!fs!)}{(2 - \\verb!f.p[i]!/\\verb!fs!)} \\quad \\text{for } i = 1, \\ldots, n
```
```math
K = f.k \\ \\mathcal{Re} \\left[ \\frac{\\prod_{i=1}^m (2*fs - f.z[i])}{\\prod_{i=1}^n (2*fs - f.p[i])} \\right]
```
Here, `m` and `n` are respectively the numbers of zeros and poles in the s-domain representation. If `m < n`,
then additional `n-m` zeros are added at `z = -1`.
"""
function bilinear(f::ZeroPoleGain{:s,Z,P,K}, fs::Real) where {Z,P,K}
    ztype = typeof(0 + zero(Z)/fs)
    z = fill(convert(ztype, -1), max(length(f.p), length(f.z)))

    ptype = typeof(0 + zero(P)/fs)
    p = Vector{typeof(zero(P)/fs)}(undef, length(f.p))

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

    ZeroPoleGain{:z}(z, p, f.k * real(num)/real(den))
end

# Pre-warp filter frequencies for digital filtering
prewarp(ftype::Lowpass, fs::Real) = Lowpass(prewarp(normalize_freq(ftype.w, fs)))
prewarp(ftype::Highpass, fs::Real) = Highpass(prewarp(normalize_freq(ftype.w, fs)))
prewarp(ftype::Bandpass, fs::Real) = Bandpass(prewarp(normalize_freq(ftype.w1, fs)), prewarp(normalize_freq(ftype.w2, fs)))
prewarp(ftype::Bandstop, fs::Real) = Bandstop(prewarp(normalize_freq(ftype.w1, fs)), prewarp(normalize_freq(ftype.w2, fs)))
# freq in half-samples per cycle
prewarp(f::Real) = 4*tan(pi*f/2)

# Digital filter design
"""
    digitalfilter(responsetype::FilterType, designmethod::FilterCoefficients[; fs::Real])

Construct a digital filter. See below for possible response and
filter types.
"""
digitalfilter(ftype::FilterType, proto::FilterCoefficients; fs::Real=2) =
    bilinear(transform_prototype(prewarp(ftype, fs), proto), 2)

#
# Special filter types
#

"""
    iirnotch(Wn::Real, bandwidth::Real[; fs=2])

Second-order digital IIR notch filter [^Orfanidis] at frequency `Wn` with
bandwidth `bandwidth`. If `fs` is not specified, `Wn` is
interpreted as a normalized frequency in half-cycles/sample.

[^Orfanidis]: Orfanidis, S. J. (1996). Introduction to signal processing.
    Englewood Cliffs, N.J: Prentice Hall, p. 370.
"""
function iirnotch(w::Real, bandwidth::Real; fs=2)
    w = normalize_freq(w, fs)
    bandwidth = normalize_freq(bandwidth, fs)

    # Eq. 8.2.23
    b = 1/(1+tan(pi*bandwidth/2))
    # Eq. 8.2.22
    cosw0 = cospi(w)
    Biquad(b, -2b*cosw0, b, -2b*cosw0, 2b-1)
end

#
# FIR filter design
#

# Get length and alpha for Kaiser window filter with specified
# transition band width and stopband attenuation in dB
function kaiserord(transitionwidth::Real, attenuation::Real=60)
    n = ceil(Int, (attenuation - 7.95)/(π*2.285*transitionwidth))+1

    if attenuation > 50
        β = 0.1102*(attenuation - 8.7)
    elseif attenuation >= 21
        β = 0.5842*(attenuation - 21)^0.4 + 0.07886*(attenuation - 21)
    else
        β = 0.0
    end

    return n, β/π
end

struct FIRWindow{T}
    window::Vector{T}
    scale::Bool
end

"""
    FIRWindow(window::Vector; scale=true)

FIR filter design using window `window`, a vector whose length
matches the number of taps in the resulting filter.

If `scale` is `true` (default), the designed FIR filter is
scaled so that the following holds:

- For [`Lowpass`](@ref) and [`Bandstop`](@ref) filters, the frequency
  response is unity at 0 (DC).
- For [`Highpass`](@ref) filters, the frequency response is unity
  at the Nyquist frequency.
- For [`Bandpass`](@ref) filters, the frequency response is unity
  in the center of the passband.
"""
FIRWindow(window::Vector; scale::Bool=true) = FIRWindow(window, scale)

# FIRWindow(n::Integer, window::Function, args...) = FIRWindow(window(n, args...))
"""
    FIRWindow(; transitionwidth::Real, attenuation::Real=60, scale::Bool=true)

Kaiser window FIR filter design. The required number of taps is
calculated based on `transitionwidth` (in half-cycles/sample)
and stopband `attenuation` (in dB). `attenuation` defaults to
60 dB.
"""
FIRWindow(; transitionwidth::Real=throw(ArgumentError("must specify transitionwidth")),
          attenuation::Real=60, scale::Bool=true) =
    FIRWindow(kaiser(kaiserord(transitionwidth, attenuation)...), scale)

# Compute coefficients for FIR prototype with specified order
function _firprototype(n::Integer, ftype::Lowpass, fs::Real, ::Type{T}) where {T<:Number}
    w = normalize_freq(ftype.w, fs)

    promote_type(typeof(w), T)[w*sinc(w*(k-(n+1)/2)) for k = 1:n]
end

firprototype(n::Integer, ftype::Lowpass, fs::Real) =
    _firprototype(n, ftype, fs, typeof(fs))

function firprototype(n::Integer, ftype::Bandpass, fs::Real)
    w1 = normalize_freq(ftype.w1, fs)
    w2 = normalize_freq(ftype.w2, fs)

    [w2*sinc(w2*(k-(n+1)/2)) - w1*sinc(w1*(k-(n+1)/2)) for k = 1:n]
end

function firprototype(n::Integer, ftype::ComplexBandpass, fs::T) where {T<:Real}
    w1 = normalize_complex_freq(ftype.w1, fs)
    w2 = normalize_complex_freq(ftype.w2, fs)
    w_center = (w2 + w1) / 2
    w_cutoff = (w2 - w1) / 2
    lp = Lowpass(w_cutoff)
    _firprototype(n, lp, 2, Complex{T}) .*= cispi.(w_center * (0:(n-1)))
end

function firprototype(n::Integer, ftype::Highpass, fs::Real)
    w = normalize_freq(ftype.w, fs)
    isodd(n) || throw(ArgumentError("FIRWindow highpass filters must have an odd number of coefficients"))

    out = [-w*sinc(w*(k-(n+1)/2)) for k = 1:n]
    out[n÷2+1] += 1
    out
end

function firprototype(n::Integer, ftype::Bandstop, fs::Real)
    w1 = normalize_freq(ftype.w1, fs)
    w2 = normalize_freq(ftype.w2, fs)
    isodd(n) || throw(ArgumentError("FIRWindow bandstop filters must have an odd number of coefficients"))

    out = [w1*sinc(w1*(k-(n+1)/2)) - w2*sinc(w2*(k-(n+1)/2)) for k = 1:n]
    out[n÷2+1] += 1
    out
end

scalefactor(coefs::Vector, ::Union{Lowpass, Bandstop}, fs::Real) = sum(coefs)
function scalefactor(coefs::Vector{T}, ::Highpass, fs::Real) where {T<:Number}
    c = zero(T)
    for k = 1:length(coefs)
        c += ifelse(isodd(k), coefs[k], -coefs[k])
    end
    c
end
function scalefactor(coefs::Vector{T}, ftype::Bandpass, fs::Real) where {T<:Number}
    n = length(coefs)
    freq = normalize_freq(middle(ftype.w1, ftype.w2), fs)
    c = zero(T)
    for k = 1:n
        c = muladd(coefs[k], cospi(freq * (k - (n + 1) / 2)), c)
    end
    c
end
function scalefactor(coefs::Vector{T}, ftype::ComplexBandpass, fs::Real) where T<:Number
    n = length(coefs)
    freq = normalize_complex_freq(middle(ftype.w1, ftype.w2), fs)
    c = zero(T)
    for k = 1:n
        c = muladd(coefs[k], cispi(-freq * (k - (n + 1) / 2)), c)
    end
    return abs(c)
end

function digitalfilter(ftype::FilterType, proto::FIRWindow; fs::Real=2)
    coefs = firprototype(length(proto.window), ftype, fs)
    @assert length(proto.window) == length(coefs)
    out = (coefs .*= proto.window)
    proto.scale ? rmul!(out, 1/scalefactor(out, ftype, fs)) : out
end


# Compute FIR coefficients necessary for arbitrary rate resampling
function resample_filter(rate::AbstractFloat, Nϕ = 32, rel_bw = 1.0, attenuation = 60)
    f_nyq       = rate >= 1.0 ? 1.0/Nϕ : rate/Nϕ
    cutoff      = f_nyq * rel_bw
    trans_width = cutoff * 0.2

    # Determine resampling filter order
    hLen, α = kaiserord(trans_width, attenuation)

    # Round the number of taps up to a multiple of Nϕ.
    # Otherwise the missing taps will be filled with 0.
    hLen = Nϕ * ceil(Int, hLen/Nϕ)

    # Ensure that the filter is an odd length
    if (iseven(hLen))
        hLen += 1
    end

    # Design filter
    h = digitalfilter(Lowpass(cutoff), FIRWindow(kaiser(hLen, α)))
    rmul!(h, Nϕ)
end

# Compute FIR coefficients necessary for rational rate resampling
function resample_filter(rate::Union{Integer,Rational}, rel_bw = 1.0, attenuation = 60)
    Nϕ          = numerator(rate)
    decimation  = denominator(rate)
    f_nyq       = min(1/Nϕ, 1/decimation)
    cutoff      = f_nyq * rel_bw
    trans_width = cutoff * 0.2

    # Determine resampling filter order
    hLen, α = kaiserord(trans_width, attenuation)

    # Round the number of taps up to a multiple of Nϕ (same as interpolation factor).
    # Otherwise the missing taps will be filled with 0.
    hLen = Nϕ * ceil(Int, hLen/Nϕ)

    # Ensure that the filter is an odd length
    if (iseven(hLen))
        hLen += 1
    end

    # Design filter
    h = digitalfilter(Lowpass(cutoff), FIRWindow(kaiser(hLen, α)))
    rmul!(h, Nϕ)
end
