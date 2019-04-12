# Filter types and conversions

abstract type FilterCoefficients{Domain} end

Base.convert(::Type{T}, f::FilterCoefficients) where {T<:FilterCoefficients} = T(f)


#
# Zero-pole gain form
#

"""
    ZeroPoleGain(z, p, k)

Filter representation in terms of zeros `z`, poles `p`, and
gain `k`:
```math
H(x) = k\\frac{(x - \\verb!z[1]!) \\ldots (x - \\verb!z[m]!)}{(x - \\verb!p[1]!) \\ldots (x - \\verb!p[n]!)}
```
"""
struct ZeroPoleGain{Domain,Z<:Number,P<:Number,K<:Number} <: FilterCoefficients{Domain}
    z::Vector{Z}
    p::Vector{P}
    k::K
end

ZeroPoleGain(f::FilterCoefficients{D}) where {D} = ZeroPoleGain{D}(f)
ZeroPoleGain(z::Vector, p::Vector, k) = ZeroPoleGain{:z}(z, p, k)
ZeroPoleGain{D,Z,P,K}(f::ZeroPoleGain) where {D,Z,P,K} = ZeroPoleGain{D,Z,P,K}(f.z, f.p, f.k)
ZeroPoleGain{D}(f::ZeroPoleGain{D,Z,P,K}) where {D,Z,P,K} = ZeroPoleGain{D,Z,P,K}(f)
ZeroPoleGain{D}(z::Vector{Z}, p::Vector{P}, k::K) where {D,Z<:Number,P<:Number,K<:Number} =
    ZeroPoleGain{D,Z,P,K}(z, p, k)

Base.promote_rule(::Type{ZeroPoleGain{D,Z1,P1,K1}}, ::Type{ZeroPoleGain{D,Z2,P2,K2}}) where {D,Z1,P1,K1,Z2,P2,K2} =
    ZeroPoleGain{D,promote_type(Z1,Z2),promote_type(P1,P2),promote_type(K1,K2)}

*(f::ZeroPoleGain{D}, g::Number) where {D} = ZeroPoleGain{D}(f.z, f.p, f.k*g)
*(g::Number, f::ZeroPoleGain{D}) where {D} = ZeroPoleGain{D}(f.z, f.p, f.k*g)
*(f1::ZeroPoleGain{D}, f2::ZeroPoleGain{D}) where {D} =
    ZeroPoleGain{D}([f1.z; f2.z], [f1.p; f2.p], f1.k*f2.k)
*(f1::ZeroPoleGain{D}, fs::ZeroPoleGain{D}...) where {D} =
    ZeroPoleGain{D}(vcat(f1.z, map(f -> f.z, fs)...), vcat(f1.p, map(f -> f.p, fs)...),
        f1.k*prod(f.k for f in fs))

#
# Transfer function form
#

struct PolynomialRatio{Domain,T<:Number} <: FilterCoefficients{Domain}
    b::Polynomial{T}
    a::Polynomial{T}

    PolynomialRatio{:z,Ti}(b::Polynomial, a::Polynomial) where {Ti<:Number} =
        new{:z,Ti}(convert(Polynomial{Ti}, b/a[end]), convert(Polynomial{Ti}, a/a[end]))
    PolynomialRatio{:s,Ti}(b::Polynomial, a::Polynomial) where {Ti<:Number} =
        new{:s,Ti}(convert(Polynomial{Ti}, b), convert(Polynomial{Ti}, a))
end
PolynomialRatio(f::FilterCoefficients{D}) where {D} = PolynomialRatio{D}(f)
"""
    PolynomialRatio(b, a)

Filter representation in terms of the coefficients of the numerator
`b` and denominator `a` of the transfer function:
```math
H(s) = \\frac{\\verb!b[1]! s^{m-1} + \\ldots + \\verb!b[m]!}{\\verb!a[1]! s^{n-1} + \\ldots + \\verb!a[n]!}
```
or equivalently:
```math
H(z) = \\frac{\\verb!b[1]! + \\ldots + \\verb!b[n]! z^{-n+1}}{\\verb!a[1]! + \\ldots + \\verb!a[n]! z^{-n+1}}
```
`b` and `a` may be specified as `Polynomial` objects or
vectors ordered from highest power to lowest.
"""
PolynomialRatio(b, a) = PolynomialRatio{:z}(b, a)
PolynomialRatio{D}(b::Polynomial{T}, a::Polynomial{T}) where {D,T<:Number} = PolynomialRatio{D,T}(b, a)
PolynomialRatio{D}(b::Polynomial, a::Polynomial) where {D} = PolynomialRatio{D}(promote(b, a)...)

# The DSP convention is highest power first. The Polynomials.jl
# convention is lowest power first.
function PolynomialRatio{:s,T}(b::Union{Number,Vector{<:Number}}, a::Union{Number,Vector{<:Number}}) where {T}
    if all(iszero, a)
        throw(ArgumentError("filter must have non-zero denominator"))
    end
    PolynomialRatio{:s,T}(Polynomial(reverse(b), :s), Polynomial(reverse(a), :s))
end
function PolynomialRatio{:z,T}(b::Union{Number,Vector{<:Number}}, a::Union{Number,Vector{<:Number}}) where {T}
    if all(iszero, a)
        throw(ArgumentError("filter must have non-zero denominator"))
    end
    b, a = pad_to_same_length(b, a)
    PolynomialRatio{:z,T}(Polynomial(reverse(b), :z), Polynomial(reverse(a), :z))
end
function pad_to_same_length(b, a)
    if length(a) < length(b)
        a = [a; zeros(eltype(a), length(b)-length(a))]
    elseif length(b) < length(a)
        b = [b; zeros(eltype(b), length(a)-length(b))]
    end
    return b, a
end
PolynomialRatio{D}(b::Union{T,Vector{T}}, a::Union{S,Vector{S}}) where {D,T<:Number,S<:Number} =
    PolynomialRatio{D,promote_type(T,S)}(b, a)

PolynomialRatio{D,T}(f::PolynomialRatio{D}) where {D,T} = PolynomialRatio{D,T}(f.b, f.a)
PolynomialRatio{D}(f::PolynomialRatio{D,T}) where {D,T} = PolynomialRatio{D,T}(f)

Base.promote_rule(::Type{PolynomialRatio{D,T}}, ::Type{PolynomialRatio{D,S}}) where {D,T,S} = PolynomialRatio{D,promote_type(T,S)}

function PolynomialRatio{D,T}(f::ZeroPoleGain{D}) where {D,T<:Real}
    b = f.k * fromroots(f.z)
    a = fromroots(f.p)
    PolynomialRatio{D,T}(Polynomial(real(coeffs(b))), Polynomial(real(coeffs(a))))
end
PolynomialRatio{D}(f::ZeroPoleGain{D,Z,P,K}) where {D,Z,P,K} =
    PolynomialRatio{D,promote_type(real(Z),real(P),K)}(f)

ZeroPoleGain{D,Z,P,K}(f::PolynomialRatio{D}) where {D,Z,P,K} =
    ZeroPoleGain{D,Z,P,K}(roots(f.b), roots(f.a), real(f.b[end]/f.a[end]))
function ZeroPoleGain{D}(f::PolynomialRatio{D,T}) where {D,T}
    z = roots(f.b)
    p = roots(f.a)
    # work around `roots` returning Any[] for polynomials of degree zero
    if isempty(z) && !(eltype(z) <: Number)
        z = typeof(one(T)/one(T))[]
    end
    if isempty(p) && !(eltype(p) <: Number)
        p = typeof(one(T)/one(T))[]
    end
    return ZeroPoleGain{D}(z, p, f.b[end]/f.a[end])
end

*(f::PolynomialRatio{D}, g::Number) where {D} = PolynomialRatio{D}(g*f.b, f.a)
*(g::Number, f::PolynomialRatio{D}) where {D} = PolynomialRatio{D}(g*f.b, f.a)
*(f1::PolynomialRatio{D}, f2::PolynomialRatio{D}) where {D} =
    PolynomialRatio{D}(f1.b*f2.b, f1.a*f2.a)
*(f1::PolynomialRatio{D}, fs::PolynomialRatio{D}...) where {D} =
    PolynomialRatio{D}(f1.b*prod(f.b for f in fs), f1.a*prod(f.a for f in fs))

"""
    coefb(f)

Coefficients of the numerator of a PolynomialRatio object, highest power
first, i.e., the `b` passed to `filt()`
"""
coefb(f::PolynomialRatio{:s}) = reverse(coeffs(f.b))
coefb(f::PolynomialRatio{:z}) = strip_trailing_zeros(reverse(pad_to_same_length(coeffs(f.b), coeffs(f.a))[1]))
coefb(f::FilterCoefficients) = coefb(PolynomialRatio(f))

"""
    coefa(f)

Coefficients of the denominator of a PolynomialRatio object, highest power
first, i.e., the `a` passed to `filt()`
"""
coefa(f::PolynomialRatio{:s}) = reverse(coeffs(f.a))
coefa(f::PolynomialRatio{:z}) = strip_trailing_zeros(reverse(pad_to_same_length(coeffs(f.b), coeffs(f.a))[2]))
coefa(f::FilterCoefficients) = coefa(PolynomialRatio(f))

function strip_trailing_zeros(x)
    last_nz = findlast(!iszero, x)
    if last_nz == length(x)
        return x
    elseif last_nz == nothing
        last_nz = 1
    end
    return x[1:last_nz]
end

#
# Biquad filter in transfer function form
# A separate immutable to improve efficiency of filtering using SecondOrderSections
#
"""
    Biquad(b0, b1, b2, a1, a2)

Filter representation in terms of the transfer function of a single
second-order section given by:
```math
H(s) = \\frac{\\verb!b0! s^2+\\verb!b1! s+\\verb!b2!}{s^2+\\verb!a1! s + \\verb!a2!}
```
or equivalently:
```math
H(z) = \\frac{\\verb!b0!+\\verb!b1! z^{-1}+\\verb!b2! z^{-2}}{1+\\verb!a1! z^{-1} + \\verb!a2! z^{-2}}
```
"""
struct Biquad{Domain,T<:Number} <: FilterCoefficients{Domain}
    b0::T
    b1::T
    b2::T
    a1::T
    a2::T
end
Biquad(f::FilterCoefficients{D}) where {D} = Biquad{D}(f)
Biquad(args...) = Biquad{:z}(args...)
Biquad{D}(b0::T, b1::T, b2::T, a1::T, a2::T) where {D,T} =
    Biquad{D,T}(b0, b1, b2, a1, a2)
Biquad{D}(b0::T, b1::T, b2::T, a0::T, a1::T, a2::T, g::Number=1) where {D,T} =
    (x = g*b0/a0; Biquad{D,typeof(x)}(x, g*b1/a0, g*b2/a0, a1/a0, a2/a0))

Biquad{D,T}(f::Biquad{D}) where {D,T} = Biquad{D,T}(f.b0, f.b1, f.b2, f.a1, f.a2)

Base.promote_rule(::Type{Biquad{D,T}}, ::Type{Biquad{D,S}}) where {D,T,S} = Biquad{D,promote_type(T,S)}

ZeroPoleGain{D,Z,P,K}(f::Biquad{D}) where {D,Z,P,K} = ZeroPoleGain{D,Z,P,K}(PolynomialRatio{D}(f))
ZeroPoleGain{D}(f::Biquad) where {D} = ZeroPoleGain{D}(convert(PolynomialRatio{D}, f))

function PolynomialRatio{D,T}(f::Biquad{D}) where {D,T}
    if f.b2 == zero(T) && f.a2 == zero(T)
        if f.b1 == zero(T) && f.a1 == zero(T)
            b = T[f.b0]
            a = T[one(T)]
        else
            b = T[f.b0, f.b1]
            a = T[one(T), f.a1]
        end
    else
        b = T[f.b0, f.b1, f.b2]
        a = T[one(T), f.a1, f.a2]
    end

    PolynomialRatio{D,T}(b, a)
end
PolynomialRatio{D}(f::Biquad{D,T}) where {D,T} = PolynomialRatio{D,T}(f)

function Biquad{D,T}(f::PolynomialRatio{D}) where {D,T}
    a, b = f.a, f.b
    xs = max(length(b), length(a))

    if xs > 0 && !isone(a[xs-1])
        throw(ArgumentError("leading denominator coefficient of a Biquad must be one"))
    end
    if xs == 3
        Biquad{D,T}(b[2], b[1], b[0], a[1], a[0])
    elseif xs == 2
        Biquad{D,T}(b[1], b[0], zero(T), a[0], zero(T))
    elseif xs == 1
        Biquad{D,T}(b[0], zero(T), zero(T), zero(T), zero(T))
    elseif xs == 0
        throw(ArgumentError("cannot convert an empty PolynomialRatio to Biquad"))
    else
        throw(ArgumentError("cannot convert a filter of length > 3 to Biquad"))
    end
end
Biquad{D}(f::PolynomialRatio{D,T}) where {D,T} = Biquad{D,T}(f)

Biquad{D,T}(f::ZeroPoleGain{D}) where {D,T} = Biquad{D,T}(convert(PolynomialRatio, f))
Biquad{D}(f::ZeroPoleGain{D}) where {D} = Biquad{D}(convert(PolynomialRatio, f))

*(f::Biquad{D}, g::Number) where {D} = Biquad{D}(f.b0*g, f.b1*g, f.b2*g, f.a1, f.a2)
*(g::Number, f::Biquad{D}) where {D} = Biquad{D}(f.b0*g, f.b1*g, f.b2*g, f.a1, f.a2)

#
# Second-order sections (array of biquads)
#
"""
    SecondOrderSections(biquads, gain)

Filter representation in terms of a cascade of second-order
sections and gain. `biquads` must be specified as a vector of
`Biquads`.
"""
struct SecondOrderSections{Domain,T,G} <: FilterCoefficients{Domain}
    biquads::Vector{Biquad{Domain,T}}
    g::G
end
SecondOrderSections(f::FilterCoefficients{D}) where {D} = SecondOrderSections{D}(f)
SecondOrderSections{D}(biquads::Vector{Biquad{D,T}}, g::G) where {D,T,G} =
    SecondOrderSections{D,T,G}(biquads, g)

Base.promote_rule(::Type{SecondOrderSections{D,T1,G1}}, ::Type{SecondOrderSections{D,T2,G2}}) where {D,T1,G1,T2,G2} =
    SecondOrderSections{D,promote_type(T1,T2),promote_type(G1,G2)}

SecondOrderSections{D,T,G}(f::SecondOrderSections) where {D,T,G} =
    SecondOrderSections{D,T,G}(f.biquads, f.g)
SecondOrderSections{D}(f::SecondOrderSections{D,T,G}) where {D,T,G} = SecondOrderSections{D,T,G}(f)

function ZeroPoleGain{D,Z,P,K}(f::SecondOrderSections{D}) where {D,Z,P,K}
    z = Z[]
    p = P[]
    k = f.g
    for biquad in f.biquads
        biquadzpk = ZeroPoleGain{D}(biquad)
        append!(z, biquadzpk.z)
        append!(p, biquadzpk.p)
        k *= biquadzpk.k
    end
    ZeroPoleGain{D,Z,P,K}(z, p, k)
end
ZeroPoleGain{D}(f::SecondOrderSections{D,T,G}) where {D,T,G} =
    ZeroPoleGain{D,complex(T),complex(T),G}(f)

function Biquad{D,T}(f::SecondOrderSections{D}) where {D,T}
    if length(f.biquads) != 1
        throw(ArgumentError("only a single second order section may be converted to a biquad"))
    end
    Biquad{D,T}(f.biquads[1]*f.g)
end
Biquad{D}(f::SecondOrderSections{D,T,G}) where {D,T,G} = Biquad{D,promote_type(T,G)}(f)

PolynomialRatio{D,T}(f::SecondOrderSections{D}) where {D,T} = PolynomialRatio{D,T}(ZeroPoleGain(f))
PolynomialRatio{D}(f::SecondOrderSections{D}) where {D} = PolynomialRatio{D}(ZeroPoleGain(f))

# Group each pole in p with its closest zero in z
# Remove paired poles from p and z
function groupzp(z, p)
    unpaired = fill(true, length(z))  # whether zeros have been mapped
    n = min(length(z), length(p))
    groupedz = similar(z, n)
    for i = 1:n
        closest_zero_idx = 0
        closest_zero_val = Inf
        for j = 1:length(z)
            val = abs(z[j] - p[i])
            if val < closest_zero_val
                closest_zero_idx = j
                closest_zero_val = val
            end
        end
        groupedz[i] = splice!(z, closest_zero_idx)
    end
    ret = (groupedz, p[1:n])
    splice!(p, 1:n)
    ret
end

# Sort zeros or poles lexicographically (so that poles are adjacent to
# their conjugates). Handle repeated values. Split real and complex
# values into separate vectors. Ensure that each value has a conjugate.
function split_real_complex(x::Vector{T}) where T
    # Get counts and store in a Dict
    d = Dict{T,Int}()
    for v in x
        # needs to be in normal form since 0.0 !== -0.0
        tonormal(x) = x == 0 ? abs(x) : x
        vn = complex(tonormal(real(v)), tonormal(imag(v)))
        d[vn] = get(d, vn, 0)+1
    end

    c = T[]
    r = typeof(real(zero(T)))[]
    for k in sort!(collect(keys(d)), by=x -> (real(x), imag(x)))
        if imag(k) != 0
            if !haskey(d, conj(k))
                # No match for conjugate
                return (c, r, false)
            elseif imag(k) > 0
                # Add key and its conjugate
                for n = 1:d[k]
                    push!(c, k, conj(k))
                end
            end
        else
            for n = 1:d[k]
                push!(r, k)
            end
        end
    end
    return (c, r, true)
end

# Convert a filter to second-order sections
# The returned sections are in ZPK form
function SecondOrderSections{D,T,G}(f::ZeroPoleGain{D,Z,P}) where {D,T,G,Z,P}
    z = f.z
    p = f.p
    nz = length(z)
    n = length(p)
    nz > n && throw(ArgumentError("ZeroPoleGain must not have more zeros than poles"))

    # Split real and complex poles
    (complexz, realz, matched) = split_real_complex(z)
    matched || throw(ArgumentError("complex zeros could not be matched to their conjugates"))
    (complexp, realp, matched) = split_real_complex(p)
    matched || throw(ArgumentError("complex poles could not be matched to their conjugates"))

    # Sort poles according to distance to unit circle (nearest first)
    sort!(complexp, by=x->abs(abs(x) - 1))
    sort!(realp, by=x->abs(abs(x) - 1))

    # Group complex poles with closest complex zeros
    z1, p1 = groupzp(complexz, complexp)
    # Group real poles with remaining complex zeros
    z2, p2 = groupzp(complexz, realp)
    # Group remaining complex poles with closest real zeros
    z3, p3 = groupzp(realz, complexp)
    # Group remaining real poles with closest real zeros
    z4, p4 = groupzp(realz, realp)

    # All zeros are now paired with a pole, but not all poles are
    # necessarily paired with a zero
    @assert isempty(complexz)
    @assert isempty(realz)
    groupedz = [z1; z2; z3; z4]::Vector{Z}
    groupedp = [p1; p2; p3; p4; complexp; realp]::Vector{P}
    @assert length(groupedz) == nz
    @assert length(groupedp) == n

    # Allocate memory for biquads
    biquads = Vector{Biquad{D,T}}(undef, (n >> 1)+(n & 1))

    # Build second-order sections in reverse
    # First do complete pairs
    npairs = length(groupedp) >> 1
    odd = isodd(n)
    for i = 1:npairs
        pairidx = 2*(npairs-i)
        biquads[odd+i] = convert(Biquad, ZeroPoleGain{D}(groupedz[pairidx+1:min(pairidx+2, length(groupedz))],
                                                         groupedp[pairidx+1:pairidx+2], one(T)))
    end

    if odd
        # Now do remaining pole and (maybe) zero
        biquads[1] = convert(Biquad, ZeroPoleGain{D}(groupedz[length(groupedp):end],
                                                     [groupedp[end]], one(T)))
    end

    SecondOrderSections{D,T,G}(biquads, f.k)
end
SecondOrderSections{D}(f::ZeroPoleGain{D,Z,P,K}) where {D,Z,P,K} =
    SecondOrderSections{D,promote_type(real(Z), real(P)), K}(f)


SecondOrderSections{D,T,G}(f::Biquad{D}) where {D,T,G} = SecondOrderSections{D,T,G}([f], one(G))
SecondOrderSections{D}(f::Biquad{D,T}) where {D,T} = SecondOrderSections{D,T,Int}(f)
SecondOrderSections{D}(f::FilterCoefficients{D}) where {D} = SecondOrderSections{D}(ZeroPoleGain(f))

*(f::SecondOrderSections{D}, g::Number) where {D} = SecondOrderSections{D}(f.biquads, f.g*g)
*(g::Number, f::SecondOrderSections{D}) where {D} = SecondOrderSections{D}(f.biquads, f.g*g)
*(f1::SecondOrderSections{D}, f2::SecondOrderSections{D}) where {D} =
    SecondOrderSections{D}([f1.biquads; f2.biquads], f1.g*f2.g)
*(f1::SecondOrderSections{D}, fs::SecondOrderSections{D}...) where {D} =
    SecondOrderSections{D}(vcat(f1.biquads, map(f -> f.biquads, fs)...), f1.g*prod(f.g for f in fs))

*(f1::Biquad{D}, f2::Biquad{D}) where {D} = SecondOrderSections{D}([f1, f2], 1)
*(f1::Biquad{D}, fs::Biquad{D}...) where {D} = SecondOrderSections{D}([f1, fs...], 1)
*(f1::SecondOrderSections{D}, f2::Biquad{D}) where {D} =
    SecondOrderSections{D}([f1.biquads; f2], f1.g)
*(f1::Biquad{D}, f2::SecondOrderSections{D}) where {D} =
    SecondOrderSections{D}([f1; f2.biquads], f2.g)
