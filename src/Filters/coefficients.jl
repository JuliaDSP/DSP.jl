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

function shiftpoly(p::LaurentPolynomial, i)
    if i > 0
        return p * LaurentPolynomial([one(eltype(p))], 1, p.var)^i
    elseif i < 0
        return p * LaurentPolynomial([one(eltype(p))], -1, p.var)^-i
    end
    return p
end

struct PolynomialRatio{Domain,T<:Number} <: FilterCoefficients{Domain}
    b::LaurentPolynomial{T}
    a::LaurentPolynomial{T}

    function PolynomialRatio{:z,Ti}(b::LaurentPolynomial, a::LaurentPolynomial) where {Ti<:Number}
        i = max(lastindex(a), lastindex(b))
        b = shiftpoly(b, -i)
        a = shiftpoly(a, -i)
        if !isone(a[0])
            if iszero(a[0])
                throw(ArgumentError("filter must have non-zero leading denominator coefficient"))
            end
            b /= a[0]
            a /= a[0]
        end
        return new{:z,Ti}(b, a)
    end
    function PolynomialRatio{:s,Ti}(b::LaurentPolynomial, a::LaurentPolynomial) where {Ti<:Number}
        if iszero(a)
            throw(ArgumentError("filter must have non-zero denominator"))
        end
        i = min(firstindex(a), firstindex(b))
        b = shiftpoly(b, -i)
        a = shiftpoly(a, -i)
        return new{:s,Ti}(b, a)
    end
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
function PolynomialRatio{:z}(b::LaurentPolynomial{T1}, a::LaurentPolynomial{T2}) where {T1,T2}
    return PolynomialRatio{:z,typeof(one(T1)/one(T2))}(b, a)
end
function PolynomialRatio{:s}(b::LaurentPolynomial{T1}, a::LaurentPolynomial{T2}) where {T1,T2}
    return PolynomialRatio{:s,promote_type(T1, T2)}(b, a)
end

# The DSP convention for Laplace domain is highest power first. The Polynomials.jl
# convention is lowest power first.
_polyprep(D::Symbol, x, T...) = LaurentPolynomial{T...}(reverse(x), D === :z ? -length(x)+1 : 0, D)
PolynomialRatio{D,T}(b::Union{Number,Vector{<:Number}}, a::Union{Number,Vector{<:Number}}) where {D,T} =
    PolynomialRatio{D,T}(_polyprep(D, b, T), _polyprep(D, a, T))
PolynomialRatio{D}(b::Union{Number,Vector{<:Number}}, a::Union{Number,Vector{<:Number}}) where {D} =
    PolynomialRatio{D}(_polyprep(D, b), _polyprep(D, a))

PolynomialRatio{D,T}(f::PolynomialRatio{D}) where {D,T} = PolynomialRatio{D,T}(f.b, f.a)
PolynomialRatio{D}(f::PolynomialRatio{D,T}) where {D,T} = PolynomialRatio{D,T}(f)

Base.promote_rule(::Type{PolynomialRatio{D,T}}, ::Type{PolynomialRatio{D,S}}) where {D,T,S} = PolynomialRatio{D,promote_type(T,S)}

function PolynomialRatio{D,T}(f::ZeroPoleGain{D}) where {D,T<:Real}
    b = convert(LaurentPolynomial{T}, real(f.k * fromroots(f.z)))
    a = convert(LaurentPolynomial{T}, real(fromroots(f.p)))
    return PolynomialRatio{D,T}(b, a)
end
PolynomialRatio{D}(f::ZeroPoleGain{D,Z,P,K}) where {D,Z,P,K} =
    PolynomialRatio{D,promote_type(real(Z),real(P),K)}(f)

ZeroPoleGain{D,Z,P,K}(f::PolynomialRatio{D}) where {D,Z,P,K} =
    ZeroPoleGain{D,Z,P,K}(ZeroPoleGain{D}(f))
function ZeroPoleGain{D}(f::PolynomialRatio{D,T}) where {D,T}
    i = -min(firstindex(f.a), firstindex(f.b), 0)
    z = roots(shiftpoly(f.b, i))
    p = roots(shiftpoly(f.a, i))
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
coefb(f::PolynomialRatio{:s}) = reverse([zeros(firstindex(f.b)); coeffs(f.b)])
coefb(f::PolynomialRatio{:z}) = reverse([coeffs(f.b); zeros(-lastindex(f.b))])
coefb(f::FilterCoefficients) = coefb(PolynomialRatio(f))

"""
    coefa(f)

Coefficients of the denominator of a PolynomialRatio object, highest power
first, i.e., the `a` passed to `filt()`
"""
coefa(f::PolynomialRatio{:s}) = reverse([zeros(firstindex(f.a)); coeffs(f.a)])
coefa(f::PolynomialRatio{:z}) = reverse([coeffs(f.a); zeros(-lastindex(f.a))])
coefa(f::FilterCoefficients) = coefa(PolynomialRatio(f))

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
    b = T[f.b0, f.b1, f.b2]
    a = T[one(T), f.a1, f.a2]
    PolynomialRatio{D,T}(b, a)
end
PolynomialRatio{D}(f::Biquad{D,T}) where {D,T} = PolynomialRatio{D,T}(f)

function Biquad{D,T}(f::PolynomialRatio{D}) where {D,T}
    a, b = f.a, f.b
    lastidx = max(lastindex(b), lastindex(a))

    if lastidx - min(firstindex(b), firstindex(a)) >= 3
        throw(ArgumentError("cannot convert a filter of length > 3 to Biquad"))
    end
    if !isone(a[lastidx])
        throw(ArgumentError("leading denominator coefficient of a Biquad must be one"))
    end
    Biquad{D,T}(b[lastidx], b[lastidx-1], b[lastidx-2], a[lastidx-1], a[lastidx-2])
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

function find_idx_of_closest(needle, haystack)
    closest_idx = 0
    closest_val = Inf
    for j = eachindex(haystack)
        val = abs(haystack[j] - needle)
        if val < closest_val
            closest_idx = j
            closest_val = val
        end
    end
    return closest_idx
end

# Convert a filter to second-order sections
# The returned sections are in ZPK form
function SecondOrderSections{D,T,G}(f::ZeroPoleGain{D,Z,P}) where {D,T,G,Z,P}
    n = length(f.p)
    length(f.z) > n && throw(ArgumentError("ZeroPoleGain must not have more zeros than poles"))

    # sort poles by distance to unit circle but complex ones first
    # TODO: for analog filters (D==:s), this should probably be the distance to the jω-axis
    p = sort(f.p, by=x->(isreal(x), abs(abs(x) - 1)))
    z = copy(f.z)

    # group poles with close zeros
    groupedp = P[]
    sizehint!(groupedp, length(p))
    groupedz = Z[]
    sizehint!(groupedz, length(z))
    while !isempty(z)
        p1 = popfirst!(p)
        push!(groupedp, p1)
        # find zero closest to current pole
        closest_zero_idx = find_idx_of_closest(p1, z)
        z1 = splice!(z, closest_zero_idx)
        push!(groupedz, z1)
        if isreal(p1)
            if !isreal(z1)
                # real pole but complex zero: find conjugate zero and a pole close to it
                # that second pole will also be real because we did the complex ones first
                closest_zero_idx′ = findfirst(==(conj(z1)), z)
                if closest_zero_idx′ === nothing
                    throw(ArgumentError("complex zeros could not be matched to their conjugates"))
                end
                closest_pole_idx′ = find_idx_of_closest(z[closest_zero_idx′], p)
                push!(groupedp, splice!(p, closest_pole_idx′))
                push!(groupedz, splice!(z, closest_zero_idx′))
            end
        else
            # complex pole: find conjugate pole and a second zero (if possible)
            closest_pole_idx′ = findfirst(==(conj(p1)), p)
            if closest_pole_idx′ === nothing
                throw(ArgumentError("complex poles could not be matched to their conjugates"))
            end
            push!(groupedp, splice!(p, closest_pole_idx′))
            if !isempty(z)
                if isreal(z1)
                    # find the second closest zero
                    closest_zero_idx′ = find_idx_of_closest(p1, z)
                    if !isreal(z[closest_zero_idx′])
                        # put back z1 and instead use the complex zero just found
                        # its conjugate will be added below
                        groupedz[end], z[closest_zero_idx′] = z[closest_zero_idx′], z1
                        z1 = groupedz[end]
                    end
                end
                if !isreal(z1)
                    # first zero is complex, so pick its conjugate as second zero
                    closest_zero_idx′ = findfirst(==(conj(z1)), z)
                    if closest_zero_idx′ === nothing
                        throw(ArgumentError("complex zeros could not be matched to their conjugates"))
                    end
                end
                push!(groupedz, splice!(z, closest_zero_idx′))
            end
        end
    end
    append!(groupedp, p)

    @assert length(groupedz) == length(f.z)
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
