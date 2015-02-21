# Filter types and conversions

abstract Filter

#
# Zero-pole gain form
#

immutable ZeroPoleGain{Z<:Number,P<:Number,K<:Number} <: Filter
    z::Vector{Z}
    p::Vector{P}
    k::K
end

#
# Transfer function form
#

immutable PolynomialRatio{T<:Number} <: Filter
    b::Poly{T}
    a::Poly{T}

    PolynomialRatio(b::Poly, a::Poly) =
        new(convert(Poly{T}, b/a[end]), convert(Poly{T}, a/a[end]))
end
PolynomialRatio{T<:Number}(b::Poly{T}, a::Poly{T}) = PolynomialRatio{T}(b, a)

# The DSP convention is highest power first. The Polynomials.jl
# convention is lowest power first.
function PolynomialRatio{T<:Number,S<:Number}(b::Union(T,Vector{T}), a::Union(S,Vector{S}))
    if findfirst(b) == 0 || findfirst(a) == 0
        error("filter must have non-zero numerator and denominator")
    end
    PolynomialRatio{promote_type(T,S)}(Poly(b[end:-1:findfirst(b)]), Poly(a[end:-1:findfirst(a)]))
end

function Base.convert(::Type{PolynomialRatio}, f::ZeroPoleGain)
    b = f.k*poly(f.z)
    a = poly(f.p)
    PolynomialRatio(Poly(real(b.a)), Poly(real(a.a)))
end

function Base.convert{T}(::Type{ZeroPoleGain}, f::PolynomialRatio{T})
    k = real(f.b[end])
    b = f.b / k
    z = convert(Vector{Complex{T}}, roots(b))
    p = convert(Vector{Complex{T}}, roots(f.a))
    ZeroPoleGain(z, p, k)
end

coefb(f::PolynomialRatio) = reverse(f.b.a)
coefa(f::PolynomialRatio) = reverse(f.a.a)

#
# Biquad filter in transfer function form
# A separate immutable to improve efficiency of filtering using SecondOrderSectionss
#

immutable Biquad{T} <: Filter
    b0::T
    b1::T
    b2::T
    a1::T
    a2::T
end
Biquad{T}(b0::T, b1::T, b2::T, a0::T, a1::T, a2::T, g::Real=1) =
    Biquad(g*b0/a0, g*b1/a0, g*b2/a0, a1/a0, a2/a0)

Base.convert(::Type{ZeroPoleGain}, f::Biquad) = convert(ZeroPoleGain, convert(PolynomialRatio, f))

function Base.convert{T}(::Type{PolynomialRatio}, f::Biquad{T})
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

    PolynomialRatio(b, a)
end

Base.convert(::Type{Biquad}, f::ZeroPoleGain) = convert(Biquad, convert(PolynomialRatio, f))

function Base.convert{T}(::Type{Biquad}, f::PolynomialRatio{T})
    a, b = f.a, f.b
    xs = max(length(b), length(a))

    if xs == 3
        Biquad(b[2], b[1], b[0], a[1], a[0])
    elseif xs == 2
        Biquad(b[1], b[0], zero(T), a[0], zero(T))
    elseif xs == 1
        Biquad(b[0], zero(T), zero(T), zero(T), zero(T))
    elseif xs == 0
        error("cannot convert an empty PolynomialRatio to Biquad")
    else
        error("cannot convert a filter of length > 3 to Biquad")
    end
end

*(f::Biquad, g::Number) = Biquad(f.b0*g, f.b1*g, f.b2*g, f.a1, f.a2)

#
# Second-order sections (array of biquads)
#

immutable SecondOrderSections{T,G} <: Filter
    biquads::Vector{Biquad{T}}
    g::G
end

realtype(x::DataType) = x
realtype{T}(::Type{Complex{T}}) = T
complextype(T::DataType) = Complex{T}
complextype{T}(::Type{Complex{T}}) = Complex{T}

function Base.convert{T}(::Type{ZeroPoleGain}, f::SecondOrderSections{T})
    t = complextype(T)
    z = t[]
    p = t[]
    k = f.g
    for biquad in f.biquads
        biquadzpk = convert(ZeroPoleGain, biquad)
        append!(z, biquadzpk.z)
        append!(p, biquadzpk.p)
        k *= biquadzpk.k
    end
    ZeroPoleGain(z, p, k)
end

Base.convert(to::Union(Type{PolynomialRatio}, Type{Biquad}), f::SecondOrderSections) =
    convert(to, convert(ZeroPoleGain, f))

# Split real and complex values in a vector into separate vectors
function split_real_complex{T<:Real}(v::Vector{Complex{T}})
    c = Complex{T}[]
    r = T[]
    for x in v
        push!(ifelse(imag(x) == 0, r, c), x)
    end
    (c, r)
end
split_real_complex{T<:Real}(v::Vector{T}) = (T[], v)

# Check that each value in a vector is followed by its complex
# conjugate
function check_conjugates(v)
    length(v) & 1 == 0 || return false
    for i = 1:2:length(v)
        v[i] == conj(v[i+1]) || return false
    end
    return true
end

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

# Convert a filter to second-order sections
# The returned sections are in ZPK form
function Base.convert{Z,P}(::Type{SecondOrderSections}, f::ZeroPoleGain{Z,P})
    z = f.z
    p = f.p
    nz = length(z)
    n = length(p)
    nz > n && error("ZeroPoleGain must not have more zeros than poles")

    # Sort poles and zeros lexicographically so that matched values are adjacent
    z = sort(z, order=Base.Order.Lexicographic)
    p = sort(p, order=Base.Order.Lexicographic)

    # Sort poles according to distance to unit circle (nearest first)
    p = sort!(p, by=x->abs(abs(x) - 1))

    # Split real and complex poles
    (complexz, realz) = split_real_complex(z)
    check_conjugates(complexz) || error("complex zeros could not be matched to their conjugates")
    (complexp, realp) = split_real_complex(p)
    check_conjugates(complexp) || error("complex poles could not be matched to their conjugates")

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
    T = promote_type(realtype(Z), realtype(P))
    biquads = Array(Biquad{T}, (n >> 1)+(n & 1))

    # Build second-order sections in reverse
    # First do complete pairs
    npairs = length(groupedp) >> 1
    odd = isodd(n)
    for i = 1:npairs
        pairidx = 2*(npairs-i)
        biquads[odd+i] = convert(Biquad, ZeroPoleGain(groupedz[pairidx+1:min(pairidx+2, length(groupedz))],
                                                         groupedp[pairidx+1:pairidx+2], one(T)))
    end

    if odd
        # Now do remaining pole and (maybe) zero
        biquads[1] = convert(Biquad, ZeroPoleGain(groupedz[length(groupedp):end],
                                                     [groupedp[end]], one(T)))
    end

    SecondOrderSections(biquads, f.k)
end

Base.convert(::Type{SecondOrderSections}, f::Filter) = convert(SecondOrderSections, convert(ZeroPoleGain, f))
