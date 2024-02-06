# Filter response functions

#
# Frequency response of a digital filter
#

using ..DSP: xcorr

"""
    H, w = freqresp(filter)

Frequency response `H` of a `filter` at (normalized) frequencies `w` in
radians/sample for a digital filter or radians/second for an analog filter
chosen as a reasonable default.
"""
function freqresp(filter::FilterCoefficients)
    w = _freqrange(filter)
    return (freqresp(filter, w), w)
end

"""
    freqresp(filter::FilterCoefficients{:z}, w)

Frequency response of digital `filter` at normalized frequency or
frequencies `w` in radians/sample.
"""
freqresp(filter::FilterCoefficients{:z}, w) = _freq(filter).(exp.(im .* w))

"""
    freqresp(filter::FilterCoefficients{:s}, w)

Frequency response of analog `filter` at frequency or frequencies `w` in
radians/second.
"""
freqresp(filter::FilterCoefficients{:s}, w) = _freq(filter).(im .* w)


_freq(filter::FilterCoefficients) = x::Number -> _freq(filter, x)

_freq(filter::FilterCoefficients, x::Number) =
    _freq(convert(PolynomialRatio, filter), x)
_freq(filter::PolynomialRatio, x::Number) = filter.b(x) / filter.a(x)

_prod_freq(f, v::Vector{<:Union{Z,Biquad{D,Z}}}, ::Type{T}) where {D,T<:Number,Z<:Number} =
    mapreduce(f, Base.mul_prod, v; init=one(promote_type(T, Z)))

_freq(filter::ZeroPoleGain, x::T) where {T<:Number} =
    filter.k * _prod_freq(z -> x - z, filter.z, T) / _prod_freq(p -> x - p, filter.p, T)
_freq(filter::Biquad, x::Number) =
    muladd(muladd(filter.b0, x, filter.b1), x, filter.b2) / muladd((x + filter.a1), x, filter.a2)
_freq(filter::SecondOrderSections, x::T) where {T<:Number} =
    filter.g * _prod_freq(b -> _freq(b, x), filter.biquads, T)


"""
    phi, w = phaseresp(filter)

Phase response `phi` of a `filter` at (normalized) frequencies `w` in
radians/sample for a digital filter or radians/second for an analog filter
chosen as a reasonable default.
"""
function phaseresp(filter::FilterCoefficients)
    w = _freqrange(filter)
    return (phaseresp(filter, w), w)
end

"""
    phaseresp(filter, w)

Phase response of a `filter` at (normalized) frequency or frequencies `w` in
radians/sample for a digital filter or radians/second for an analog filter.
"""
function phaseresp(filter::FilterCoefficients, w)
    h = freqresp(filter, w)
    unwrap(angle.(h); dims=ndims(h))
end

"""
    tau, w = grpdelay(filter)

Group delay `tau` of a `filter` at (normalized) frequencies `w` in
radians/sample for a digital filter or radians/second for an analog filter
chosen as a reasonable default.
"""
function grpdelay(filter::FilterCoefficients)
    w = _freqrange(filter)
    return (grpdelay(filter, w), w)
end

"""
    grpdelay(fliter, w)

Group delay of a digital 'filter' at normalized frequency
or frequencies 'w' in radians/sample.
"""
function grpdelay(filter::FilterCoefficients{:z}, w)
    filter = convert(PolynomialRatio, filter)
    b, a = coefb(filter), coefa(filter)

    # Linear Phase FIR
    if (length(a) == 1) & (_is_sym(b) | _is_anti_sym(b))
        return fill((length(b)-1)/2, length(w))
    end

    c = xcorr(b, a; padmode = :none)
    cr = range(0, stop=length(c)-1) .* c
    ejw = exp.(-im .* w)
    num = Polynomial(cr).(ejw)
    den = Polynomial(c).(ejw)
    return real.(num ./ den) .- (length(a) - 1)
end

function grpdelay(filter::FilterCoefficients{:s}, w)
    filter = convert(PolynomialRatio, filter)
    b, a = filter.b, filter.a
    bd = derivative(b)
    ad = derivative(a)
    s = im .* w
    return real.((bd*a - ad*b).(s) ./ (a * b).(s))
end

"""
    impresp(filter, n=100)

Impulse response of a digital `filter` with `n` points.
"""
function impresp(filter::FilterCoefficients{:z}, n=100)
  i = [1; zeros(n-1)]
  filt(filter, i)
end

"""
    stepresp(filter, n=100)

Step response of a digital `filter` with `n` points.
"""
function stepresp(filter::FilterCoefficients{:z}, n=100)
  i = ones(n)
  filt(filter, i)
end



#
# Helper functions
#

function _is_sym(x::AbstractArray)
    n = length(x) ÷ 2
    return all(x[1+i] == x[end-i] for i in 0:n-1)
end

function _is_anti_sym(x::AbstractArray)
    n = length(x) ÷ 2
    return all(x[1+i] == -x[end-i] for i in 0:n)
end

_freqrange(::FilterCoefficients{:z}) = range(0, stop=π, length=257)
function _freqrange(filter::FilterCoefficients{:s})
    filter = convert(ZeroPoleGain, filter)
    w_interesting = sort!(Float64.(abs.([filter.p; filter.z])))
    include_zero = !isempty(w_interesting) && iszero(w_interesting[1])
    w_interesting = collect(Compat.Iterators.dropwhile(iszero, w_interesting))
    if isempty(w_interesting) # no non-zero poles or zeros
        if !include_zero || !isfinite(1/filter.k)
            return [0; 10 .^ (0:6)] # fallback
        end
        # include the point where |H|=1 (if any) and go further by factor 10
        return range(0.0, stop=10 * Float64(max(filter.k, 1/filter.k)), length=200)
    end
    # normal case: go from smalles to largest pole/zero, extended by factor 10
    w_min, w_max = w_interesting[[1,end]]
    w = 10 .^ range(log10(w_min)-1, stop=log10(w_max)+1, length=200)
    return include_zero ? [0.0; w] : w
end
