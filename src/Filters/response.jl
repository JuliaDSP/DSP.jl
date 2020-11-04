# Filter response functions

#
# Frequency response of a digital filter
#

using ..DSP: xcorr

"""
    freqresp(filter)

Frequency response of a digital `filter` at normalized frequencies
`w = range(0, stop=π, length=250)` in radians/sample.
"""
freqresp(filter::FilterCoefficients{:z}) = freqresp(filter, range(0, stop=π, length=250))

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

"""
    freqresp(filter, hz, fs)

Frequency response of a `filter` at frequency or frequencies `hz` with sampling rate `fs`
for a digital filter or frequencies `hz/fs` for an analog filter.
"""
freqresp(filter::FilterCoefficients, hz::Union{Number, AbstractVector}, fs::Number) =
    freqresp(filter,  hz * ((2 * pi) / fs))


_freq(filter::FilterCoefficients) = x::Number -> _freq(filter, x)

_freq(filter::FilterCoefficients, x::Number) =
    _freq(convert(PolynomialRatio, filter), x)
_freq(filter::PolynomialRatio, x::Number) = filter.b(x) ./ filter.a(x)
_freq(filter::ZeroPoleGain, x::Number) =
    filter.k * prod([x - z for z in filter.z]) / prod([x - p for p in filter.p])
_freq(filter::Biquad, x::Number) =
    ((filter.b0*x + filter.b1)*x + filter.b2) / ((x + filter.a1)*x  + filter.a2)
_freq(filter::SecondOrderSections, x::Number) =
    filter.g * prod([_freq(b, x) for b in filter.biquads])


"""
    phasez(filter, w = range(0, stop=π, length=250))

Phase response of a digital `filter` at normalised frequency
or frequencies `w` in radians/sample.
"""
function phasez(filter::FilterCoefficients{:z}, w = range(0, stop=π, length=250))
    h = freqresp(filter, w)
    unwrap(angle.(h); dims=ndims(h))
end


"""
    grpdelayz(fliter, w = range(0, stop=π, length=250))

Group delay of a digital 'filter' at normalized frequency
or frequencies 'w' in radians/sample.
"""
function grpdelayz(filter::FilterCoefficients{:z}, w = range(0, stop=π, length=250))
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


"""
    impz(filter, n=100)

Impulse response of a digital `filter` with `n` points.
"""
function impz(filter::FilterCoefficients{:z}, n=100)
  i = [1; zeros(n-1)]
  filt(filter, i)
end

"""
    stepz(filter, n=100)

Step response of a digital `filter` with `n` points.
"""
function stepz(filter::FilterCoefficients{:z}, n=100)
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
