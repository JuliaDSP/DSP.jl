# Filter response functions

#
# Frequency response of a digital filter
#

function freqz(filter::FilterCoefficients, w::Number)
    filter = convert(PolynomialRatio, filter)
    ejw = exp(-im * w)
    polyval(filter.b, ejw) ./ polyval(filter.a, ejw)
end

function freqz(filter::ZeroPoleGain, w::Number)
    ejw = exp(im * w)
    filter.k * prod([ejw - z for z in filter.z]) / prod([ejw - p for p in filter.p])
end

function freqz(filter::Biquad, w::Number)
    ejw = exp(-im * w)
    ejw2 = ejw*ejw
    (filter.b0 + filter.b1*ejw + filter.b2*ejw2) / (1 + filter.a1*ejw  + filter.a2*ejw2)
end

function freqz(filter::SecondOrderSections, w::Number)
    filter.g * prod([freqz(b, w) for b in filter.biquads])
end

"""
    freqz(filter, w = range(0, stop=π, length=250))

Frequency response of a digital `filter` at normalised frequency
or frequencies `w` in radians/sample.
"""
function freqz(filter::FilterCoefficients, w = Compat.range(0, stop=π, length=250))
    [freqz(filter, i) for i = w]
end

"""
    freqz(filter, hz, fs)

Frequency response of a digital `filter` at frequency or
frequencies `hz` with sampling rate `fs`.
"""
function freqz(filter::FilterCoefficients, hz::Union{Number, AbstractVector}, fs::Number)
    freqz(filter, hz_to_radians_per_second(hz, fs))
end


"""
    phasez(filter, w = range(0, stop=π, length=250))

Phase response of a digital `filter` at normalised frequency
or frequencies `w` in radians/sample.
"""
function phasez(filter::FilterCoefficients, w = Compat.range(0, stop=π, length=250))
    h = freqz(filter, w)
    unwrap(-atan.(imag(h), real(h)); dims=ndims(h))
end


"""
    impz(filter, n=100)

Impulse response of a digital `filter` with `n` points.
"""
function impz(filter::FilterCoefficients, n=100)
  i = [1; zeros(n-1)]
  filt(filter, i)
end

"""
    stepz(filter, n=100)

Step response of a digital `filter` with `n` points.
"""
function stepz(filter::FilterCoefficients, n=100)
  i = ones(n)
  filt(filter, i)
end


"""
    freqs(filter, w)

Frequency response of an analog `filter` at normalised frequency
or frequencies `w` in radians/sample.
"""
function freqs(filter::FilterCoefficients, w::Number)
    filter = convert(PolynomialRatio, filter)
    s = im * w
    polyval(filter.b, s) ./ polyval(filter.a, s)
end

function freqs(filter::ZeroPoleGain, w::Number)
    s = im * w
    filter.k * prod([s - z for z in filter.z]) / prod([s - p for p in filter.p])
end

function freqs(filter::Biquad, w::Number)
    s = im * w
    s2 = s*s
    (filter.b0*s2 + filter.b1*s + filter.b2) / (s2 + filter.a1*s  + filter.a2)
end

function freqs(filter::SecondOrderSections, w::Number)
    filter.g * prod([freqs(b, w) for b in filter.biquads])
end

function freqs(filter::FilterCoefficients, w::AbstractVector)
    [freqs(filter, i) for i = w]
end

"""
    freqs(filter, hz, fs)

Frequency response of an analog `filter` at frequency or
frequencies `hz` with sampling rate `fs`.
"""
function freqs(filter::FilterCoefficients, hz::Union{Number, AbstractVector}, fs::Number)
    freqs(filter, hz_to_radians_per_second(hz, fs))
end


#
# Helper functions
#

function hz_to_radians_per_second(hz, fs)
    hz * ((2 * pi) / fs)
end
