# Filter response functions

#
# Frequency response of a digital filter
#

function freqz(filter::FilterCoefficients, w::Number)
    filter = convert(PolynomialRatio, filter)
    ejw = exp(-im * w)
    polyval(filter.b, ejw) ./ polyval(filter.a, ejw)
end

function freqz(filter::FilterCoefficients, w = linspace(0, π, 250))
    filter = convert(PolynomialRatio, filter)
    [freqz(filter, i) for i = w]
end

function freqz(filter::FilterCoefficients, hz::@compat(Union{Number, AbstractVector}), fs::Number)
    filter = convert(PolynomialRatio, filter)
    freqz(filter, hz_to_radians_per_second(hz, fs))
end


#
# Phase response of a digital filter
#

function phasez(filter::FilterCoefficients, w = linspace(0, π, 250))
    h = freqz(filter, w)
    unwrap(-atan2.(imag(h), real(h)))
end


#
# Impulse response of a digital filter
#

function impz(filter::FilterCoefficients, n=100)
  i = [1; zeros(n-1)]
  filt(filter, i)
end

#
# Step response of a digital filter
#

function stepz(filter::FilterCoefficients, n=100)
  i = ones(n)
  filt(filter, i)
end


#
# Frequency response of an analog filter
#

function freqs(filter::FilterCoefficients, w::Number)
    filter = convert(PolynomialRatio, filter)
    s = im * w
    polyval(filter.b, s) ./ polyval(filter.a, s)
end

function freqs(filter::FilterCoefficients, w::AbstractVector)
    filter = convert(PolynomialRatio, filter)
    [freqs(filter, i) for i = w]
end

function freqs(filter::FilterCoefficients, hz::@compat(Union{Number, AbstractVector}), fs::Number)
    filter = convert(PolynomialRatio, filter)
    freqs(filter, hz_to_radians_per_second(hz, fs))
end


#
# Helper functions
#

function hz_to_radians_per_second(hz, fs)
    hz * ((2 * pi) / fs)
end
