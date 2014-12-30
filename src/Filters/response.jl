# Filter response functions

#
# Frequency response of a digital filter
#

function freqz(filter::Filter, w::Number)
    filter = convert(TFFilter, filter)
    ejw = exp(-im * w)
    polyval(filter.b, ejw) ./ polyval(filter.a, ejw)
end

function freqz(filter::Filter, w::AbstractVector)
    filter = convert(TFFilter, filter)
    [freqz(filter, i) for i = w]
end

function freqz(filter::Filter, hz::Union(Number, AbstractVector), fs::Number)
    filter = convert(TFFilter, filter)
    freqz(filter, hz_to_radians_per_second(hz, fs))
end

function freqz(filter::Filter)
  filter = convert(TFFilter, filter)
  w = linspace(0, π, 250)
  [freqz(filter, i) for i = w]
end


#
# Phase response of a digital filter
#

function phasez(filter::Filter, w::AbstractVector)
    h = freqz(filter, w)
    unwrap(-atan2(imag(h), real(h)))
end

function phasez(filter::Filter)
  h = freqz(filter)
  unwrap(-atan2(imag(h), real(h)))
end

#
# Impulse response of a digital filter
#

function impz(filter::Filter, n=100)
  i = [1, zeros(n-1)]
  filt(filter, i)
end

#
# Step response of a digital filter
#

function stepz(filter::Filter, n=100)
  i = [0, ones(n-1)]
  filt(filter, i)
end


#
# Frequency response of an analog filter
#

function freqs(filter::Filter, w::Number)
    filter = convert(TFFilter, filter)
    s = im * w
    polyval(filter.b, s) ./ polyval(filter.a, s)
end

function freqs(filter::Filter, w::AbstractVector)
    filter = convert(TFFilter, filter)
    [freqs(filter, i) for i = w]
end

function freqs(filter::Filter, hz::Union(Number, AbstractVector), fs::Number)
    filter = convert(TFFilter, filter)
    freqs(filter, hz_to_radians_per_second(hz, fs))
end


#
# Helper functions
#

function hz_to_radians_per_second(hz, fs)
    hz * ((2 * pi) / fs)
end
