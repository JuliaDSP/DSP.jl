# Filter response functions for Julia
# Created by Robert Luke (Robert.Luke@med.kuleuven.be)

module FilterResponse

export freqz, freqs

using Polynomial
using ..FilterDesign


#######################################
#
# Frequency response of a digital filter
#
#######################################


function freqz(filter::Filter, w::Number)

    filter = convert(TFFilter, filter)

    ejw = exp(-im * w)
    polyval(filter.b, ejw) / polyval(filter.a, ejw)
end


function freqz(filter::Filter, w::AbstractVector)

    filter = convert(TFFilter, filter)

    [freqz(filter, i) for i = w]
end


function freqz(filter::Filter, hz::Union(Number, AbstractVector), fs::Integer)

    filter = convert(TFFilter, filter)

    freqz(filter, hz_to_radians_per_second(hz, fs))
end


#######################################
#
# Frequency response of an analog filter
#
#######################################

function freqs(filter::Filter, w::Number)

    filter = convert(TFFilter, filter)

    s = im * w
    polyval(filter.b, s) / polyval(filter.a, s)
end


function freqs(filter::Filter, w::AbstractVector)

    filter = convert(TFFilter, filter)

    [freqs(filter, i) for i = w]
end


function freqs(filter::Filter, hz::Union(Number, AbstractVector), fs::Integer)

    filter = convert(TFFilter, filter)

    freqs(filter, hz_to_radians_per_second(hz, fs))
end


#######################################
#
# Helper functions
#
#######################################

function hz_to_radians_per_second(hz, fs)
    hz * ((2 * pi) / fs)
end


end
