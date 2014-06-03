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

    zml = exp(-im * w)
    h = polyval(filter.b, zml) / polyval(filter.a, zml)

    return h
end


function freqz(filter::Filter, w::AbstractVector)

    filter = convert(TFFilter, filter)

    h = [freqz(filter, i) for i = w]

    return h
end


function freqz(filter::Filter, hz::Union(Number, AbstractVector), fs::Integer)

    filter = convert(TFFilter, filter)

    w = hz_to_radians_per_second(hz, fs)

    h = freqz(filter, w)

    return h
end


#######################################
#
# Frequency response of an analog filter
#
#######################################

function freqs(filter::Filter, w::Number)

    filter = convert(TFFilter, filter)

    s = im * w
    h = polyval(filter.b, s) / polyval(filter.a, s)

    return h
end


function freqs(filter::Filter, w::AbstractVector)

    filter = convert(TFFilter, filter)

    h = [freqs(filter, i) for i = w]

    return h
end


function freqs(filter::Filter, hz::Union(Number, AbstractVector), fs::Integer)

    filter = convert(TFFilter, filter)

    w = hz_to_radians_per_second(hz, fs)

    h = freqs(filter, w)

    return h
end


#######################################
#
# Helper functions
#
#######################################

function hz_to_radians_per_second(hz, fs)
    return hz * ((2*pi)/fs)
end


end
