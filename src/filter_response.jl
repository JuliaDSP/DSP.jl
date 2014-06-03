# Filter response functions for Julia
# Created by Robert Luke (Robert.Luke@med.kuleuven.be)

module FilterResponse

export response

using Polynomial
using ..FilterDesign


# Filter frequency response for normalised frequency in radians per sample
function response(filter::Filter, is_digital::Bool, w::Number)

    filter = convert(TFFilter, filter)

    if is_digital
        zml = exp(-im * w)
        h = polyval(filter.b, zml) / polyval(filter.a, zml)
    else
        s = im * w
        h = polyval(filter.b, s) / polyval(filter.a, s)
    end

    return h
end


# Filter frequency response for frequencies in radians per sample
function response(filter::Filter, is_digital::Bool, w::AbstractVector)

    filter = convert(TFFilter, filter)

    h = [response(filter, is_digital, i) for i = w]

    return h
end


# Filter frequency response for frequency in Hz
function response(filter::Filter, is_digital::Bool, hz::Number, fs::Integer)

    filter = convert(TFFilter, filter)

    w = hz_to_radians_per_second(hz, fs)
    h = response(filter, is_digital, )

    return h
end


# Filter frequency response for frequencies in Hz
function response(filter::Filter, is_digital::Bool, hz::AbstractVector, fs::Integer)

    filter = convert(TFFilter, filter)

    h = [response(filter, is_digital, i, fs) for i = hz]

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
