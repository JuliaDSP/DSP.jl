# Filter response functions for Julia
# Created by Robert Luke (Robert.Luke@med.kuleuven.be)

module FilterResponse

export response

using Polynomial

include("filter_design.jl")
using Reexport
@reexport using .FilterDesign


# Filter frequency response for normalised frequency in radians per sample
function response(tff::TFFilter, is_digital::Bool, w::Number)

        if is_digital
            zml = exp(-im * w)
            h = polyval(tff.b, zml) / polyval(tff.a, zml)
        else
            s = im * w 
            h = polyval(tff.b, s) / polyval(tff.a, s)
        end

    return h, w
end


# Filter frequency response for frequencies in radians per sample
function response(tff::TFFilter, is_digital::Bool, w::Array)

    h = Array(Complex, size(w))
    for i = 1:length(w)
        h[i], w[i] = response(tff, is_digital, w[i])
    end

    return h, w
end


# Filter frequency response for frequency in Hz
function response(tff::TFFilter, is_digital::Bool, hz::Number, fs::Integer)
    
    w = hz_to_radians_per_second(hz, fs)

    h, w = response(tff, is_digital, w)

    hz = radians_per_second_to_hz(w, fs)

    return h, hz
end


# Filter frequency response for frequencies in Hz
function response(tff::TFFilter, is_digital::Bool, hz::Array, fs::Integer)

    h = Array(Complex, size(hz))
    hz_return = Array(Float64, size(hz))

    for i = 1:length(hz)
        h[i], hz_return[i] = response(tff, is_digital, hz[i], fs)
    end

    return h, hz_return
end


# Filter response for an array of frequencies in Hz
function response(filter_type::Filter, is_digital::Bool, hz::Array, fs::Integer)

    return response(convert(TFFilter, filter_type), is_digital, convert(Array, hz), fs)
end


# Filter response for a range of frequencies in radians per sample
function response(filter_type::Filter, is_digital::Bool, w::Range)

    return response(filter_type, is_digital, convert(Array, w))
end


# Filter response for a range of frequencies in Hz
function response(filter_type::Filter, is_digital::Bool, hz::Range, fs::Integer)

    return response(filter_type, is_digital, convert(Array, hz), fs)
end



#######################################
#
# Helper functions
#
#######################################

function hz_to_radians_per_second(hz, fs)
    return hz * ((2*pi)/fs)
end


function radians_per_second_to_hz(w, fs)
    return w * (fs/(2*pi))
end


end
