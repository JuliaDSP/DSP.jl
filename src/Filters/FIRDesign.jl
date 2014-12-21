#==============================================================================#
#               ____ ____ _  _ ____ ___ ____ _  _ ___ ____                     #
#               |    |  | |\ | [__   |  |__| |\ |  |  [__                      #
#               |___ |__| | \| ___]  |  |  | | \|  |  ___]                     #
#==============================================================================#

@enum( FIRResponse, LOWPASS, BANDPASS, HIGHPASS, BANDSTOP )




#==============================================================================#
#          _  _ ____ _ ____ ____ ____    ___  ____ ____ _ ____ _  _            #
#          |_/  |__| | [__  |___ |__/    |  \ |___ [__  | | __ |\ |            #
#          | \_ |  | | ___] |___ |  \    |__/ |___ ___] | |__] | \|            #
#==============================================================================#

function kaiserlength( transition::Real, attenuation::Real = 60; samplerate = 1.0 )

    transition = transition./samplerate
    numtaps      = iceil(( attenuation - 7.95 )/( 2*π*2.285*transition ))

    if attenuation > 50
        β = 0.1102*( attenuation - 8.7 )
    elseif attenuation >= 21
        β = 0.5842*( attenuation - 21 )^( 0.4 ) + 0.07886*( attenuation - 21 )
    else
        β = 0.0
    end

    return numtaps, β
end




#==============================================================================#
#          ____ _ ____    ___  ____ ____ ___ ____ ___ _   _ ___  ____          #
#          |___ | |__/    |__] |__/ |  |  |  |  |  |   \_/  |__] |___          #
#          |    | |  \    |    |  \ |__|  |  |__|  |    |   |    |___          #
#==============================================================================#
# numtaps       = Desired number of filter taps
#             = If FIRType is HIGHPASS, may return  samples to make it
#               a type 1 filter.
# F           = Cutoff frequency. Real for high-pass & lowpass, Vector{Real} for
#               band-pass & band-reject
# FIRResponse = The response of the filter: LOWPASS, BANDPASS, HIGHPASS, BANDSTOP

function firprototype( numtaps::Integer, F::Union(Real, Vector{Real}); response::FIRResponse = LOWPASS )
    M = numtaps-1
    if response == LOWPASS
        prototype = [ 2*F*sinc(2*F*(n-M/2)) for n = 0:M ]
    elseif response == BANDPASS
        prototype = [ 2*(F[1]*sinc(2*F[1]*(n-M/2)) - F[2]*sinc(2*F[2]*(n-M/2))) for n = 0:M ]
    elseif response == HIGHPASS
        M = isodd( M ) ? M+1 : M
        prototype = [ sinc(n-M/2) - 2*F*sinc(2*F*(n-M/2)) for n = 0:M ]
    elseif response == BANDSTOP
        prototype = [ 2*(F[2]*sinc(2*F[2]*(n-M/2)) - F[1]*sinc(2*F[1]*(n-M/2))) for n = 0:M ]
    else
        error("Not a valid FIR_TYPE")
    end

    return prototype
end




#==============================================================================#
#                        ____ _ ____ ___  ____ ____                            #
#                        |___ | |__/ |  \ |___ [__                             #
#                        |    | |  \ |__/ |___ ___]                            #
#==============================================================================#

function firdes( numtaps::Integer, cutoff::Union(Real, Vector), windowfunction::Function; response::FIRResponse = LOWPASS, samplerate = 1.0, beta = 6.75 )

    cutoff    = cutoff ./ samplerate
    prototype = firprototype( numtaps, cutoff, response=response )
    numtaps   = length( prototype )

    if windowfunction == kaiser
        return prototype .* kaiser( numtaps, beta )
    else
        return prototype .*  windowfunction( numtaps )
    end

end

function firdes( cutoff::Union(Real, Vector{Real}), transitionwidth::Real, stopbandAttenuation::Real = 60; response::FIRResponse = LOWPASS, samplerate = 1.0 )

    ( numtaps, β ) = kaiserlength( transitionwidth, stopbandAttenuation; samplerate = samplerate )
    firdes( numtaps, cutoff, kaiser, response = response, samplerate=samplerate, beta=β )

end
