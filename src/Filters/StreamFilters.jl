#==============================================================================#
#       ___ _   _ ___  ____ ____      /    ____ ____ _  _ ____ ___ ____        #
#        |   \_/  |__] |___ [__      /     |    |  | |\ | [__   |  |__/        #
#        |    |   |    |___ ___]    /      |___ |__| | \| ___]  |  |  \ .      #
#==============================================================================#

typealias PFB{T} Matrix{T}          # polyphase filter bank
typealias PNFB{T} Vector{Poly{T}}   # polynomial filter bank (used for farrow filter)

abstract Filter
abstract FIRKernel
# TODO: all kernels: add field original taps

# Single rate FIR kernel
type FIRStandard{T} <: FIRKernel
    h::Vector{T}
    hLen::Int
end

function FIRStandard( h::Vector )
    h    = flipud( h )
    hLen = length( h )
    FIRStandard( h, hLen )
end


# Interpolator FIR kernel
type FIRInterpolator{T} <: FIRKernel
    pfb::PFB{T}
    interpolation::Int
    N𝜙::Int
    tapsPer𝜙::Int
end

function FIRInterpolator( h::Vector, interpolation::Integer )
    pfb           = taps2pfb( h, interpolation )
    tapsPer𝜙      = size( pfb )[1]
    N𝜙            = size( pfb )[2]
    interpolation = interpolation
    FIRInterpolator( pfb, interpolation, N𝜙, tapsPer𝜙 )
end


# Decimator FIR kernel
type FIRDecimator{T} <: FIRKernel
    h::Vector{T}
    hLen::Int
    decimation::Int
    inputDeficit::Int
end

function FIRDecimator( h::Vector, decimation::Integer )
    h            = flipud( h )
    hLen         = length( h )
    decimation   = decimation
    inputDeficit = 1
    FIRDecimator( h, hLen, decimation, inputDeficit )
end


# Rational resampler FIR kernel
type FIRRational{T}  <: FIRKernel
    pfb::PFB{T}
    ratio::Rational{Int}
    N𝜙::Int
    tapsPer𝜙::Int
    criticalYidx::Int
    𝜙Idx::Int
    inputDeficit::Int
end

function FIRRational( h::Vector, ratio::Rational )
    pfb          = taps2pfb( h, num(ratio) )
    N𝜙           = size( pfb )[2]
    tapsPer𝜙     = size( pfb )[1]
    criticalYidx = ifloor( tapsPer𝜙 * ratio )
    𝜙Idx         = 1
    inputDeficit = 1
    FIRRational( pfb, ratio, N𝜙, tapsPer𝜙, criticalYidx, 𝜙Idx, inputDeficit )
end


# Arbitrary resampler FIR kernel
# This kernel is different from the others in that it has two polyphase filtlter banks.
# The the second filter bank, dpfb, is the derivative of pfb. The purpose of this is to
# allow us to compute two y values, yLower & yUpper, whitout having to advance the input
# index by 1. It makes the kernel simple by not having to store extra state in the case
# when where's at the last polphase branch and the last available input sample. By using
# a derivitive filter, we can always compute the output in that scenario.
# See section 7.6.1 in [1] for a better explanation.
type FIRArbitrary{T} <: FIRKernel # TODO: since farrow is also arbitrary, find a new name
    rate::Float64
    pfb::PFB{T}
    dpfb::PFB{T}
    N𝜙::Int
    tapsPer𝜙::Int
    𝜙Accumulator::Float64
    𝜙Idx::Int
    α::Float64
    Δ::Float64
    inputDeficit::Int
    xIdx::Int
end

function FIRArbitrary( h::Vector, rate::Real, N𝜙::Integer )
    dh           = [ diff( h ), 0 ]
    pfb          = taps2pfb( h,  N𝜙 )
    dpfb         = taps2pfb( dh, N𝜙 )
    tapsPer𝜙     = size( pfb )[1]
    𝜙Accumulator = 1.0
    𝜙Idx         = 1
    α            = 0.0
    Δ            = N𝜙/rate
    inputDeficit = 1
    xIdx         = 1
    FIRArbitrary( rate, pfb, dpfb, N𝜙, tapsPer𝜙, 𝜙Accumulator, 𝜙Idx, α, Δ, inputDeficit, xIdx )
end


# Farrow filter kernel.
# Takes a polyphase filterbank and converts each row of taps into a polynomial.
# That we can calculate filter tap values for any arbitrary 𝜙Idx, not just integers between 1 and N𝜙
type FIRFarrow{T} <: FIRKernel
    rate::Float64
    pfb::PFB{T}
    pnfb::PNFB{T}
    polyorder::Int
    currentTaps::Vector{T}
    N𝜙::Int
    tapsPer𝜙::Int
    𝜙Idx::Float64
    Δ::Float64
    inputDeficit::Int
    xIdx::Int
end

function FIRFarrow{T}( h::Vector{T}, rate::Real, N𝜙::Integer, polyorder::Integer )
    pfb          = taps2pfb( h,  N𝜙 )
    pnfb         = pfb2pnfb( pfb, polyorder )
    tapsPer𝜙     = size( pfb )[1]
    𝜙Idx         = 1.0
    Δ            = N𝜙/rate
    inputDeficit = 1
    xIdx         = 1
    currentTaps  = T[ polyval( pnfb[tapIdx], 𝜙Idx ) for tapIdx in 1:tapsPer𝜙 ]
    FIRFarrow( rate, pfb, pnfb, polyorder, currentTaps, N𝜙, tapsPer𝜙, 𝜙Idx, Δ, inputDeficit, xIdx )
end


# FIRFilter - the kernel does the heavy lifting
type FIRFilter{Tk<:FIRKernel} <: Filter
    kernel::Tk
    history::Vector
    historyLen::Int
end

# Constructor for single-rate, decimating, interpolating, and rational resampling filters
function FIRFilter( h::Vector, resampleRatio::Rational = 1//1 )
    interpolation = num( resampleRatio )
    decimation    = den( resampleRatio )
    historyLen    = 0

    if resampleRatio == 1                                     # single-rate
        kernel     = FIRStandard( h )
        historyLen = kernel.hLen - 1
    elseif interpolation == 1                                 # decimate
        kernel     = FIRDecimator( h, decimation )
        historyLen = kernel.hLen - 1
    elseif decimation == 1                                    # interpolate
        kernel     = FIRInterpolator( h, interpolation )
        historyLen = kernel.tapsPer𝜙 - 1
    else                                                      # rational
        kernel     = FIRRational( h, resampleRatio )
        historyLen = kernel.tapsPer𝜙 - 1
    end

    history = zeros( historyLen )

    FIRFilter( kernel, history, historyLen )
end

# Constructor for arbitrary resampling filter (polyphase interpolator w/ intra-phase linear interpolation )
function FIRFilter( h::Vector, rate::FloatingPoint, N𝜙::Integer = 32 )
    rate > 0.0 || error( "rate must be greater than 0" )
    kernel     = FIRArbitrary( h, rate, N𝜙 )
    historyLen = kernel.tapsPer𝜙 - 1
    history    = zeros( historyLen )
    FIRFilter( kernel, history, historyLen )
end

# Constructor for farrow filter (polyphase interpolator w/ polynomial genrated intra-phase taps )
function FIRFilter( h::Vector, rate::FloatingPoint, N𝜙::Integer, polyorder::Integer )
    rate > 0.0 || error( "rate must be greater than 0" )
    kernel     = FIRFarrow( h, rate, N𝜙, polyorder )
    historyLen = kernel.tapsPer𝜙 - 1
    history    = zeros( historyLen )
    FIRFilter( kernel, history, historyLen )
end



#==============================================================================#
#                    ____ ____ ___ ___  _  _ ____ ____ ____                    #
#                    [__  |___  |  |__] |__| |__| [__  |___                    #
#                    ___] |___  |  |    |  | |  | ___] |___                    #
#==============================================================================#
# Sets the kernel's phase (𝜙Idx+α).
#   Valid input is [0, 1]

function setphase( kernel::Union(FIRInterpolator, FIRRational), 𝜙::Number )
    @assert zero(𝜙) <= 𝜙 <= one(𝜙)
    kernel.𝜙Idx = int(𝜙Idx)
    return kernel.𝜙Idx
end

function setphase( kernel::FIRArbitrary, 𝜙::Number )
    @assert zero(𝜙) <= 𝜙 <= one(𝜙)
    (α, 𝜙Idx)   = modf( 𝜙 * kernel.N𝜙 )
    kernel.𝜙Idx = int(𝜙Idx)
    kernel.α    = α
    return 𝜙Idx, α
end

function setphase( kernel::FIRFarrow, 𝜙::Number )
    @assert zero(𝜙) <= 𝜙 <= one(𝜙)
    kernel.𝜙Idx = 𝜙*(kernel.N𝜙-1)+1
    tapsforphase!( kernel.currentTaps, kernel, kernel.𝜙Idx  )
    return kernel.𝜙Idx
end


setphase( self::FIRFilter, 𝜙::Number ) = setphase( self.kernel, 𝜙 )


#==============================================================================#
#                            ____ ____ ____ ____ ___                           #
#                            |__/ |___ [__  |___  |                            #
#                            |  \ |___ ___] |___  |                            #
#==============================================================================#

# Resets filter and its kernel to an initial state

# Does nothing for non-rational kernels
reset( self::FIRKernel ) = self

# For rational kernel, set 𝜙Idx back to 1
reset( self::FIRRational ) = self.𝜙Idx = 1

# For rational kernel, set 𝜙Idx back to 1
function reset( self::FIRArbitrary )
    self.yCount = 0
    update( self )
end

# For FIRFilter, set history vector to zeros of same type and required length
function reset( self::FIRFilter )
    self.history = zeros( eltype( self.history ), self.historyLen )
    reset( self.kernel )
    return self
end




#==============================================================================#
#                      _  _    ___ ____    ___  ____ ___                       #
#                      |__|     |  |  |    |__] |___ |__]                      #
#                      |  |     |  |__|    |    |    |__]                      #
#==============================================================================#

# Converts a vector of coefficients to a matrix. Each column is a filter.
# NOTE: also flips the matrix up/down so computing the dot product of a
#       column with a signal vector is more efficient (since filter is convolution)
# Appends zeros if necessary.
# Example:
#   julia> taps2pfb( [1:9], 4 )
#   3x4 Array{Int64,2}:
#    9  0  0  0
#    5  6  7  8
#    1  2  3  4
#
#  In this example, the first phase, or 𝜙, is [9, 5, 1].

function taps2pfb{T}( h::Vector{T}, N𝜙::Integer )
    hLen     = length( h )
    tapsPer𝜙 = iceil( hLen/N𝜙 )
    pfbSize  = tapsPer𝜙 * N𝜙
    pfb      = Array( T, tapsPer𝜙, N𝜙 )
    hIdx     = 1

    for rowIdx in tapsPer𝜙:-1:1, colIdx in 1:N𝜙
        tap = hIdx > hLen ? zero(T) : h[hIdx]
        @inbounds pfb[rowIdx,colIdx] = tap
        hIdx += 1
    end

    return pfb
end




#==============================================================================#
#        ___  ____ _    _   _ _  _ ____ _  _ _ ____ _       ____ ___           #
#        |__] |  | |     \_/  |\ | |  | |\/| | |__| |       |___ |__]          #
#        |    |__| |___   |   | \| |__| |  | | |  | |___    |    |__]          #
#==============================================================================#

# Convert a polyphase filterbank into a polynomial filterbank

function pfb2pnfb{T}( pfb::PFB{T}, polyorder::Integer )
    (tapsPer𝜙, N𝜙) = size( pfb )
    result         = Array( Poly{T}, tapsPer𝜙 )

    for i in 1:tapsPer𝜙
        row = vec( pfb[i,:] )
        result[i] = polyfit( row, polyorder )
    end

    return result
end

function taps2pnfb{T}( h::Vector{T}, N𝜙::Integer, polyorder::Integer )
    hLen     = length( h )
    tapsPer𝜙 = iceil( hLen/N𝜙 )
    pnfb     = Array( Poly{T}, tapsPer𝜙 )
    pfbSize  = N𝜙 * tapsPer𝜙
    h        = hLen < pfbSize + 1 ? [ h, zeros( T, pfbSize+1-hLen ) ] : h

    pnfbIdx = tapsPer𝜙
    for startIdx in 0:N𝜙:hLen-N𝜙
        row           = h[startIdx+1:startIdx+N𝜙+1]
        pnfb[pnfbIdx] = polyfit( row, polyorder )
        pnfbIdx      -= 1
    end

    return pnfb
end


#==============================================================================#
#               ____ _  _ ___ ___  _  _ ___    _    ____ _  _                  #
#               |  | |  |  |  |__] |  |  |     |    |___ |\ |                  #
#               |__| |__|  |  |    |__|  |     |___ |___ | \|                  #
#==============================================================================#

# Calculates the resulting length of a multirate filtering operation, given a
#   FIRFilter{FIRRational} object and an input vector
#
# ( It's hard to explain how this works without a diagram )

function outputlength( inputlength::Integer, ratio::Rational, initial𝜙::Integer )
    interpolation = num( ratio )
    decimation    = den( ratio )
    outLen        = (( inputlength * interpolation ) - initial𝜙 + 1 ) / decimation
    iceil(  outLen  )
end

function outputlength( kernel::FIRStandard, inputlength::Integer )
    inputlength
end

function outputlength( kernel::FIRInterpolator, inputlength::Integer )
    kernel.interpolation * inputlength
end

function outputlength( kernel::FIRDecimator, inputlength::Integer )
    outputlength( inputlength-kernel.inputDeficit+1, 1//kernel.decimation, 1 )
end

function outputlength( kernel::FIRRational, inputlength::Integer )
    outputlength( inputlength-kernel.inputDeficit+1, kernel.ratio, kernel.𝜙Idx )
end

function outputlength( kernel::FIRArbitrary, inputlength::Integer )
    iceil( (inputlength-kernel.inputDeficit+1) * kernel.rate )
end

function outputlength( kernel::FIRFarrow, inputlength::Integer )
    iceil( (inputlength-kernel.inputDeficit+1) * kernel.rate )
end

function outputlength( self::FIRFilter, inputlength::Integer )
    outputlength( self.kernel, inputlength )
end




#==============================================================================#
#                 _ _  _ ___  _  _ ___    _    ____ _  _                       #
#                 | |\ | |__] |  |  |     |    |___ |\ |                       #
#                 | | \| |    |__|  |     |___ |___ | \|                       #
#==============================================================================#

function inputlength( outputlength::Int, ratio::Rational, initial𝜙::Integer )
    interpolation = num( ratio )
    decimation    = den( ratio )
    inLen         = ( outputlength * decimation + initial𝜙 - 1 ) / interpolation
    iceil( inLen )
end

function inputlength( self::FIRFilter{FIRStandard}, outputlength::Integer )
    outputlength
end

function inputlength( self::FIRFilter{FIRInterpolator}, outputlength::Integer )
    kernel = self.kernel
    inputlength( outputlength, kernel.interpolation//1, 1 )
end

function inputlength( self::FIRFilter{FIRDecimator}, outputlength::Integer )
    kernel = self.kernel
    inLen  = inputlength( outputlength, 1//kernel.decimation, 1 )
    inLen  = inLen + kernel.inputlength - 1
end

function inputlength( self::FIRFilter{FIRRational}, outputlength::Integer )
    kernel = self.kernel
    inLen = inputlength( outputlength, kernel.ratio, kernel.𝜙Idx )
    inLen = inLen + kernel.inputDeficit - 1
end




#==============================================================================#
#              _  _ ____ _  _ ___    ___  _  _ ____ ____ ____                  #
#              |\ | |___  \/   |     |__] |__| |__| [__  |___                  #
#              | \| |___ _/\_  |     |    |  | |  | ___] |___                  #
#==============================================================================#

function nextphase( currentphase::Integer, ratio::Rational )
    interpolation = num( ratio )
    decimation    = den( ratio )
    𝜙Step         = mod( decimation, interpolation )
    𝜙Next         = currentphase + 𝜙Step
    𝜙Next         = 𝜙Next > interpolation ? 𝜙Next - interpolation : 𝜙Next
end




#==============================================================================#
#               ____ _ _  _ ____ _    ____    ____ ____ ___ ____               #
#               [__  | |\ | | __ |    |___    |__/ |__|  |  |___               #
#               ___] | | \| |__] |___ |___    |  \ |  |  |  |___               #
#==============================================================================#

function filt!{Tb,Th,Tx}( buffer::Vector{Tb}, self::FIRFilter{FIRStandard{Th}}, x::Vector{Tx} )
    kernel              = self.kernel
    history::Vector{Tx} = self.history
    hLen                = kernel.hLen
    historyLen          = self.historyLen
    bufLen              = length( buffer )
    xLen                = length( x )
    outLen              = xLen
    criticalYidx        = min( hLen, outLen )

    bufLen >= xLen || error( "buffer length must be >= x length" )

    for yIdx in 1:criticalYidx        # this first loop takes care of filter ramp up and previous history
        @inbounds buffer[yIdx] = unsafe_dot( kernel.h, history, x, yIdx )
    end

    for yIdx in criticalYidx+1:xLen
        @inbounds buffer[yIdx] = unsafe_dot( kernel.h, x, yIdx )
    end

    self.history = shiftin!( history, x )

    return buffer
end

function filt{Th,Tx}( self::FIRFilter{FIRStandard{Th}}, x::Vector{Tx} )
    buffer = Array( promote_type(Th, Tx), length(x) )
    filt!( buffer, self, x )
end




#==============================================================================#
#               _ _  _ ___ ____ ____ ___  _    ____ ____ ___ ____              #
#               | |\ |  |  |___ |__/ |__] |    |  | |__|  |  |___              #
#               | | \|  |  |___ |  \ |    |___ |__| |  |  |  |___              #
#==============================================================================#

function filt!{Tb,Th,Tx}( buffer::Vector{Tb}, self::FIRFilter{FIRInterpolator{Th}}, x::Vector{Tx} )
    kernel              = self.kernel
    history::Vector{Tx} = self.history
    interpolation       = kernel.interpolation
    N𝜙                  = kernel.N𝜙
    tapsPer𝜙            = kernel.tapsPer𝜙
    xLen                = length( x )
    bufLen              = length( buffer )
    historyLen          = self.historyLen
    outLen              = outputlength( self, xLen )
    criticalYidx        = min( historyLen*interpolation, outLen )
    inputIdx            = 1
    𝜙                   = 1

    bufLen >= outLen || error( "length( buffer ) must be >= interpolation * length(x)")

    for yIdx in 1:criticalYidx
        @inbounds buffer[yIdx] = unsafe_dot( kernel.pfb, 𝜙, history, x, inputIdx )
        (𝜙, inputIdx) = 𝜙 == N𝜙 ? ( 1, inputIdx+1 ) : ( 𝜙+1, inputIdx )
    end
    for yIdx in criticalYidx+1:outLen
        @inbounds buffer[yIdx] = unsafe_dot( kernel.pfb, 𝜙, x, inputIdx )
        (𝜙, inputIdx) = 𝜙 == N𝜙 ? ( 1, inputIdx+1 ) : ( 𝜙+1, inputIdx )
    end

    self.history = shiftin!( history, x )

    return buffer
end

function filt{Th,Tx}( self::FIRFilter{FIRInterpolator{Th}}, x::Vector{Tx} )
    xLen   = length( x )
    outlen = outputlength( self, xLen )
    buffer = Array( promote_type(Th,Tx), outlen )
    filt!( buffer, self, x )
    return buffer
end




#==============================================================================#
#           ____ ____ ___     ____ ____ ____ ____ _  _ ___  _    ____          #
#           |__/ |__|  |      |__/ |___ [__  |__| |\/| |__] |    |___          #
#           |  \ |  |  |  .   |  \ |___ ___] |  | |  | |    |___ |___          #
#==============================================================================#

function filt!{Tb,Th,Tx}( buffer::Vector{Tb}, self::FIRFilter{FIRRational{Th}}, x::Vector{Tx} )
    kernel              = self.kernel
    history::Vector{Tx} = self.history
    xLen                = length( x )
    bufLen              = length( buffer )
    bufIdx              = 0

    if xLen < kernel.inputDeficit
        self.history = shiftin!( history, x )
        kernel.inputDeficit -= xLen
        return bufIdx
    end

    outLen = outputlength( xLen-kernel.inputDeficit+1, kernel.ratio, kernel.𝜙Idx )
    bufLen >= outLen || error( "buffer is too small" )

    interpolation       = num( kernel.ratio )
    decimation          = den( kernel.ratio )
    𝜙IdxStepSize        = mod( decimation, interpolation )
    critical𝜙Idx        = kernel.N𝜙 - 𝜙IdxStepSize
    inputIdx            = kernel.inputDeficit

    while inputIdx <= xLen
        bufIdx += 1
        if inputIdx < kernel.tapsPer𝜙
            accumulator = unsafe_dot( kernel.pfb, kernel.𝜙Idx, history, x, inputIdx )
        else
            accumulator = unsafe_dot( kernel.pfb, kernel.𝜙Idx, x, inputIdx )
        end

        buffer[ bufIdx ] = accumulator
        inputIdx      += ifloor( ( kernel.𝜙Idx + decimation - 1 ) / interpolation )
        kernel.𝜙Idx    = nextphase( kernel.𝜙Idx, kernel.ratio )
    end

    kernel.inputDeficit = inputIdx - xLen
    self.history        = shiftin!( history, x )

    return bufIdx
end

function filt{Th,Tx}( self::FIRFilter{FIRRational{Th}}, x::Vector{Tx} )
    kernel         = self.kernel
    xLen           = length( x )
    bufLen         = outputlength( self, xLen )
    buffer         = Array( promote_type(Th,Tx), bufLen )
    samplesWritten = filt!( buffer, self, x )

    samplesWritten == bufLen || resize!( buffer, samplesWritten)

    return buffer
end




#==============================================================================#
#                      ___  ____ ____ _ _  _ ____ ___ ____                     #
#                      |  \ |___ |    | |\/| |__|  |  |___                     #
#                      |__/ |___ |___ | |  | |  |  |  |___                     #
#==============================================================================#

function filt!{Tb,Th,Tx}( buffer::Vector{Tb}, self::FIRFilter{FIRDecimator{Th}}, x::Vector{Tx} )
    kernel = self.kernel
    xLen   = length( x )

    if xLen < kernel.inputDeficit
        self.history = shiftin!( history, x )
        kernel.inputDeficit -= xLen
        return Tx[]
    end

    outLen              = outputlength( self, xLen )
    history::Vector{Tx} = self.history
    inputIdx            = kernel.inputDeficit
    yIdx                = 0

    while inputIdx <= xLen
        accumulator = zero( Tb )
        yIdx       += 1

        if inputIdx < kernel.hLen
            accumulator = unsafe_dot( kernel.h, history, x, inputIdx )
        else
            accumulator = unsafe_dot( kernel.h, x, inputIdx )
        end

        buffer[ yIdx ] = accumulator
        inputIdx      += kernel.decimation
    end

    kernel.inputDeficit = inputIdx - xLen
    self.history        = shiftin!( history, x )

    return yIdx
end

function filt{Th,Tx}( self::FIRFilter{FIRDecimator{Th}}, x::Vector{Tx} )
    kernel = self.kernel
    xLen   = length( x )
    Tb     = promote_type( Th, Tx)

    if xLen < kernel.inputDeficit
        history::Vector{Tx} = self.history
        self.history       = shiftin!( history, x )
        kernel.inputDeficit -= xLen
        return Tb[]
    end

    outLen = outputlength( self, xLen )
    buffer = Array( Tb, outLen )
    filt!( buffer, self, x )

    return buffer
end




#==============================================================================#
#        ____ ____ ___      ____ ____ ____ ____ _  _ ___  _    ____ ____       #
#        |__| |__/ |__]     |__/ |___ [__  |__| |\/| |__] |    |___ |__/       #
#        |  | |  \ |__] .   |  \ |___ ___] |  | |  | |    |___ |___ |  \       #
#==============================================================================#

# Updates FIRArbitrary state. See Section 7.5.1 in [1].
#   [1] uses a phase accumilator that increments by Δ (N𝜙/rate)
function update( kernel::FIRArbitrary )
    kernel.𝜙Accumulator += kernel.Δ

    if kernel.𝜙Accumulator > kernel.N𝜙
        kernel.xIdx        += ifloor( (kernel.𝜙Accumulator-1) / kernel.N𝜙 )
        kernel.𝜙Accumulator = mod( (kernel.𝜙Accumulator-1), kernel.N𝜙 ) + 1
    end

    kernel.𝜙Idx = ifloor( kernel.𝜙Accumulator )
    kernel.α    = kernel.𝜙Accumulator - kernel.𝜙Idx
end


# Generates a vector of filter taps for an arbitrary phase index.
function tapsforphase!{T}( buffer::Vector{T}, kernel::FIRArbitrary{T}, phase::Real )
    0 <= phase <= kernel.N𝜙 + 1         || error( "phase must be >= 0 and <= N𝜙+1" )
    length( buffer ) >= kernel.tapsPer𝜙 || error( "buffer is too small" )

    (α, 𝜙Idx) = modf( phase )
    𝜙Idx      = int( 𝜙Idx )

    for tapIdx in 1:kernel.tapsPer𝜙
        buffer[tapIdx] = kernel.pfb[tapIdx,𝜙Idx] + α*kernel.dpfb[tapIdx,𝜙Idx]
    end
    buffer
end

tapsforphase{T}( kernel::FIRArbitrary{T}, phase::Real ) = tapsforphase!( Array(T,kernel.tapsPer𝜙), kernel, phase )


function filt!{Tb,Th,Tx}( buffer::Vector{Tb}, self::FIRFilter{FIRArbitrary{Th}}, x::Vector{Tx} )
    kernel              = self.kernel
    pfb                 = kernel.pfb
    dpfb                = kernel.dpfb
    xLen                = length( x )
    bufIdx              = 0
    history::Vector{Tx} = self.history
    # TODO: Remove when arb and farrow filters are rock-solid.
    # db_vec_phi          = Array(Float64, bufLen)
    # db_vec_xidx         = Array(Int, bufLen)

    # Do we have enough input samples to produce one or more output samples?
    if xLen < kernel.inputDeficit
        self.history = shiftin!( history, x )
        kernel.inputDeficit -= xLen
        return bufIdx
    end

    # Skip over input samples that are not needed to produce output results.
    # We do this by seting inputIdx to inputDeficit which was calculated in the previous run.
    # InputDeficit is set to 1 when instantiation the FIRArbitrary kernel, that way the first
    #   input always produces an output.
    kernel.xIdx = kernel.inputDeficit

    while kernel.xIdx <= xLen
        # TODO: Remove when arb and farrow filters are rock-solid.
        # db_vec_xidx[bufIdx] = kernel.xIdx
        # db_vec_phi[bufIdx]  = kernel.𝜙Idx + kernel.α
        bufIdx += 1

        if kernel.xIdx < kernel.tapsPer𝜙
            yLower = unsafe_dot( pfb,  kernel.𝜙Idx, history, x, kernel.xIdx )
            yUpper = unsafe_dot( dpfb, kernel.𝜙Idx, history, x, kernel.xIdx )
        else
            yLower = unsafe_dot( pfb,  kernel.𝜙Idx, x, kernel.xIdx )
            yUpper = unsafe_dot( dpfb, kernel.𝜙Idx, x, kernel.xIdx )
        end
        buffer[bufIdx] = yLower + yUpper * kernel.α
        update( kernel )
    end

    kernel.inputDeficit = kernel.xIdx - xLen
    self.history        = shiftin!( history, x )

    # TODO: Remove when arb and farrow filters are rock-solid.
    # resize!( db_vec_phi, length(buffer) )
    # resize!( db_vec_xidx, length(buffer) )
    # return buffer, db_vec_xidx, db_vec_phi
    return bufIdx
end

function filt{Th,Tx}( self::FIRFilter{FIRArbitrary{Th}}, x::Vector{Tx} )
    bufLen         = outputlength( self, length(x) )
    buffer         = Array( promote_type(Th,Tx), bufLen )
    samplesWritten = filt!( buffer, self, x )

    samplesWritten == bufLen || resize!( buffer, samplesWritten)

    return buffer
end




#==============================================================================#
#              ____ ____ ____ ____ ____ _ _ _    ____ _ _    ___               #
#              |___ |__| |__/ |__/ |  | | | |    |___ | |     |                #
#              |    |  | |  \ |  \ |__| |_|_|    |    | |___  |                #
#==============================================================================#

# Generates a vector of filter taps for an arbitray (non-integer) phase index using polynomials
function tapsforphase!{T}( buffer::Vector{T}, kernel::FIRFarrow{T}, phase::Real )
    0 <= phase <= kernel.N𝜙 + 1         || error( "phase must be >= 0 and <= N𝜙+1" )
    length( buffer ) >= kernel.tapsPer𝜙 || error( "buffer is too small" )

    for tapIdx in 1:kernel.tapsPer𝜙
        buffer[tapIdx] = polyval( kernel.pnfb[tapIdx], phase  )
    end

    return buffer
end

tapsforphase{T}( kernel::FIRFarrow{T}, phase::Real ) = tapsforphase!( Array(T,kernel.tapsPer𝜙), kernel, phase )


# Updates farrow filter state.
# Generates new taps.
function update( kernel::FIRFarrow )
    kernel.𝜙Idx += kernel.Δ

    if kernel.𝜙Idx > kernel.N𝜙
        kernel.xIdx += ifloor( (kernel.𝜙Idx-1) / kernel.N𝜙 )
        kernel.𝜙Idx  = mod( (kernel.𝜙Idx-1), kernel.N𝜙 ) + 1
    end

    # tapsforphase!( kernel.currentTaps, kernel, kernel.𝜙Idx ) # TODO: why does this produce worse results than below?
    for tapIdx in 1:kernel.tapsPer𝜙
        @inbounds kernel.currentTaps[tapIdx] = polyval( kernel.pnfb[tapIdx], kernel.𝜙Idx )
    end
end


function filt!{Tb,Th,Tx}( buffer::Vector{Tb}, self::FIRFilter{FIRFarrow{Th}}, x::Vector{Tx} )
    kernel              = self.kernel
    xLen                = length( x )
    bufIdx              = 0
    history::Vector{Tx} = self.history
    # TODO: Remove when arb and farrow filters are rock-solid.
    # db_vec_phi          = Array(Float64, bufLen)
    # db_vec_xidx         = Array(Int, bufLen)

    # Do we have enough input samples to produce one or more output samples?
    if xLen < kernel.inputDeficit
        self.history = shiftin!( history, x )
        kernel.inputDeficit -= xLen
        return bufIdx
    end

    # Skip over input samples that are not needed to produce output results.
    kernel.xIdx = kernel.inputDeficit

    while kernel.xIdx <= xLen
        bufIdx        += 1
        # TODO: Remove when arb and farrow filters are rock-solid.
        # db_vec_xidx[bufIdx] = kernel.xIdx
        # db_vec_phi[bufIdx]  = kernel.𝜙Idx
        if kernel.xIdx < kernel.tapsPer𝜙
            y = unsafe_dot( kernel.currentTaps, history, x, kernel.xIdx )
        else
            y = unsafe_dot( kernel.currentTaps, x, kernel.xIdx )
        end
        buffer[bufIdx] = y
        update( kernel )
    end

    kernel.inputDeficit = kernel.xIdx - xLen
    self.history        = shiftin!( history, x )

    # TODO: Remove when arb and farrow filters are rock-solid.
    # resize!( db_vec_phi, length(buffer) )
    # resize!( db_vec_xidx, length(buffer) )
    # return buffer, db_vec_xidx, db_vec_phi
    return bufIdx
end

function filt{Th,Tx}( self::FIRFilter{FIRFarrow{Th}}, x::Vector{Tx} )
    bufLen         = outputlength( self, length(x) )
    buffer         = Array( promote_type(Th,Tx), bufLen )
    samplesWritten = filt!( buffer, self, x )

    samplesWritten == bufLen || resize!( buffer, samplesWritten)

    return buffer
end




#==============================================================================#
#       ____ ___ ____ ___ ____ _    ____ ____ ____    ____ _ _    ___          #
#       [__   |  |__|  |  |___ |    |___ [__  [__     |___ | |     |           #
#       ___]  |  |  |  |  |___ |___ |___ ___] ___]    |    | |___  |           #
#==============================================================================#

# Single-rate, decimation, interpolation, and rational resampling.
function filt( h::Vector, x::Vector, ratio::Rational = 1//1 )
    self = FIRFilter( h, ratio )
    filt( self, x )
end

# Arbitrary resampling with polyphase interpolation and two neighbor lnear interpolation.
function filt( h::Vector, x::Vector, rate::FloatingPoint, N𝜙::Integer = 32 )
    self = FIRFilter( h, rate, N𝜙 )
    filt( self, x )
end

# Arbitrary resampling with polyphase interpolation and polynomial generated intra-phase taps.
function filt( h::Vector, x::Vector, rate::FloatingPoint, N𝜙::Integer, polyorder::Integer )
    self = FIRFilter( h, rate, N𝜙, polyorder )
    filt( self, x )
end




#==============================================================================#
#               ____ ____ ____ ____ ____ ____ _  _ ____ ____ ____              #
#               |__/ |___ |___ |___ |__/ |___ |\ | |    |___ [__               #
#               |  \ |___ |    |___ |  \ |___ | \| |___ |___ ___]              #
#==============================================================================#

# [1] F.J. Harris, *Multirate Signal Processing for Communication Systems*. Prentice Hall, 2004
# [2] Dick, C.; Harris, F., "Options for arbitrary resamplers in FPGA-based modulators," Signals, Systems and Computers, 2004. Conference Record of the Thirty-Eighth Asilomar Conference on , vol.1, no., pp.777,781 Vol.1, 7-10 Nov. 2004
# [3] Kim, S.C.; Plishker, W.L.; Bhattacharyya, S.S., "An efficient GPU implementation of an arbitrary resampling polyphase channelizer," Design and Architectures for Signal and Image Processing (DASIP), 2013 Conference on, vol., no., pp.231,238, 8-10 Oct. 2013
# [4] Horridge, J.P.; Frazer, Gordon J., "Accurate arbitrary resampling with exact delay for radar applications," Radar, 2008 International Conference on , vol., no., pp.123,127, 2-5 Sept. 2008
# [5] Blok, M., "Fractional delay filter design for sample rate conversion," Computer Science and Information Systems (FedCSIS), 2012 Federated Conference on , vol., no., pp.701,706, 9-12 Sept. 2012
