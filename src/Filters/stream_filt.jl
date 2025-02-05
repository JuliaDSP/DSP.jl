using ..DSP: _zeropad

const PFB{T} = Matrix{T}          # polyphase filter bank

abstract type FIRKernel{T} end

# Single rate FIR kernel
struct FIRStandard{T} <: FIRKernel{T}
    h::Vector{T}
    hLen::Int
end

function FIRStandard(h::Vector)
    h    = reverse(h)
    hLen = length(h)
    FIRStandard(h, hLen)
end


# Interpolator FIR kernel
mutable struct FIRInterpolator{T} <: FIRKernel{T}
    const pfb::PFB{T}
    const interpolation::Int
    const Nϕ::Int
    const tapsPerϕ::Int
    const hLen::Int
    inputDeficit::Int
    ϕIdx::Int
end

function FIRInterpolator(h::Vector, interpolation::Int)
    pfb          = taps2pfb(h, interpolation)
    tapsPerϕ, Nϕ = size(pfb)
    inputDeficit = 1
    ϕIdx         = 1
    hLen         = length(h)

    FIRInterpolator(pfb, interpolation, Nϕ, tapsPerϕ, hLen, inputDeficit, ϕIdx)
end
FIRInterpolator(h::Vector, interpolation::Integer) = FIRInterpolator(h, Int(interpolation))

# Decimator FIR kernel
mutable struct FIRDecimator{T} <: FIRKernel{T}
    const h::Vector{T}
    const hLen::Int
    const decimation::Int
    inputDeficit::Int
end

function FIRDecimator(h::Vector, decimation::Int)
    h            = reverse(h)
    hLen         = length(h)
    inputDeficit = 1
    FIRDecimator(h, hLen, decimation, inputDeficit)
end
FIRDecimator(h::Vector, decimation::Integer) = FIRDecimator(h, Int(decimation))

# Rational resampler FIR kernel
mutable struct FIRRational{T}  <: FIRKernel{T}
    const pfb::PFB{T}
    const ratio::Rational{Int}
    const Nϕ::Int
    const ϕIdxStepSize::Int
    const tapsPerϕ::Int
    const hLen::Int
    ϕIdx::Int
    inputDeficit::Int
end

function FIRRational(h::Vector, ratio::Rational{Int})
    pfb          = taps2pfb(h, numerator(ratio))
    tapsPerϕ, Nϕ = size(pfb)
    ϕIdxStepSize = mod(denominator(ratio), numerator(ratio))
    ϕIdx         = 1
    inputDeficit = 1
    hLen         = length(h)
    FIRRational(pfb, ratio, Nϕ, ϕIdxStepSize, tapsPerϕ, hLen, ϕIdx, inputDeficit)
end
FIRRational(h::Vector, ratio::Union{Integer,Rational}) = FIRRational(h, convert(Rational{Int}, ratio))

#
# Arbitrary resampler FIR kernel
#
# This kernel is different from the others in that it has two polyphase filter banks.
# The second filter bank, dpfb, is the derivative of pfb. The purpose of this is to
# allow us to compute two y values, yLower & yUpper, without having to advance the input
# index by 1. It makes the kernel simpler by not having to store extra state in the case
# when where's at the last polyphase branch and the last available input sample. By using
# a derivative filter, we can always compute the output in that scenario.
# See section 7.6.1 in [1] for a better explanation.

mutable struct FIRArbitrary{T} <: FIRKernel{T}
    const rate::Float64
    const pfb::PFB{T}
    const dpfb::PFB{T}
    const Nϕ::Int
    const tapsPerϕ::Int
    ϕAccumulator::Float64
    ϕIdx::Int
    α::Float64
    const Δ::Float64
    inputDeficit::Int
    xIdx::Int
    const hLen::Int
end

function FIRArbitrary(h::Vector, rate::Float64, Nϕ::Int)
    dh           = [diff(h); zero(eltype(h))]
    pfb          = taps2pfb(h,  Nϕ)
    dpfb         = taps2pfb(dh, Nϕ)
    tapsPerϕ     = size(pfb, 1)
    ϕAccumulator = 0.0
    ϕIdx         = 1
    α            = 0.0
    Δ            = Nϕ/rate
    inputDeficit = 1
    xIdx         = 1
    hLen         = length(h)
    FIRArbitrary(
        rate,
        pfb,
        dpfb,
        Nϕ,
        tapsPerϕ,
        ϕAccumulator,
        ϕIdx,
        α,
        Δ,
        inputDeficit,
        xIdx,
        hLen
    )
end
FIRArbitrary(h::Vector, rate::Real, Nϕ::Integer) = FIRArbitrary(h, convert(Float64, rate), convert(Int, Nϕ))

# FIRFilter - the kernel does the heavy lifting
mutable struct FIRFilter{Tk<:FIRKernel,T}
    const kernel::Tk
    const h::Vector{T}
    const historyLen::Int
    history::Vector
end

# Constructor for single-rate, decimating, interpolating, and rational resampling filters
"""
    FIRFilter(h::Vector, ratio::Union{Integer,Rational}=1)
    FIRFilter(ratio::Union{Integer,Rational}, args...)

Construct a stateful FIRFilter object from the vector of filter taps `h`.
`ratio` is an optional rational integer which specifies
the input to output sample rate relationship (e.g. `147//160` for
converting recorded audio from 48 KHz to 44.1 KHz).

Constructs `h` with `resample_filter(ratio, args...)` if it is not provided.
"""
function FIRFilter(h::Vector, resampleRatio::Union{Integer,Rational} = 1)
    interpolation = numerator(resampleRatio)
    decimation    = denominator(resampleRatio)
    historyLen    = 0

    if resampleRatio == 1                                     # single-rate
        kernel     = FIRStandard(h)
        historyLen = kernel.hLen - 1
    elseif decimation == 1                                    # interpolate
        kernel     = FIRInterpolator(h, interpolation)
        historyLen = kernel.tapsPerϕ - 1
    elseif interpolation == 1                                 # decimate
        kernel     = FIRDecimator(h, decimation)
        historyLen = kernel.hLen - 1
    else                                                      # rational
        kernel     = FIRRational(h, resampleRatio)
        historyLen = kernel.tapsPerϕ - 1
    end

    history = zeros(historyLen)

    FIRFilter(kernel, h, historyLen, history)
end

# Constructor for arbitrary resampling filter (polyphase interpolator w/ intra-phase linear interpolation)
"""
    FIRFilter(h::Vector, rate::AbstractFloat, Nϕ::Integer=32)
    FIRFilter(rate::AbstractFloat, Nϕ::Integer=32, args...)

Returns a polyphase FIRFilter object from the vector of filter taps `h`.
`rate` is a floating point number that specifies the input to output
sample-rate relationship ``\\frac{fs_{out}}{fs_{in}}``. `Nϕ` is an
optional parameter which specifies the number of *phases* created from
`h`. `Nϕ` defaults to 32.

Constructs `h` with `resample_filter(rate, Nϕ, args...)` if it is not provided.
"""
function FIRFilter(h::Vector, rate::AbstractFloat, Nϕ::Integer=32)
    rate > 0.0 || throw(DomainError(rate, "rate must be greater than 0"))
    kernel     = FIRArbitrary(h, rate, Nϕ)
    historyLen = kernel.tapsPerϕ - 1
    history    = zeros(historyLen)
    FIRFilter(kernel, h, historyLen, history)
end

# Constructor for a resampling FIR filter, where the user needs only to set the sampling rate
function FIRFilter(rate::AbstractFloat, Nϕ::Integer=32, args...)
    h = resample_filter(rate, Nϕ, args...)
    FIRFilter(h, rate, Nϕ)
end

function FIRFilter(ratio::Union{Integer,Rational}, args...)
    h = resample_filter(ratio, args...)
    FIRFilter(h, ratio)
end


#
# setphase! set's filter kernel phase index
#
function setphase!(kernel::FIRDecimator, ϕ::Real)
    ϕ >= zero(ϕ) || throw(DomainError(ϕ, "ϕ must be >= 0"))
    xThrowaway = round(Int, ϕ)
    kernel.inputDeficit += xThrowaway
    nothing
end

function setphase!(kernel::Union{FIRInterpolator, FIRRational}, ϕ::Real)
    ϕ >= zero(ϕ) || throw(DomainError(ϕ, "ϕ must be >= 0"))
    xThrowaway, ϕIdx = divrem(round(Int, ϕ * kernel.Nϕ), kernel.Nϕ)
    kernel.inputDeficit += xThrowaway
    kernel.ϕIdx = ϕIdx + 1
    nothing
end

function setphase!(kernel::FIRArbitrary, ϕ::Real)
    ϕ >= zero(ϕ) || throw(DomainError(ϕ, "ϕ must be >= 0"))
    (ϕ, xThrowaway) = modf(ϕ)
    kernel.inputDeficit += round(Int, xThrowaway)
    kernel.ϕAccumulator = ϕ * kernel.Nϕ
    kernel.ϕIdx         = 1 + floor(Int, kernel.ϕAccumulator)
    kernel.α            = modf(kernel.ϕAccumulator)[1]
    nothing
end

setphase!(self::FIRFilter, ϕ::Real) = setphase!(self.kernel, ϕ)

#
# reset! filter and its kernel to an initial state
#

# Generic case for FIRInterpolator and FIRStandard
function reset!(kernel::FIRKernel)
    kernel
end

function reset!(kernel::FIRRational)
    kernel.ϕIdx         = 1
    kernel.inputDeficit = 1
    kernel
end

function reset!(kernel::FIRDecimator)
    kernel.inputDeficit = 1
    kernel
end

function reset!(kernel::FIRArbitrary)
    kernel.ϕAccumulator = 0.0
    kernel.ϕIdx         = 1
    kernel.α            = 0.0
    kernel.inputDeficit = 1
    kernel.xIdx         = 1
    kernel
end

# For FIRFilter, set history vector to zeros of same type and required length
function reset!(self::FIRFilter)
    fill!(self.history, 0)
    reset!(self.kernel)
    self
end

#
# taps2pfb
#
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
#  In this example, the first phase, or ϕ, is [9, 5, 1].

function taps2pfb(h::Vector{T}, Nϕ::Integer) where T
    hLen     = length(h)
    tapsPerϕ = ceil(Int, hLen / Nϕ)
    pfb      = Matrix{T}(undef, tapsPerϕ, Nϕ)
    hIdx     = 1

    for rowIdx in tapsPerϕ:-1:1, colIdx in 1:Nϕ
        tap = hIdx > hLen ? zero(T) : h[hIdx]
        pfb[rowIdx, colIdx] = tap
        hIdx += 1
    end

    return pfb
end


#
# Calculates the resulting length of a multirate filtering operation, given a
#   FIRFilter{FIRRational} object and an input vector
#
# ( It's hard to explain how this works without a diagram )
#

function outputlength(inputlength::Integer, ratio::Union{Integer,Rational}, initialϕ::Integer)
    interpolation = numerator(ratio)
    decimation    = denominator(ratio)
    outLen        = ((inputlength * interpolation) - initialϕ + 1) / decimation
    ceil(Int, outLen)
end

function outputlength(::FIRStandard, inputlength::Integer)
    inputlength
end

function outputlength(kernel::FIRInterpolator, inputlength::Integer)
    outputlength(inputlength-kernel.inputDeficit+1, kernel.interpolation, kernel.ϕIdx)
end

function outputlength(kernel::FIRDecimator, inputlength::Integer)
    outputlength(inputlength-kernel.inputDeficit+1, 1//kernel.decimation, 1)
end

function outputlength(kernel::FIRRational, inputlength::Integer)
    outputlength(inputlength-kernel.inputDeficit+1, kernel.ratio, kernel.ϕIdx)
end

function outputlength(kernel::FIRArbitrary, inputlength::Integer)
    ceil(Int, (inputlength-kernel.inputDeficit+1) * kernel.rate - kernel.ϕAccumulator / kernel.Δ)
end

function outputlength(self::FIRFilter, inputlength::Integer)
    outputlength(self.kernel, inputlength)
end


#
# Calculates the input length of a multirate filtering operation,
# given the output length
# With RoundDown, inputlength returns the largest input length such that the
# actual output length will be at most the given one.
# With RoundUp, inputlength returns the shortest input length such that the
# actual output length will be at least the given one.
#

function inputlength(outputlength::Int, ratio::Union{Integer,Rational}, initialϕ::Integer, r::RoundingMode=RoundDown)
    interpolation = numerator(ratio)
    decimation    = denominator(ratio)
    d             = r == RoundUp || r == RoundFromZero ? decimation : 1
    inLen         = (outputlength * decimation + initialϕ - d) / interpolation
    round(Int, inLen, r)
end

function inputlength(::FIRStandard, outputlength::Integer, ::RoundingMode=RoundDown)
    outputlength
end

function inputlength(kernel::FIRInterpolator, outputlength::Integer, r::RoundingMode=RoundDown)
    inLen = inputlength(outputlength, kernel.interpolation, kernel.ϕIdx, r)
    inLen += kernel.inputDeficit - 1
end

function inputlength(kernel::FIRDecimator, outputlength::Integer, r::RoundingMode=RoundDown)
    inLen  = inputlength(outputlength, 1//kernel.decimation, 1, r)
    inLen += kernel.inputDeficit - 1
end

function inputlength(kernel::FIRRational, outputlength::Integer, r::RoundingMode=RoundDown)
    inLen  = inputlength(outputlength, kernel.ratio, kernel.ϕIdx, r)
    inLen += kernel.inputDeficit - 1
end

function inputlength(kernel::FIRArbitrary, outputlength::Integer, r::RoundingMode=RoundDown)
    d      = r == RoundUp || r == RoundFromZero ? 1 : 0
    inLen  = floor(Int, (outputlength - d + kernel.ϕAccumulator / kernel.Δ)/kernel.rate) + d
    inLen += kernel.inputDeficit - 1
end

function inputlength(self::FIRFilter, outputlength::Integer, r::RoundingMode=RoundDown)
    inputlength(self.kernel, outputlength, r)
end


#
# Calculates the delay caused by the FIR filter in # samples, at the input sample rate, caused by the filter process
#

timedelay(kernel::Union{FIRRational,FIRInterpolator,FIRArbitrary}) =
    (kernel.hLen - 1) / (2 * kernel.Nϕ)
timedelay(kernel::Union{FIRStandard,FIRDecimator}) = (kernel.hLen - 1) / 2
timedelay(self::FIRFilter) = timedelay(self.kernel)

#
# Single rate filtering
#

function filt!(buffer::AbstractVector{Tb}, self::FIRFilter{FIRStandard{Th}}, x::AbstractVector{Tx}) where {Tb,Th,Tx}
    kernel              = self.kernel
    history::Vector{Tx} = self.history
    bufLen              = length(buffer)
    xLen                = length(x)

    bufLen >= xLen || throw(ArgumentError("buffer length must be >= length(x)"))

    h = kernel.h
    for i = 1:min(kernel.hLen-1, xLen)
        buffer[i] = unsafe_dot(h, history, x, i)
    end
    for i = kernel.hLen:xLen
        buffer[i] = unsafe_dot(h, x, i)
    end

    self.history = shiftin!(history, x)

    return xLen
end


#
# Interpolation
#

function filt!(buffer::AbstractVector{Tb}, self::FIRFilter{FIRInterpolator{Th}}, x::AbstractVector{Tx}) where {Tb,Th,Tx}
    kernel              = self.kernel
    history::Vector{Tx} = self.history
    interpolation       = kernel.interpolation
    xLen                = length(x)
    bufLen              = length(buffer)
    bufIdx              = 0

    if xLen < kernel.inputDeficit
        self.history = shiftin!(history, x)
        kernel.inputDeficit -= xLen
        return bufIdx
    end

    inputIdx = kernel.inputDeficit
    bufLen >= outputlength(self, xLen) || throw(ArgumentError("length(buffer) must be >= interpolation * length(x)"))

    while inputIdx <= xLen
        bufIdx += 1

        if inputIdx < kernel.tapsPerϕ
            accumulator = unsafe_dot(kernel.pfb, kernel.ϕIdx, history, x, inputIdx)
        else
            accumulator = unsafe_dot(kernel.pfb, kernel.ϕIdx, x, inputIdx)
        end

        buffer[bufIdx]          = accumulator
        (kernel.ϕIdx, inputIdx) = kernel.ϕIdx == kernel.Nϕ ? (1, inputIdx+1) : (kernel.ϕIdx+1, inputIdx)
    end

    kernel.inputDeficit = 1
    self.history        = shiftin!(history, x)

    return bufIdx
end


#
# Rational resampling
#

function filt!(buffer::AbstractVector{Tb}, self::FIRFilter{FIRRational{Th}}, x::AbstractVector{Tx}) where {Tb,Th,Tx}
    kernel              = self.kernel
    history::Vector{Tx} = self.history
    xLen                = length(x)
    bufLen              = length(buffer)
    bufIdx              = 0

    if xLen < kernel.inputDeficit
        self.history = shiftin!(history, x)
        kernel.inputDeficit -= xLen
        return bufIdx
    end

    outLen = outputlength(xLen-kernel.inputDeficit+1, kernel.ratio, kernel.ϕIdx)
    bufLen >= outLen || throw(ArgumentError("buffer is too small"))

    interpolation       = numerator(kernel.ratio)
    decimation          = denominator(kernel.ratio)
    inputIdx            = kernel.inputDeficit

    while inputIdx <= xLen
        bufIdx += 1

        if inputIdx < kernel.tapsPerϕ
            accumulator = unsafe_dot(kernel.pfb, kernel.ϕIdx, history, x, inputIdx)
        else
            accumulator = unsafe_dot(kernel.pfb, kernel.ϕIdx, x, inputIdx)
        end

        buffer[bufIdx]  = accumulator
        inputIdx       += div(kernel.ϕIdx + decimation - 1, interpolation)
        ϕIdx            = kernel.ϕIdx + kernel.ϕIdxStepSize
        kernel.ϕIdx     = ϕIdx > interpolation ? ϕIdx - interpolation : ϕIdx
    end

    kernel.inputDeficit = inputIdx - xLen
    self.history        = shiftin!(history, x)

    return bufIdx
end


#
# Decimation
#

function filt!(buffer::AbstractVector{Tb}, self::FIRFilter{FIRDecimator{Th}}, x::AbstractVector{Tx}) where {Tb,Th,Tx}
    kernel              = self.kernel
    bufLen              = length(buffer)
    xLen                = length(x)
    history::Vector{Tx} = self.history
    bufIdx              = 0

    if xLen < kernel.inputDeficit
        self.history = shiftin!(history, x)
        kernel.inputDeficit -= xLen
        return bufIdx
    end

    outLen              = outputlength(self, xLen)
    inputIdx            = kernel.inputDeficit

    nbufout = fld(xLen - inputIdx, kernel.decimation) + 1
    bufLen >= nbufout || throw(ArgumentError("buffer length insufficient"))

    while inputIdx <= xLen
        bufIdx += 1

        if inputIdx < kernel.hLen
            accumulator = unsafe_dot(kernel.h, history, x, inputIdx)
        else
            accumulator = unsafe_dot(kernel.h, x, inputIdx)
        end

        buffer[bufIdx] = accumulator
        inputIdx      += kernel.decimation
    end

    kernel.inputDeficit = inputIdx - xLen
    self.history        = shiftin!(history, x)

    return bufIdx
end


#
# Arbitrary resampling
#
# Updates FIRArbitrary state. See Section 7.5.1 in [1].
#   [1] uses a phase accumulator that increments by Δ (Nϕ/rate)

function update!(kernel::FIRArbitrary)
    kernel.ϕAccumulator += kernel.Δ

    if kernel.ϕAccumulator >= kernel.Nϕ
        Δx, kernel.ϕAccumulator = divrem(kernel.ϕAccumulator, kernel.Nϕ)
        kernel.xIdx += Int(Δx)
    end

    kernel.α, foffset = modf(kernel.ϕAccumulator)
    kernel.ϕIdx = 1 + Int(foffset)
end

function filt!(
    buffer::AbstractVector{Tb},
    self::FIRFilter{FIRArbitrary{Th}},
    x::AbstractVector{Tx}
) where {Tb,Th,Tx}
    kernel              = self.kernel
    pfb                 = kernel.pfb
    dpfb                = kernel.dpfb
    xLen                = length(x)
    bufIdx              = 0
    history::Vector{Tx} = self.history

    # Do we have enough input samples to produce one or more output samples?
    if xLen < kernel.inputDeficit
        self.history = shiftin!(history, x)
        kernel.inputDeficit -= xLen
        return bufIdx
    end

    # Skip over input samples that are not needed to produce output results.
    # We do this by seting inputIdx to inputDeficit which was calculated in the previous run.
    # InputDeficit is set to 1 when instantiation the FIRArbitrary kernel, that way the first
    #   input always produces an output.
    kernel.xIdx = kernel.inputDeficit

    while kernel.xIdx <= xLen
        bufIdx += 1

        if kernel.xIdx < kernel.tapsPerϕ
            yLower = unsafe_dot(pfb,  kernel.ϕIdx, history, x, kernel.xIdx)
            yUpper = unsafe_dot(dpfb, kernel.ϕIdx, history, x, kernel.xIdx)
        else
            yLower = unsafe_dot(pfb,  kernel.ϕIdx, x, kernel.xIdx)
            yUpper = unsafe_dot(dpfb, kernel.ϕIdx, x, kernel.xIdx)
        end

        # Used to have @inbounds. Restore @inbounds if buffer length
        # can be verified prior to access.
        buffer[bufIdx] = muladd(yUpper, kernel.α, yLower)
        update!(kernel)
    end

    kernel.inputDeficit = kernel.xIdx - xLen
    self.history        = shiftin!(history, x)

    return bufIdx
end

function filt(self::FIRFilter{Tk}, x::AbstractVector) where Tk<:FIRKernel
    buffer = allocate_output(self, x)
    bufLen = length(buffer)
    samplesWritten = filt!(buffer, self, x)
    if Tk <: FIRArbitrary
        samplesWritten == bufLen || resize!(buffer, samplesWritten)
    else
        samplesWritten == bufLen || throw(AssertionError("Length of resampled output different from expectation."))
    end
    return buffer
end

function allocate_output(sf::FIRFilter{Tk}, x::AbstractVector{Tx}) where {Th,Tx,Tk<:FIRKernel{Th}}
    # In some cases when `filt(::FIRFilter{FIRArbitrary}, x)` is called
    # with certain values of `x`, `filt!(buffer, ::FIRFilter{FIRArbitrary}, x)`
    # tries to write one sample too many to the buffer and a `BoundsError`
    # is thrown. Add one extra sample to catch these exceptional cases.
    #
    # See https://github.com/JuliaDSP/DSP.jl/issues/317
    #
    # FIXME: Remove this if and when the code in
    #        `filt!(buffer, ::FIRFilter{FIRArbitrary}, x)`
    #        is updated to properly account for pathological arbitrary rates.
    outLen = outputlength(sf, length(x))
    if Tk <: FIRArbitrary
        outLen += 1
    end
    return Vector{promote_type(Th, Tx)}(undef, outLen)
end


#
# Stateless filt implementations
#

# Single-rate, decimation, interpolation, and rational resampling.
function filt(h::Vector, x::AbstractVector, ratio::Union{Integer,Rational})
    self = FIRFilter(h, ratio)
    filt(self, x)
end

# Arbitrary resampling with polyphase interpolation and two neighbor linear interpolation.
function filt(h::Vector, x::AbstractVector, rate::AbstractFloat, Nϕ::Integer=32)
    self = FIRFilter(h, rate, Nϕ)
    filt(self, x)
end

"""
    resample(x::AbstractVector, rate::Real[, coef::Vector])

Resample `x` at rational or arbitrary `rate`.
`coef` is an optional vector of FIR filter taps. If `coef`
is not provided, the taps will be computed using a Kaiser window.

Internally, `resample` uses a polyphase `FIRFilter` object,
but performs additional operations to make resampling a signal easier.
It compensates for the `FIRFilter`'s delay (ramp-up), and appends
zeros to `x`. The result is that when the input and output signals
are plotted on top of each other, they correlate very well, but one
signal will have more samples than the other.
"""
function resample(x::AbstractVector, rate::Union{Integer,Rational}, h::Vector)
    _resample!(x, rate, FIRFilter(h, rate))
end

function resample(x::AbstractVector, rate::AbstractFloat, h::Vector, Nϕ::Integer=32)
    _resample!(x, rate, FIRFilter(h, rate, Nϕ))
end

function _resample!(x::AbstractVector{T}, rate::Real, sf::FIRFilter) where T
    undelay!(sf)
    outLen  = ceil(Int, length(x) * rate)
    xPadded = _zeropad(x, inputlength(sf, outLen, RoundUp))

    buffer = allocate_output(sf, xPadded)
    samplesWritten = filt!(buffer, sf, xPadded)
    return checked_resample_output!(buffer, outLen, samplesWritten, sf)
end

function undelay!(sf::FIRFilter)
    # Get delay, in # of samples at the output rate, caused by filtering processes
    τ = timedelay(sf)

    # Use setphase! to
    #   a) adjust the input samples to skip over before producing and output (integer part of τ)
    #   b) set the ϕ index of the PFB (fractional part of τ)
    setphase!(sf, τ)
end

function checked_resample_output!(y::AbstractVector, outLen, samplesWritten, ::FIRFilter{Tk}) where Tk<:FIRKernel
    if !(Tk <: FIRArbitrary)
        samplesWritten == length(y) || throw(AssertionError("Length of resampled output different from expectation."))
    end
    samplesWritten >= outLen || throw(AssertionError("Resample output shorter than expected."))
    length(y) == outLen || resize!(y, outLen)
    return y
end

"""
    resample(x::AbstractVector, rate::Real, args::Real...)

Constructs a filter with `resample_filter` using the optional arguments `args`,
and resamples the signal `x` with it.
"""
function resample(x::AbstractVector, rate::Real, args::Real...)
    sf = FIRFilter(rate, args...)
    _resample!(x, rate, sf)
end

"""
    resample(x::AbstractArray, rate::Real, h::Vector = resample_filter(rate); dims)
    resample(x::AbstractArray, rate::AbstractFloat, h::Vector, Nϕ=32; dims)

Resample an array `x` along dimension `dims`.
If `rate` is an `AbstractFloat`, the number of phases `Nϕ`
to be used in constructing `FIRArbitrary` can be supplied
as an optional argument, which defaults to 32.
"""
function resample(x::AbstractArray, rate::Union{Integer,Rational}, h::Vector; dims)
    sf = FIRFilter(h, rate)
    return _resample!(x, rate, sf; dims)
end
function resample(x::AbstractArray, rate::AbstractFloat, h::Vector, Nϕ=32; dims)
    sf = FIRFilter(h, rate, Nϕ)
    _resample!(x, rate, sf; dims)
end

resample(x::AbstractArray, rate::Real, args::Real...; dims) =
    _resample!(x, rate, FIRFilter(rate, args...); dims)

function _resample!(x::AbstractArray{T}, rate::Real, sf::FIRFilter; dims::Int) where T
    undelay!(sf)
    size_v  = size(x, dims)
    outLen  = ceil(Int, size_v * rate)
    xPadded = Vector{T}(undef, inputlength(sf, outLen, RoundUp))
    xPadded[length(x)+1:end] .= zero(T)
    buffer  = allocate_output(sf, xPadded)
    bufLen  = length(buffer)

    mapslices(x; dims) do v::AbstractVector
        undelay!(reset!(sf))
        length(buffer) == bufLen || resize!(buffer, bufLen)
        copyto!(xPadded, v)
        samplesWritten = filt!(buffer, sf, xPadded)
        return checked_resample_output!(buffer, outLen, samplesWritten, sf)
    end
end

#
# References
#

# [1] F.J. Harris, *Multirate Signal Processing for Communication Systems*. Prentice Hall, 2004
# [2] Dick, C.; Harris, F., "Options for arbitrary resamplers in FPGA-based modulators," Signals, Systems and Computers, 2004. Conference Record of the Thirty-Eighth Asilomar Conference on , vol.1, no., pp.777,781 Vol.1, 7-10 Nov. 2004
# [3] Kim, S.C.; Plishker, W.L.; Bhattacharyya, S.S., "An efficient GPU implementation of an arbitrary resampling polyphase channelizer," Design and Architectures for Signal and Image Processing (DASIP), 2013 Conference on, vol., no., pp.231,238, 8-10 Oct. 2013
# [4] Horridge, J.P.; Frazer, Gordon J., "Accurate arbitrary resampling with exact delay for radar applications," Radar, 2008 International Conference on , vol., no., pp.123,127, 2-5 Sept. 2008
# [5] Blok, M., "Fractional delay filter design for sample rate conversion," Computer Science and Information Systems (FedCSIS), 2012 Federated Conference on , vol., no., pp.701,706, 9-12 Sept. 2012
