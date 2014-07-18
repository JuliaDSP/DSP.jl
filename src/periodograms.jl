# The periodogram module contains functions which compute non-parametric
# estimates of the periodogram P[s] of a signal s.  An overview of some 
# of the methods is available at:
# http://www.ee.lamar.edu/gleb/adsp/Lecture%2008%20-%20Nonparametric%20SE.pdf
module Periodograms
using ..Util
export arraysplit, nextfastfft, periodogram, welch_pgram, spectrogram, power,
       freq, time

## ARRAY SPLITTER

immutable ArraySplit{T<:AbstractVector,S,W} <: AbstractVector{Vector{S}}
    s::T
    buf::Vector{S}
    n::Int
    noverlap::Int
    window::W
    k::Int

    function ArraySplit(s, n, noverlap, nfft, window)
        # n = noverlap is a problem - the algorithm will not terminate.
        (0 <= noverlap < n) || error("noverlap must be between zero and n")
        nfft >= n || error("nfft must be >= n")
        new(s, zeros(S, nfft), n, noverlap, window, div((length(s) - n), n - noverlap)+1)
    end
end
ArraySplit(s::AbstractVector, n, noverlap, nfft, window) =
    ArraySplit{typeof(s),fftintype(eltype(s)),typeof(window)}(s, n, noverlap, nfft, window)

function Base.getindex{T,S}(x::ArraySplit{T,S,Nothing}, i::Int)
    (i >= 1 && i <= x.k) || throw(BoundsError())
    copy!(x.buf, 1, x.s, (i-1)*(x.n-x.noverlap) + 1, x.n)
end
function Base.getindex{T,S,W}(x::ArraySplit{T,S,W}, i::Int)
    (i >= 1 && i <= x.k) || throw(BoundsError())
    offset = (i-1)*(x.n-x.noverlap)
    window = x.window
    for i = 1:x.n
        @inbounds x.buf[i] = x.s[offset+i]*window[i]
    end
    x.buf
end
Base.start(x::ArraySplit) = 1
Base.next(x::ArraySplit, i::Int) = (x[i], i+1)
Base.done(x::ArraySplit, i::Int) = i > x.k
Base.size(x::ArraySplit) = (x.k,)
Base.similar(x::ArraySplit, T::Type, args...) = Array(T, args...)

# Split an array into subarrays of length N, with overlapping regions
# of length noverlap. To avoid allocation, the returned AbstractVector
# always returns the same vector (with different contents) at all
# indices.
arraysplit(s, n, noverlap, nfft=n, window=nothing) = ArraySplit(s, n, noverlap, nfft, window)

## UTILITY FUNCTIONS

# Convert the output of an FFT to a PSD and add it to out
function fft2pow!{T}(out::Array{T}, s_fft::Vector{Complex{T}}, nfft::Int, r::Real, onesided::Bool, offset::Int=0)
    m1 = convert(T, 1/r)
    n = length(s_fft)
    if onesided
        m2 = convert(T, 2/r)
        out[offset+1] += abs2(s_fft[1])*m1
        for i = 2:n-1
            @inbounds out[offset+i] += abs2(s_fft[i])*m2
        end
        out[offset+n] += abs2(s_fft[end])*ifelse(iseven(nfft), m1, m2)
    else
        if n == nfft
            for i = 1:length(s_fft)
                @inbounds out[offset+i] += abs2(s_fft[i])*m1
            end
        else
            # Convert real FFT to two-sided
            out[offset+1] += abs2(s_fft[1])*m1
            @inbounds for i = 2:length(s_fft)-1
                v = abs2(s_fft[i])*m1
                out[offset+i] += v
                out[offset+nfft-i+2] += v
            end
            out[offset+n] += abs2(s_fft[n])*m1
            if isodd(nfft)
                out[offset+nfft] += abs2(s_fft[n])*m1
            end
        end
    end
    out
end

# Calculate sum of abs2
# Remove this once we drop support for Julia 0.2
if isdefined(Base, :sumabs2)
    sumabs2(x) = Base.sumabs2(x)
else
    function sumabs2(s)
        x = zero(eltype(s))
        for i = 1:length(s)
            @inbounds x += abs2(s[i])
        end
        x
    end
end

# Evaluate a window function at n points, returning both the window
# (or nothing if no window) and the squared L2 norm of the window
compute_window(::Nothing, n::Int) = (nothing, n)
function compute_window(window::Function, n::Int)
    win = window(n)::Vector{Float64}
    norm2 = sumabs2(win)
    (win, norm2)
end
function compute_window(window::AbstractVector, n::Int)
    length(window) == n || error("length of window must match input")
    (window, sumabs2(window))
end

## PERIODOGRAMS
abstract TFR{T}
immutable Periodogram{T,F<:Union(Frequencies,Range)} <: TFR{T}
    power::Vector{T}
    freq::F
end
power(p::TFR) = p.power
freq(p::TFR) = p.freq
Base.fftshift{T,F<:Frequencies}(p::Periodogram{T,F}) =
    Periodogram(p.freq.nreal == p.freq.n ? p.power : fftshift(p.power), fftshift(p.freq))
Base.fftshift{T,F<:Range}(p::Periodogram{T,F}) = p

# Compute the periodogram of a signal S, defined as 1/N*X[s(n)]^2, where X is the
# DTFT of the signal S.
function periodogram{T<:Number}(s::AbstractVector{T}; onesided::Bool=eltype(s)<:Real,
                                nfft::Int=nextfastfft(length(s)), fs::Real=1,
                                window::Union(Function,AbstractVector,Nothing)=nothing)
    onesided && T <: Complex && error("cannot compute one-sided FFT of a complex signal")
    nfft >= length(s) || error("nfft must be >= n")

    win, norm2 = compute_window(window, length(s))
    if nfft == length(s) && win == nothing && isa(s, StridedArray)
        input = s # no need to pad
    else
        input = zeros(fftintype(T), nfft)
        if win != nothing
            for i = 1:length(s)
                @inbounds input[i] = s[i]*win[i]
            end
        else
            copy!(input, s)
        end
    end

    s_fft = T <: Real ? rfft(input) : fft(input)
    Periodogram(fft2pow!(zeros(fftabs2type(T), onesided ? (nfft >> 1)+1 : nfft),
                         s_fft, nfft, fs*norm2, onesided),
                onesided ? rfftfreq(nfft, fs) : fftfreq(nfft, fs))
end

# Compute an estimate of the power spectral density of a signal s via Welch's
# method.  The resulting periodogram has length N and is computed with an overlap
# region of length M.  The method is detailed in "The Use of Fast Fourier Transform
# for the Estimation of Power Spectra: A Method based on Time Averaging over Short,
# Modified Periodograms."  P. Welch, IEEE Transactions on Audio and Electroacoustics,
# vol AU-15, pp 70-73, 1967.
function welch_pgram{T<:Number}(s::AbstractVector{T}, n::Int=length(s)>>3, noverlap::Int=n>>1;
                                onesided::Bool=eltype(s)<:Real,
                                nfft::Int=nextfastfft(n), fs::Real=1,
                                window::Union(Function,AbstractVector,Nothing)=nothing)
    onesided && T <: Complex && error("cannot compute one-sided FFT of a complex signal")

    win, norm2 = compute_window(window, n)
    sig_split = arraysplit(s, n, noverlap, nfft, win)
    out = zeros(fftabs2type(T), onesided ? (nfft >> 1)+1 : nfft)
    r = fs*norm2*length(sig_split)

    tmp = Array(fftouttype(T), T<:Real ? (nfft >> 1)+1 : nfft)
    plan = FFTW.Plan(sig_split.buf, tmp, 1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT).plan
    for sig in sig_split
        FFTW.execute(plan, sig, tmp)
        fft2pow!(out, tmp, nfft, r, onesided)
    end

    Periodogram(out, onesided ? rfftfreq(nfft, fs) : fftfreq(nfft, fs))
end

## SPECTROGRAM

if !isdefined(Base, :FloatRange)
    typealias FloatRange{T} Range{T}
end
immutable Spectrogram{T,F<:Union(Frequencies,Range)} <: TFR{T}
    power::Matrix{T}
    freq::F
    time::FloatRange{Float64}
end
Base.fftshift{T,F<:Frequencies}(p::Spectrogram{T,F}) =
    Spectrogram(p.freq.nreal == p.freq.n ? p.power : fftshift(p.power, 1), fftshift(p.freq), p.time)
Base.fftshift{T,F<:Range}(p::Spectrogram{T,F}) = p
time(p::Spectrogram) = p.time

function spectrogram{T}(s::AbstractVector{T}, n::Int=length(s)>>3, noverlap::Int=n>>1; 
                        onesided::Bool=eltype(s)<:Real,
                        nfft::Int=nextfastfft(n), fs::Real=1,
                        window::Union(Function,AbstractVector,Nothing)=nothing)
    onesided && T <: Complex && error("cannot compute one-sided FFT of a complex signal")

    win, norm2 = compute_window(window, n)
    sig_split = arraysplit(s, n, noverlap, nfft, win)
    nout = onesided ? (nfft >> 1)+1 : nfft
    out = zeros(fftabs2type(T), nout, length(sig_split))
    tmp = Array(fftouttype(T), T<:Real ? (nfft >> 1)+1 : nfft)
    r = fs*norm2

    plan = FFTW.Plan(sig_split.buf, tmp, 1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT).plan
    offset = 0
    for sig in sig_split
        FFTW.execute(plan, sig, tmp)
        fft2pow!(out, tmp, nfft, r, onesided, offset)
        offset += nout
    end

    Spectrogram(out, onesided ? rfftfreq(nfft, fs) : fftfreq(nfft, fs),
                ((0:length(sig_split)-1)*(n-noverlap)+n/2)/fs)
end

end # end module definition
