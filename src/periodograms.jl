# The periodogram module contains functions which compute non-parametric
# estimates of the periodogram P[s] of a signal s.
module Periodograms
using LinearAlgebra: mul!
using ..Util, ..Windows
using Statistics: mean!
export arraysplit, nextfastfft, periodogram, welch_pgram,
       spectrogram, power, freq, stft,
       MTConfig, mt_pgram, mt_pgram!,
       MTSpectrogramConfig, mt_spectrogram, mt_spectrogram!,
       MTCrossSpectraConfig, mt_cross_power_spectra, mt_cross_power_spectra!,
       MTCoherenceConfig, mt_coherence, mt_coherence!,
       coherence
import ..DSP: allocate_output
using FFTW
import FFTW: Frequencies, fftfreq, rfftfreq
## ARRAY SPLITTER

struct ArraySplit{T<:AbstractVector,S,W} <: AbstractVector{Vector{S}}
    s::T
    buf::Vector{S}
    n::Int
    noverlap::Int
    window::W
    k::Int

    function ArraySplit{Ti,Si,Wi}(s, n, noverlap, nfft, window) where {Ti<:AbstractVector,Si,Wi}
        # n = noverlap is a problem - the algorithm will not terminate.
        (0 ≤ noverlap < n) || throw(DomainError((noverlap=noverlap, n=n), "noverlap must be between zero and n"))
        nfft >= n || throw(DomainError((nfft=nfft, n=n), "nfft must be >= n"))
        new{Ti,Si,Wi}(s, zeros(Si, nfft), n, noverlap, window, length(s) >= n ? div((length(s) - n), n - noverlap)+1 : 0)
    end
end
ArraySplit(s::AbstractVector, n, noverlap, nfft, window) =
    ArraySplit{typeof(s),fftintype(eltype(s)),typeof(window)}(s, n, noverlap, nfft, window)

function Base.getindex(x::ArraySplit{T,S,Nothing}, i::Int) where {T,S}
    (1 <= i <= x.k) || throw(BoundsError(x, i))
    copyto!(x.buf, 1, x.s, (i-1)*(x.n-x.noverlap) + firstindex(x.s), x.n)
end
function Base.getindex(x::ArraySplit{T,S,W}, i::Int) where {T,S,W}
    (1 <= i <= x.k) || throw(BoundsError(x, i))
    offset = (i-1)*(x.n-x.noverlap) + firstindex(x.s) - 1
    window = x.window
    for i = 1:x.n
        @inbounds x.buf[i] = x.s[offset+i]*window[i]
    end
    x.buf
end

function Base.iterate(x::ArraySplit, i::Int = 1)
    i > x.k ? nothing : (x[i], i+1)
end
Base.size(x::ArraySplit) = (x.k,)

"""
    arraysplit(s, n, m)

Split an array into arrays of length `n` with overlapping regions
of length `m`. Iterating or indexing the returned AbstractVector
always yields the same Vector with different contents.
"""
arraysplit(s, n, noverlap, nfft=n, window=nothing) = ArraySplit(s, n, noverlap, nfft, window)

## Make collect() return the correct split arrays rather than repeats of the last computed copy
Base.collect(x::ArraySplit) = collect(copy(a) for a in x)

## UTILITY FUNCTIONS

# Convert the output of an FFT to a PSD and add it to out
function fft2pow!(out::AbstractArray{T}, s_fft::AbstractVector{Complex{T}}, nfft::Int, r::Real, onesided::Bool, offset::Int=0) where T
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
            for i in eachindex(s_fft)
                @inbounds out[offset+i] += abs2(s_fft[i])*m1
            end
        else
            # Convert real FFT to two-sided
            out[offset+1] += abs2(s_fft[1])*m1
            @inbounds for i = 2:n-1
                v = abs2(s_fft[i])*m1
                out[offset+i] += v
                out[offset+nfft-i+2] += v
            end
            out[offset+n] += abs2(s_fft[n])*m1
            if isodd(nfft)
                out[offset+n+1] += abs2(s_fft[n])*m1
            end
        end
    end
    out
end

# Convert the output of a 2-d FFT to a 2-d PSD and add it to out
function fft2pow2!(out::Matrix{T}, s_fft::Matrix{Complex{T}}, n1::Int, n2::Int, r::Real) where T
    m1 = convert(T, 1/r)
    for j = 1:n2, i = 1:n1
        @inbounds out[i,j] += abs2(s_fft[i,j])*m1
    end
    out
end
# Convert the output of a 2-d FFT to a radial PSD and add it to out
function fft2pow2radial!(out::Array{T}, s_fft::Matrix{Complex{T}}, n1::Int, n2::Int, r::Real, ptype::Int) where T
    nmin = min(n1,n2)
    n1max = n1>>1 + 1    # since rfft is used
    n1max != size(s_fft,1) && error("fft size incorrect")
    m1 = convert(T, 1/r)
    m2 = convert(T, 2/r)
    wavenum = 0          # wavenumber index
    kmax = length(out)   # the highest wavenumber
    wc = zeros(Int,kmax) # wave count for radial average
    if n1 == nmin        # scale the wavevector for non-square s_fft
        c1 = 1
        c2 = n1/n2
    else
        c1 = n2/n1
        c2 = 1
    end

    @inbounds begin
        for j = 1:n2
            kj2 = ifelse(j <= n2>>1 + 1, j-1, -n2+j-1)
            kj2 = (kj2*c2)^2

            wavenum = round(Int, sqrt( (c1*(1-1))^2 + kj2 )) + 1
            if wavenum<=kmax
                out[wavenum] += abs2(s_fft[1,j])*m1
                wc[wavenum] += 1
            end
            for i = 2:n1max-1
                wavenum = round(Int, sqrt( (c1*(i-1))^2 + kj2 )) + 1
                if wavenum<=kmax
                    out[wavenum] += abs2(s_fft[i,j])*m2
                    wc[wavenum] += 2
                end
            end
            wavenum = round(Int, sqrt( (c1*(n1max-1))^2 + kj2 )) + 1
            if wavenum<=kmax
                out[wavenum] += abs2(s_fft[n1max,j])*ifelse(iseven(n1), m1, m2)
                wc[wavenum] += ifelse(iseven(n1), 1, 2)
            end
        end
    end
    if ptype == 2
        for i = 1:kmax
            @inbounds out[i] /= wc[i]
        end
    end
    out
end

function fft2oneortwosided!(out::Array{Complex{T}}, s_fft::Vector{Complex{T}}, nfft::Int, onesided::Bool, offset::Int=0) where T
    n = length(s_fft)
    copyto!(out, offset+1, s_fft, 1, n)
    if !onesided && n != nfft
        # Convert real FFT to two-sided
        @inbounds for i = 2:n-1
            out[offset+nfft-i+2] = conj(s_fft[i])
        end
        if isodd(nfft)
            out[offset+n+1] = conj(s_fft[n])
        end
    end
    out
end

# Evaluate a window function at n points, returning both the window
# (or nothing if no window) and the squared L2 norm of the window
compute_window(::Nothing, n::Int) = (nothing, n)
function compute_window(window::Function, n::Int)
    win = window(n)::Vector{Float64}
    norm2 = sum(abs2, win)
    (win, norm2)
end
function compute_window(window::AbstractVector, n::Int)
    length(window) == n || error("length of window must match input")
    (window, sum(abs2, window))
end

## PERIODOGRAMS
abstract type TFR{T} end
struct Periodogram{T,F<:Union{Frequencies,AbstractRange}, V <: AbstractVector{T}} <: TFR{T}
    power::V
    freq::F
end
struct Periodogram2{T,F1<:Union{Frequencies,AbstractRange},F2<:Union{Frequencies,AbstractRange}, M<:AbstractMatrix{T}} <: TFR{T}
    power::M
    freq1::F1
    freq2::F2
end

"""
    power(p)

For a `Periodogram`, returns the computed power at each frequency as
a Vector.

For a `Spectrogram`, returns the computed power at each frequency and
time bin as a Matrix. Dimensions are frequency × time.

For a `CrossPowerSpectra`, returns the pairwise power between each pair
of channels at each frequency. Dimensions are channel x channel x frequency.
"""
power(p::TFR) = p.power

"""
    freq(p)

Returns the frequency bin centers for a given `Periodogram`,
`Spectrogram`, `CrossPowerSpectra`, or `Coherence` object.

Returns a tuple of frequency bin centers for a given `Periodogram2`
object.
"""
freq(p::TFR) = p.freq
freq(p::Periodogram2) = (p.freq1, p.freq2)
FFTW.fftshift(p::Periodogram{T,F}) where {T,F<:Frequencies} =
    Periodogram(p.freq.n_nonnegative == p.freq.n ? p.power : fftshift(p.power), fftshift(p.freq))
FFTW.fftshift(p::Periodogram{T,F}) where {T,F<:AbstractRange} = p
# 2-d
FFTW.fftshift(p::Periodogram2{T,F1,F2}) where {T,F1<:Frequencies,F2<:Frequencies} =
    Periodogram2(p.freq1.n_nonnegative == p.freq1.n ? fftshift(p.power,2) : fftshift(p.power), fftshift(p.freq1), fftshift(p.freq2))
FFTW.fftshift(p::Periodogram2{T,F1,F2}) where {T,F1<:AbstractRange,F2<:Frequencies} =
    Periodogram2(fftshift(p.power,2), p.freq1, fftshift(p.freq2))
FFTW.fftshift(p::Periodogram2{T,F1,F2}) where {T,F1<:AbstractRange,F2<:AbstractRange} = p

# Compute the periodogram of a signal S, defined as 1/N*X[s(n)]^2, where X is the
# DTFT of the signal S.
"""
    periodogram(s; onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, window=nothing)

Computes periodogram of a signal by FFT and returns a
Periodogram object.

For real signals, the two-sided periodogram is symmetric and this
function returns a one-sided (real only) periodogram by default. A
two-sided periodogram can be obtained by setting `onesided=false`.

`nfft` specifies the number of points to use for the Fourier
transform. If `length(s)` < `nfft`, then the input is padded
with zeros. By default, `nfft` is the closest size for which the
Fourier transform can be computed with maximal efficiency.

`fs` is the sample rate of the original signal, and `window` is
an optional window function or vector to be applied to the original
signal before computing the Fourier transform. The computed
periodogram is normalized so that the area under the periodogram is
equal to the uncentered variance (or average power) of the original
signal.
"""
function periodogram(s::AbstractVector{T}; onesided::Bool=eltype(s)<:Real,
                     nfft::Int=nextfastfft(length(s)), fs::Real=1,
                     window::Union{Function,AbstractVector,Nothing}=nothing) where T<:Number
    onesided && T <: Complex && throw(ArgumentError("cannot compute one-sided FFT of a complex signal"))
    nfft >= length(s) || throw(DomainError((nfft=nfft, n=length(s)), "nfft must be >= n = length(s)"))

    win, norm2 = compute_window(window, length(s))
    if nfft == length(s) && win === nothing && isa(s, StridedArray)
        input = s # no need to pad
    else
        input = zeros(fftintype(T), nfft)
        if win !== nothing
            for i in eachindex(s, win)
                @inbounds input[i] = s[i]*win[i]
            end
        else
            copyto!(input, s)
        end
    end

    s_fft = T <: Real ? rfft(input) : fft(input)
    Periodogram(fft2pow!(zeros(fftabs2type(T), onesided ? (nfft >> 1)+1 : nfft),
                         s_fft, nfft, fs*norm2, onesided),
                onesided ? rfftfreq(nfft, fs) : fftfreq(nfft, fs))
end

# Compute the periodogram of a 2-d signal S. Returns 1/N*X[s(n)]^2, where X is the 2-d
# DTFT of the signal S if radialsum and radialavg are both false (default),
# a radial sum if radialsum=true, or a radial averave if radialavg=true
"""
    periodogram(s::AbstractMatrix; nfft=nextfastfft(size(s)), fs=1, radialsum=false, radialavg=false)

Computes periodogram of a 2-d signal by FFT and returns a
Periodogram2 object.

Returns a 2-d periodogram by default. A radially summed or
averaged periodogram is returned as a Periodogram object
if `radialsum` or  `radialavg` is true, respectively.

`nfft` specifies the number of points to use for the Fourier
transform. If `size(s)` < `nfft`, then the input is padded
with zeros. By default, `nfft` is the closest size for which the
Fourier transform can be computed with maximal efficiency. `fs`
is the sample rate of the original signal in both directions.

For `radialsum=true` the value of `power[k]` is proportional to
``\\frac{1}{N}\\sum_{k\\leq |k'|<k+1} |X[k']|^2``.
For `radialavg=true` it is proportional to
``\\frac{1}{N \\#\\{k\\leq |k'|<k+1\\}} \\sum_{k\\leq |k'|<k+1} |X[k']|^2``.
The computation of `|k'|` takes into account non-square signals
by scaling the coordinates of the wavevector accordingly.
"""
function periodogram(s::AbstractMatrix{T};
                     nfft::NTuple{2,Int}=nextfastfft(size(s)),
                     fs::Real=1,
                     radialsum::Bool=false, radialavg::Bool=false) where T<:Real
    size(s,1)<=nfft[1] && size(s,2)<=nfft[2] || throw(ArgumentError("nfft must be >= size(s)"))
    size(s,1)>1 && size(s,2)>1 || throw(ArgumentError("dimensions of s must be > 1"))
    if radialsum && radialavg
        throw(ArgumentError("radialsum and radialavg are mutually exclusive"))
    elseif !radialsum && !radialavg
        ptype = 0
    elseif radialsum
        ptype = 1
    elseif radialavg
        ptype = 2
    end
    norm2 = length(s)
    nmin = minimum(nfft)

    if prod(nfft) == length(s) && isa(s, StridedArray)
        input = s # no need to pad
    else
        input = zeros(fftintype(T), nfft)
        input[1:size(s,1), 1:size(s,2)] = s
    end

    if ptype == 0
        s_fft = fft(input)
        out = zeros(fftabs2type(T), nfft)
        fft2pow2!(out,s_fft,nfft...,fs*norm2)
        return Periodogram2(out, fftfreq(nfft[1],fs), fftfreq(nfft[2],fs))
    else
        s_fft = rfft(input)
        out = zeros(fftabs2type(T), nmin>>1 + 1)
        fft2pow2radial!(out,s_fft,nfft...,fs*norm2, ptype)
        return Periodogram(out, Frequencies(length(out), length(out), fs/nmin))
    end
end

forward_plan(X::AbstractArray{T}, ::AbstractArray{Complex{T}}) where {T<:Union{Float32, Float64}} =
    plan_rfft(X)
forward_plan(X::AbstractArray{T}, ::AbstractArray{T}) where {T<:Union{ComplexF32, ComplexF64}} =
    plan_fft(X)

# Compute an estimate of the power spectral density of a signal s via Welch's
# method.  The resulting periodogram has length N and is computed with an overlap
# region of length M.  The method is detailed in "The Use of Fast Fourier Transform
# for the Estimation of Power Spectra: A Method based on Time Averaging over Short,
# Modified Periodograms."  P. Welch, IEEE Transactions on Audio and Electroacoustics,
# vol AU-15, pp 70-73, 1967.
"""
    welch_pgram(s, n=div(length(s), 8), noverlap=div(n, 2); onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, window=nothing)

Computes the Welch periodogram of a signal `s` based on segments with `n` samples
with overlap of `noverlap` samples, and returns a Periodogram
object. For a Bartlett periodogram, set `noverlap=0`. See
[`periodogram`](@ref) for description of optional keyword arguments.
"""
function welch_pgram(s::AbstractVector{T}, n::Int=length(s)>>3, noverlap::Int=n>>1;
                     onesided::Bool=eltype(s)<:Real,
                     nfft::Int=nextfastfft(n), fs::Real=1,
                     window::Union{Function,AbstractVector,Nothing}=nothing) where T<:Number
    onesided && T <: Complex && throw(ArgumentError("cannot compute one-sided FFT of a complex signal"))
    nfft >= n || throw(DomainError((nfft=nfft, n=n), "nfft must be >= n"))

    win, norm2 = compute_window(window, n)
    sig_split = arraysplit(s, n, noverlap, nfft, win)
    out = zeros(fftabs2type(T), onesided ? (nfft >> 1)+1 : nfft)
    r = fs*norm2*length(sig_split)

    tmp = Vector{fftouttype(T)}(undef, T<:Real ? (nfft >> 1)+1 : nfft)
    plan = forward_plan(sig_split.buf, tmp)
    for sig in sig_split
        mul!(tmp, plan, sig)
        fft2pow!(out, tmp, nfft, r, onesided)
    end

    Periodogram(out, onesided ? rfftfreq(nfft, fs) : fftfreq(nfft, fs))
end

## SPECTROGRAM

const Float64Range = typeof(range(0.0, step=1.0, length=2))

struct Spectrogram{T,F<:Union{Frequencies,AbstractRange}, M<:AbstractMatrix{T}} <: TFR{T}
    power::M
    freq::F
    time::Float64Range
end
FFTW.fftshift(p::Spectrogram{T,F}) where {T,F<:Frequencies} =
    Spectrogram(p.freq.n_nonnegative == p.freq.n ? p.power : fftshift(p.power, 1), fftshift(p.freq), p.time)
FFTW.fftshift(p::Spectrogram{T,F}) where {T,F<:AbstractRange} = p

"""
    time(p)

Returns the time bin centers for a given Spectrogram object.
"""
Base.time(p::Spectrogram) = p.time

"""
    spectrogram(s, n=div(length(s), 8), noverlap=div(n, 2); onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, window=nothing)

Computes the spectrogram of a signal `s` based on segments with `n` samples
with overlap of `noverlap` samples, and returns a Spectrogram object. See
[`periodogram`](@ref) for description of optional keyword arguments.
"""
function spectrogram(s::AbstractVector{T}, n::Int=length(s)>>3, noverlap::Int=n>>1;
                     onesided::Bool=eltype(s)<:Real,
                     nfft::Int=nextfastfft(n), fs::Real=1,
                     window::Union{Function,AbstractVector,Nothing}=nothing) where T

    out = stft(s, n, noverlap, PSDOnly(); onesided=onesided, nfft=nfft, fs=fs, window=window)
    Spectrogram(out, onesided ? rfftfreq(nfft, fs) : fftfreq(nfft, fs),
                (n/2 : n-noverlap : (size(out,2)-1)*(n-noverlap)+n/2) / fs)

end

struct PSDOnly end
stfttype(T::Type, psdonly::PSDOnly) = fftabs2type(T)
stfttype(T::Type, psdonly::Nothing) = fftouttype(T)

"""
    stft(s, n=div(length(s), 8), noverlap=div(n, 2); onesided=eltype(s)<:Real, nfft=nextfastfft(n), fs=1, window=nothing)

Computes the STFT of a signal `s` based on segments with `n` samples
with overlap of `noverlap` samples, and returns a matrix containing the STFT
coefficients. See [`periodogram`](@ref) for description of optional
keyword arguments.
"""
function stft(s::AbstractVector{T}, n::Int=length(s)>>3, noverlap::Int=n>>1,
              psdonly::Union{Nothing,PSDOnly}=nothing;
              onesided::Bool=eltype(s)<:Real, nfft::Int=nextfastfft(n), fs::Real=1,
              window::Union{Function,AbstractVector,Nothing}=nothing) where T
    onesided && T <: Complex && throw(ArgumentError("cannot compute one-sided FFT of a complex signal"))

    win, norm2 = compute_window(window, n)
    sig_split = arraysplit(s, n, noverlap, nfft, win)
    nout = onesided ? (nfft >> 1)+1 : nfft
    out = zeros(stfttype(T, psdonly), nout, length(sig_split))
    tmp = Vector{fftouttype(T)}(undef, T<:Real ? (nfft >> 1)+1 : nfft)
    r = fs*norm2

    plan = forward_plan(sig_split.buf, tmp)
    offset = 0
    for sig in sig_split
        mul!(tmp, plan, sig)
        if isa(psdonly, PSDOnly)
            fft2pow!(out, tmp, nfft, r, onesided, offset)
        else
            fft2oneortwosided!(out, tmp, nfft, onesided, offset)
        end
        offset += nout
    end
    out
end

include("multitaper.jl")

end # end module definition
