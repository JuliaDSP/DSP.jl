# Implementations of filt, filtfilt, and fftfilt

#
# filt and filt!
#

## PolynomialRatio
_zerosi{T,S}(f::PolynomialRatio{T}, x::AbstractArray{S}) =
    zeros(promote_type(T, S), max(length(f.a), length(f.b))-1)

Base.filt!{T,S}(out, f::PolynomialRatio{T}, x::AbstractArray{S}, si=_zerosi(f, x)) =
    filt!(out, coefb(f), coefa(f), x, si)
Base.filt(f::PolynomialRatio, x, si=_zerosi(f, x)) = filt(coefb(f), coefa(f), x, si)

## SecondOrderSections
_zerosi{T,G,S}(f::SecondOrderSections{T,G}, x::AbstractArray{S}) =
    zeros(promote_type(T, G, S), 2, length(f.biquads))

# filt! algorithm (no checking, returns si)
function _filt!{S,N}(out::AbstractArray, si::AbstractArray{S,N}, f::SecondOrderSections,
                     x::AbstractArray, col::Int)
    g = f.g
    biquads = f.biquads
    n = length(biquads)
    @inbounds for i = 1:size(x, 1)
        yi = x[i, col]
        for fi = 1:n
            biquad = biquads[fi]
            xi = yi
            yi = si[1, fi] + biquad.b0*xi
            si[1, fi] = si[2, fi] + biquad.b1*xi - biquad.a1*yi
            si[2, fi] = biquad.b2*xi - biquad.a2*yi
        end
        out[i, col] = yi*g
    end
    si
end

function Base.filt!{S,N}(out::AbstractArray, f::SecondOrderSections, x::AbstractArray,
                         si::AbstractArray{S,N}=_zerosi(f, x))
    biquads = f.biquads
    ncols = Base.trailingsize(x, 2)

    size(x) != size(out) && error("out size must match x")
    (size(si, 1) != 2 || size(si, 2) != length(biquads) || (N > 2 && Base.trailingsize(si, 3) != ncols)) &&
        error("si must be 2 x nbiquads or 2 x nbiquads x nsignals")

    initial_si = si
    g = f.g
    n = length(biquads)
    for col = 1:ncols
        _filt!(out, initial_si[:, :, N > 2 ? col : 1], f, x, col)
    end
    out
end

Base.filt{T,G,S<:Number}(f::SecondOrderSections{T,G}, x::AbstractArray{S}, si=_zerosi(f, x)) =
    filt!(Array{promote_type(T, G, S)}(size(x)), f, x, si)

## Biquad
_zerosi{T,S}(f::Biquad{T}, x::AbstractArray{S}) =
    zeros(promote_type(T, S), 2)

# filt! algorithm (no checking, returns si)
function _filt!(out::AbstractArray, si1::Number, si2::Number, f::Biquad,
                x::AbstractArray, col::Int)
    @inbounds for i = 1:size(x, 1)
        xi = x[i, col]
        yi = si1 + f.b0*xi
        si1 = si2 + f.b1*xi - f.a1*yi
        si2 = f.b2*xi - f.a2*yi
        out[i, col] = yi
    end
    (si1, si2)
end

# filt! variant that preserves si
function Base.filt!{S,N}(out::AbstractArray, f::Biquad, x::AbstractArray,
                         si::AbstractArray{S,N}=_zerosi(f, x))
    ncols = Base.trailingsize(x, 2)

    size(x) != size(out) && error("out size must match x")
    (size(si, 1) != 2 || (N > 1 && Base.trailingsize(si, 2) != ncols)) &&
        error("si must have two rows and 1 or nsignals columns")

    initial_si = si
    for col = 1:ncols
        _filt!(out, initial_si[1, N > 1 ? col : 1], initial_si[2, N > 1 ? col : 1], f, x, col)
    end
    out
end

Base.filt{T,S<:Number}(f::Biquad{T}, x::AbstractArray{S}, si=_zerosi(f, x)) =
    filt!(Array{promote_type(T, S)}(size(x)), f, x, si)

## For arbitrary filters, convert to SecondOrderSections
Base.filt(f::FilterCoefficients, x) = filt(convert(SecondOrderSections, f), x)
Base.filt!(out, f::FilterCoefficients, x) = filt!(out, convert(SecondOrderSections, f), x)

#
# Direct form II transposed filter with state
#

immutable DF2TFilter{T<:FilterCoefficients,S<:Array}
    coef::T
    state::S

    function DF2TFilter(coef::PolynomialRatio, state::Vector)
        length(state) == length(coef.a)-1 == length(coef.b)-1 ||
            throw(ArgumentError("length of state vector must match filter order"))
        new(coef, state)
    end
    function DF2TFilter(coef::SecondOrderSections, state::Matrix)
        (size(state, 1) == 2 && size(state, 2) == length(coef.biquads)) ||
            throw(ArgumentError("state must be 2 x nbiquads"))
        new(coef, state)
    end
    function DF2TFilter(coef::Biquad, state::Vector)
        length(state) == 2 || throw(ArgumentError("length of state must be 2"))
        new(coef, state)
    end
end

## PolynomialRatio
DF2TFilter{T,S}(coef::PolynomialRatio{T},
                state::Vector{S}=zeros(T, max(length(coef.a), length(coef.b))-1)) =
    DF2TFilter{PolynomialRatio{T}, Vector{S}}(coef, state)

function Base.filt!{T,S}(out::AbstractVector, f::DF2TFilter{PolynomialRatio{T},Vector{S}}, x::AbstractVector)
    length(x) != length(out) && throw(ArgumentError("out size must match x"))

    si = f.state
    # Note: these are in the Polynomials.jl convention, which is
    # reversed WRT the usual DSP convention
    b = f.coef.b.a
    a = f.coef.a.a
    n = length(b)
    if n == 1
        scale!(out, x, b[1])
    else
        @inbounds for i=1:length(x)
            xi = x[i]
            val = si[1] + b[n]*xi
            for j=1:n-2
                si[j] = si[j+1] + b[n-j]*xi - a[n-j]*val
            end
            si[n-1] = b[1]*xi - a[1]*val
            out[i] = val
        end
    end
    out
end

## SecondOrderSections
DF2TFilter{T,G,S}(coef::SecondOrderSections{T,G},
                  state::Matrix{S}=zeros(promote_type(T, G), 2, length(coef.biquads))) =
    DF2TFilter{SecondOrderSections{T,G}, Matrix{S}}(coef, state)

function Base.filt!{T,G,S}(out::AbstractVector, f::DF2TFilter{SecondOrderSections{T,G},Matrix{S}}, x::AbstractVector)
    length(x) != length(out) && throw(ArgumentError("out size must match x"))
    _filt!(out, f.state, f.coef, x, 1)
    out
end

## Biquad
DF2TFilter{T,S}(coef::Biquad{T}, state::Vector{S}=zeros(T, 2)) =
    DF2TFilter{Biquad{T}, Vector{S}}(coef, state)
function Base.filt!{T,S}(out::AbstractVector, f::DF2TFilter{Biquad{T},Vector{S}}, x::AbstractVector)
    length(x) != length(out) && throw(ArgumentError("out size must match x"))
    si = f.state
    (si[1], si[2]) =_filt!(out, si[1], si[2], f.coef, x, 1)
    out
end

# Variant that allocates the output
Base.filt{T,S<:Array}(f::DF2TFilter{T,S}, x::AbstractVector) =
    filt!(Array{eltype(S)}(length(x)), f, x)

# Fall back to SecondOrderSections
DF2TFilter(coef::FilterCoefficients) = DF2TFilter(convert(SecondOrderSections, coef))

#
# filtfilt
#

# Extrapolate the beginning of a signal for use by filtfilt. This
# computes:
#
# [(2 * x[1]) .- x[pad_length+1:-1:2],
#  x,
#  (2 * x[end]) .- x[end-1:-1:end-pad_length]]
#
# in place in output. The istart and n parameters determine the portion
# of the input signal x to extrapolate.
function extrapolate_signal!(out, ostart, sig, istart, n, pad_length)
    length(out) >= n+2*pad_length || error("output is incorrectly sized")
    x = 2*sig[istart]
    for i = 1:pad_length
        out[ostart+i-1] = x - sig[istart+pad_length+1-i]
    end
    copy!(out, ostart+pad_length, sig, istart, n)
    x = 2*sig[istart+n-1]
    for i = 1:pad_length
        out[ostart+n+pad_length+i-1] = x - sig[istart+n-1-i]
    end
    out
end

# Zero phase digital filtering by processing data in forward and reverse direction
function iir_filtfilt(b::AbstractVector, a::AbstractVector, x::AbstractArray)
    zi = filt_stepstate(b, a)
    zitmp = copy(zi)
    pad_length = 3 * (max(length(a), length(b)) - 1)
    t = Base.promote_eltype(b, a, x)
    extrapolated = Array{t}(size(x, 1)+pad_length*2)
    out = similar(x, t)

    istart = 1
    for i = 1:Base.trailingsize(x, 2)
        extrapolate_signal!(extrapolated, 1, x, istart, size(x, 1), pad_length)
        reverse!(filt!(extrapolated, b, a, extrapolated, scale!(zitmp, zi, extrapolated[1])))
        filt!(extrapolated, b, a, extrapolated, scale!(zitmp, zi, extrapolated[1]))
        for j = 1:size(x, 1)
            @inbounds out[j, i] = extrapolated[end-pad_length+1-j]
        end
        istart += size(x, 1)
    end

    out
end

# Zero phase digital filtering with an FIR filter in a single pass
function filtfilt(b::AbstractVector, x::AbstractArray)
    nb = length(b)
    # Only need as much padding as the order of the filter
    t = Base.promote_eltype(b, x)
    extrapolated = similar(x, t, size(x, 1)+2nb-2, Base.trailingsize(x, 2))

    # Convolve b with its reverse
    newb = filt(b, reverse(b))
    resize!(newb, 2nb-1)
    for i = 1:nb-1
        newb[nb+i] = newb[nb-i]
    end

    # Extrapolate each column
    for i = 1:size(extrapolated, 2)
        extrapolate_signal!(extrapolated, (i-1)*size(extrapolated, 1)+1, x, (i-1)*size(x, 1)+1, size(x, 1), nb-1)
    end

    # Filter
    out = filt(newb, extrapolated)

    # Drop garbage at start
    reshape(out[2nb-1:end, :], size(x))
end

# Choose whether to use FIR or iir_filtfilt depending on
# length of a
function filtfilt(b::AbstractVector, a::AbstractVector, x::AbstractArray)
    if length(a) == 1
        if a[1] != 1
            b /= a[1]
        end
        filtfilt(b, x)
    else
        iir_filtfilt(b, a, x)
    end
end

# Zero phase digital filtering for second order sections
function filtfilt{T,G,S}(f::SecondOrderSections{T,G}, x::AbstractArray{S})
    zi = filt_stepstate(f)
    zitmp = similar(zi)
    pad_length = 6 * length(f.biquads)
    t = Base.promote_type(T, G, S)
    extrapolated = Array{t}(size(x, 1)+pad_length*2)
    out = similar(x, t)

    istart = 1
    for i = 1:Base.trailingsize(x, 2)
        extrapolate_signal!(extrapolated, 1, x, istart, size(x, 1), pad_length)
        reverse!(filt!(extrapolated, f, extrapolated, scale!(zitmp, zi, extrapolated[1])))
        filt!(extrapolated, f, extrapolated, scale!(zitmp, zi, extrapolated[1]))
        for j = 1:size(x, 1)
            @inbounds out[j, i] = extrapolated[end-pad_length+1-j]
        end
        istart += size(x, 1)
    end

    out
end

# Support for other filter types
filtfilt(f::FilterCoefficients, x) = filtfilt(convert(SecondOrderSections, f), x)
filtfilt(f::PolynomialRatio, x) = filtfilt(coefb(f), coefa(f), x)

## Initial filter state

# Compute an initial state for filt with coefficients (b,a) such that its
# response to a step function is steady state.
function filt_stepstate{T<:Number}(b::@compat(Union{AbstractVector{T}, T}), a::@compat(Union{AbstractVector{T}, T}))
    scale_factor = a[1]
    if scale_factor != 1.0
        a = a ./ scale_factor
        b = b ./ scale_factor
    end

    bs = length(b)
    as = length(a)
    sz = max(bs, as)
    sz > 0 || error("a and b must have at least one element each")
    sz == 1 && return T[]

    # Pad the coefficients with zeros if needed
    bs<sz && (b = copy!(zeros(eltype(b), sz), b))
    as<sz && (a = copy!(zeros(eltype(a), sz), a))

    # construct the companion matrix A and vector B:
    A = [-a[2:end] [eye(T, sz-2); zeros(T, 1, sz-2)]]
    B = b[2:end] - a[2:end] * b[1]
    # Solve si = A*si + B
    # (I - A)*si = B
    scale_factor \ (I - A) \ B
 end

function filt_stepstate{T}(f::SecondOrderSections{T})
    biquads = f.biquads
    si = Array{T}(2, length(biquads))
    y = one(T)
    for i = 1:length(biquads)
        biquad = biquads[i]

        # At steady state, we have:
        #  y = s1 + b0*x
        # s1 = s2 + b1*x - a1*y
        # s2 = b2*x - a2*y
        # where x is the input and y is the output. Solving these
        # equations yields the following.
        si[1, i] = (-(biquad.a1 + biquad.a2)*biquad.b0 + biquad.b1 + biquad.b2)/
                   (1 + biquad.a1 + biquad.a2)*y
        si[2, i] = (biquad.a1*biquad.b2 - biquad.a2*(biquad.b0 + biquad.b1) + biquad.b2)/
                   (1 + biquad.a1 + biquad.a2)*y
        y *= (biquad.b0 + biquad.b1 + biquad.b2)/(1 + biquad.a1 + biquad.a2)
    end
    si
end

#
# filt implementation for FIR filters (faster than Base)
#

for n = 2:15
    silen = n-1
    si = [@compat(Symbol("si$i")) for i = 1:silen]
    @eval function Base.filt!{T}(out, b::NTuple{$n,T}, x)
        size(x) != size(out) && error("out size must match x")
        ncols = Base.trailingsize(x, 2)
        for col = 0:ncols-1
            $(Expr(:block, [:($(si[i]) = zero(b[$i])*zero(x[1])) for i = 1:silen]...))
            offset = col*size(x, 1)
            @inbounds for i=1:size(x, 1)
                xi = x[i+offset]
                val = $(si[1]) + b[1]*xi
                $(Expr(:block, [:($(si[j]) = $(si[j+1]) + b[$(j+1)]*xi) for j = 1:(silen-1)]...))
                $(si[silen]) = b[$(silen+1)]*xi
                out[i+offset] = val
            end
        end
        out
    end
end

chain = :(throw(ArgumentError("invalid tuple size")))
for n = 15:-1:2
    chain = quote
        if length(h) == $n
            filt!(out, ($([:(h[$i]) for i = 1:n]...),), x)
        else
            $chain
        end
    end
end

@eval function small_filt!{T}(out::AbstractArray, h::AbstractVector{T}, x::AbstractArray)
    $chain
end

function Base.filt!(out::AbstractArray, h::AbstractVector, x::AbstractArray)
    if length(h) == 1
        return scale!(out, h[1], x)
    elseif length(h) <= 15
        return small_filt!(out, h, x)
    end

    h = reverse(h)
    hLen = length(h)
    xLen = size(x, 1)
    size(x) != size(out) && error("out size must match x")
    ncols = Base.trailingsize(x, 2)

    for col = 0:ncols-1
        offset = col*size(x, 1)
        for i = 1:min(hLen-1, xLen)
            dotprod = zero(eltype(h))*zero(eltype(x))
            @simd for j = 1:i
                @inbounds dotprod += h[hLen-i+j] * x[j+offset]
            end
            @inbounds out[i+offset] = dotprod
        end
        for i = hLen:xLen
            @inbounds out[i+offset] = unsafe_dot(h, x, i+offset)
        end
    end
    out
end

Base.filt(h::AbstractArray, x::AbstractArray) =
    Base.filt!(Array{eltype(x)}(size(x)), h, x)

#
# fftfilt and filt
#

const FFT_LENGTHS = 2.^(1:28)
# FFT times computed on a Core i7-3930K @4.4GHz
# The real time doesn't matter, just the relative difference
const FFT_TIMES = [6.36383e-7, 6.3779e-7 , 6.52212e-7, 6.65282e-7, 7.12794e-7, 7.63172e-7,
                   7.91914e-7, 1.02289e-6, 1.37939e-6, 2.10868e-6, 4.04436e-6, 9.12889e-6,
                   2.32142e-5, 4.95576e-5, 0.000124927, 0.000247771, 0.000608867, 0.00153119,
                   0.00359037, 0.0110568, 0.0310893, 0.065813, 0.143516, 0.465745, 0.978072,
                   2.04371, 4.06017, 8.77769]

# Determine optimal length of the FFT for fftfilt
function optimalfftfiltlength(nb, nx)
    nfft = 0
    if nb > FFT_LENGTHS[end] || nb >= nx
        nfft = nextfastfft(nx+nb-1)
    else
        fastestestimate = Inf
        firsti = max(1, searchsortedfirst(FFT_LENGTHS, nb))
        lasti = max(1, searchsortedfirst(FFT_LENGTHS, nx+nb-1))
        L = 0
        for i = firsti:lasti
            curL = FFT_LENGTHS[i] - (nb - 1)
            estimate = ceil(Int, nx/curL)*FFT_TIMES[i]
            if estimate < fastestestimate
                nfft = FFT_LENGTHS[i]
                fastestestimate = estimate
                L = curL
            end
        end

        if L > nx
            # If L > nx, better to find next fast power
            nfft = nextfastfft(nx+nb-1)
        end
    end
    nfft
end

# Filter x using FIR filter b by overlap-save method
function fftfilt{T<:Real}(b::AbstractVector{T}, x::AbstractArray{T},
                          nfft=optimalfftfiltlength(length(b), length(x)))
    nb = length(b)
    nx = size(x, 1)
    normfactor = 1/nfft

    L = min(nx, nfft - (nb - 1))
    tmp1 = Array{T}(nfft)
    tmp2 = Array{Complex{T}}(nfft >> 1 + 1)
    out = Array{T}(size(x))

    @julia_newer_than v"0.4.0-dev+6068" begin
        p1 = plan_rfft(tmp1)
        p2 = plan_brfft(tmp2, nfft)
    end begin
        p1 = FFTW.Plan(tmp1, tmp2, 1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT)
        p2 = FFTW.Plan(tmp2, tmp1, 1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT)
    end

    # FFT of filter
    filterft = similar(tmp2)
    copy!(tmp1, b)
    tmp1[nb+1:end] = zero(T)
    @julia_newer_than(v"0.4.0-dev+6068",
                      A_mul_B!(filterft, p1, tmp1),
                      FFTW.execute(p1.plan, tmp1, filterft))

    # FFT of chunks
    for colstart = 0:nx:length(x)-1
        off = 1
        while off <= nx
            npadbefore = max(0, nb - off)
            xstart = off - nb + npadbefore + 1
            n = min(nfft - npadbefore, nx - xstart + 1)

            tmp1[1:npadbefore] = zero(T)
            tmp1[npadbefore+n+1:end] = zero(T)

            copy!(tmp1, npadbefore+1, x, colstart+xstart, n)
            @julia_newer_than(v"0.4.0-dev+6068",
                              A_mul_B!(tmp2, p1, tmp1),
                              FFTW.execute(T, p1.plan))
            broadcast!(*, tmp2, tmp2, filterft)
            @julia_newer_than(v"0.4.0-dev+6068",
                              A_mul_B!(tmp1, p2, tmp2),
                              FFTW.execute(T, p2.plan))

            # Copy to output
            for j = 0:min(L - 1, nx - off)
                @inbounds out[colstart+off+j] = tmp1[nb+j]*normfactor
            end

            off += L
        end
    end

    out
end

# Filter x using FIR filter b, heuristically choosing to perform
# convolution in the time domain using filt or in the frequency domain
# using fftfilt
function Base.filt{T<:Number}(b::AbstractVector{T}, x::AbstractArray{T})
    nb = length(b)
    nx = size(x, 1)

    filtops = length(x) * min(nx, nb)
    if filtops <= 500000
        # 65536 is apprximate cutoff where FFT-based algorithm may be
        # more effective (due to overhead for allocation, plan
        # creation, etc.)
        filt!(Array{eltype(x)}(size(x)), b, x)
    else
        # Estimate number of multiplication operations for fftfilt()
        # and filt()
        nfft = optimalfftfiltlength(nb, nx)
        L = min(nx, nfft - (nb - 1))
        nchunk = ceil(Int, nx/L)*div(length(x), nx)
        fftops = (2*nchunk + 1) * nfft * log2(nfft)/2 + nchunk * nfft + 500000

        filtops > fftops ? fftfilt(b, x, nfft) : filt!(Array{eltype(x)}(size(x)), b, x)
    end
end
