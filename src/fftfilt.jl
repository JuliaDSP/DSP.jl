module FFTFilt
export fftfilt, firfilt

const FFT_LENGTHS = 2.^(1:28)
# FFT times computed on a Core i7-3930K @4.4GHz
# The real time doesn't matter, just the relative difference
const FFT_TIMES = [6.36383e-7, 6.3779e-7 , 6.52212e-7, 6.65282e-7, 7.12794e-7, 7.63172e-7,
                   7.91914e-7, 1.02289e-6, 1.37939e-6, 2.10868e-6, 4.04436e-6, 9.12889e-6,
                   2.32142e-5, 4.95576e-5, 0.000124927, 0.000247771, 0.000608867, 0.00153119,
                   0.00359037, 0.0110568, 0.0310893, 0.065813, 0.143516, 0.465745, 0.978072,
                   2.04371, 4.06017, 8.77769]
const FAST_FFT_FACTORS = [2, 3, 5, 7]

# Determine optimal length of the FFT for fftfilt
function optimalfftfiltlength(nb, nx)
    nfft = 0
    if nb > FFT_LENGTHS[end] || nb >= nx
        nfft = nextprod(FAST_FFT_FACTORS, nx+nb-1)
    else
        fastestestimate = Inf
        firsti = max(1, searchsortedfirst(FFT_LENGTHS, nb))
        lasti = max(1, searchsortedfirst(FFT_LENGTHS, nx+nb-1))
        L = 0
        for i = firsti:lasti
            curL = FFT_LENGTHS[i] - (nb - 1)
            estimate = iceil(nx/curL)*FFT_TIMES[i]
            if estimate < fastestestimate
                nfft = FFT_LENGTHS[i]
                fastestestimate = estimate
                L = curL
            end
        end

        if L > nx
            # If L > nx, better to find next fast power
            nfft = nextprod(FAST_FFT_FACTORS, nx+nb-1)
        end
    end
    nfft
end

# Filter x using FIR filter b by overlap-save method
function fftfilt{T<:Real}(b::Vector{T}, x::Vector{T},
                          nfft=optimalfftfiltlength(length(b), length(x)))
    nb = length(b)
    nx = length(x)

    L = min(nx, nfft - (nb - 1))
    tmp1 = Array(T, nfft)
    tmp2 = Array(Complex{T}, nfft << 1 + 1)
    out = zeros(T, nx)

    p1 = FFTW.Plan(tmp1, tmp2, 1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT)
    p2 = FFTW.Plan(tmp2, tmp1, 1, FFTW.ESTIMATE, FFTW.NO_TIMELIMIT)

    # FFT of filter
    filterft = similar(tmp2)
    copy!(tmp1, b)
    tmp1[nb+1:end] = zero(T)
    FFTW.execute(p1.plan, tmp1, filterft)

    # FFT of chunks
    off = 1
    while off <= nx
        npadbefore = max(0, nb - off)
        xstart = off - nb + npadbefore + 1
        n = min(nfft - npadbefore, nx - xstart + 1)

        tmp1[1:npadbefore] = zero(T)
        tmp1[npadbefore+n+1:end] = zero(T)

        copy!(tmp1, npadbefore+1, x, xstart, n)
        FFTW.execute(T, p1.plan)
        broadcast!(*, tmp2, tmp2, filterft)
        FFTW.execute(T, p2.plan)
        copy!(out, off, tmp1, nb, min(L, nx - off + 1))

        off += L
    end

    # Normalize
    scale!(out, 1/nfft)
end

# Filter x using FIR filter b, heuristically choosing to perform
# convolution in the time domain using filt or in the frequency domain
# using fftfilt
function firfilt{T<:Number}(b::AbstractVector{T}, x::AbstractVector{T})
    nb = length(b)
    nx = length(x)

    filtops = nx * min(nx, nb)
    if filtops <= 100000
        # 65536 is apprximate cutoff where FFT-based algorithm may be
        # more effective (due to overhead for allocation, plan
        # creation, etc.)
        filt(b, [one(T)], x)
    else
        # Estimate number of multiplication operations for fftfilt()
        # and filt()
        nfft = optimalfftfiltlength(nb, nx)
        L = min(nx, nfft - (nb - 1))
        nchunk = iceil(nx/L)
        fftops = (2*nchunk + 1) * nfft * log2(nfft)/2 + nchunk * nfft + 100000

        filtops > fftops ? fftfilt(b, x, nfft) : filt(b, [one(T)], x)
    end
end
end