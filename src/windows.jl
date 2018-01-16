module Windows
using ..DSP: @importffts
using ..Util
import SpecialFunctions: besseli
using Compat: copyto!, uninitialized
@importffts

export  rect,
        hanning,
        hamming,
        tukey,
        cosine,
        lanczos,
        triang,
        bartlett,
        gaussian,
        bartlett_hann,
        blackman,
        kaiser,
        dpss,
        dpsseig
#
# Window functions
#

"""
    rect(n)

Rectangular window function of length `n`.
"""
function rect(n::Integer)
    ones(n)
end

"""
    hanning(n)

Hanning window of length `n`.
"""
function hanning(n::Integer)
    [0.5*(1 - cos(2*pi*k/(n-1))) for k=0:(n-1)]
end

"""
    hamming(n)

Hamming window of length `n`.
"""
function hamming(n::Integer)
    [0.54 - 0.46*cos(2*pi*k/(n-1)) for k=0:(n-1)]
end

"""
    tukey(n, alpha)

Tukey window of length `n`, parameterized by `alpha`. For
`alpha` = 0, the window is equivalent to a rectangular window.
For `alpha` = 1, the window is a Hann window.
"""
function tukey(n::Integer, alpha::Real)
    # check that alpha is reasonable
    !(0 <= alpha <= 1) && error("tukey window alpha parameter must be 0 <= alpha <= 1.")

    # if alpha is less than machine precision, call it zero and return the
    # rectangular window for this length.  if we don't short circuit this
    # here, it will blow up below.
    t = zeros(n)
    if abs(alpha) <= eps()
        t = rect(n)
    else
        m = alpha*(n-1)/2
        for k=0:(n-1)
            if k <= m
                t[k+1] = 0.5*(1 + cos(pi*(k/m - 1)))
            elseif k <= n-1-m
                t[k+1] = 1
            else
                t[k+1] = 0.5*(1 + cos(pi*(k/m - 2/alpha + 1)))
            end
        end
    end

    return t
end

"""
    cosine(n)

Cosine window of length `n`. Also called the sine window for
obvious reasons.
"""
function cosine(n::Integer)
    [sin(pi*k/(n-1)) for k=0:(n-1)]
end

"""
    lanczos(n)

Lanczos window of length `n`.
"""
function lanczos(n::Integer)
    [sinc(2*k/(n-1) - 1) for k=0:(n-1)]
end

"""
    triang(n)

Triangular window of length `n`.
"""
function triang(n::Integer)
    [1 - abs((k - (n-1)/2))/(n/2) for k=0:(n-1)]
end

"""
    bartlett(n)

Bartlett window of length `n`.
"""
function bartlett(n::Integer)
    [2/(n-1)*((n-1)/2 - abs(k - (n-1)/2)) for k=0:(n-1)]
end

"""
    gaussian(n, sigma)

Gaussian window of length `n` parameterized by the standard
deviation `sigma`.
"""
function gaussian(n::Integer, sigma::Real)
    if !(0 < sigma <= 0.5)
        error("sigma must be greater than 0 and less than or equal to 0.5.")
    end
    [exp(-0.5*((k-(n-1)/2)/(sigma*(n-1/2)))^2) for k=0:(n-1)]
end

"""
    bartlett_hann(n)

Bartlett-Hann window of length `n`.
"""
function bartlett_hann(n::Integer)
    a0, a1, a2 = 0.62, 0.48, 0.38
    t = 2*pi/(n-1)
    [a0 - a1*abs(k/(n-1) - 0.5) - a2*cos(t*k) for k=0:(n-1)]
end

"""
    blackman(n)

"Exact" Blackman window, alpha = 0.16.
"""
function blackman(n::Integer)
    a0, a1, a2 = 0.42, 0.5, 0.08
    t = 2*pi/(n-1)
    [a0 - a1*cos(t*k) + a2*cos(t*k*2) for k=0:(n-1)]
end

"""
    kaiser(n, alpha)

Kaiser window of length `n` parameterized by `alpha`.
"""
function kaiser(n::Integer, alpha::Real)
    pf = 1.0/besseli(0,pi*alpha)
    [pf*besseli(0, pi*alpha*(sqrt(1 - (2*k/(n-1) - 1)^2))) for k=0:(n-1)]
end

# Discrete prolate spheroid sequences (Slepian tapers)
#
# See Gruenbacher, D. M., & Hummels, D. R. (1994). A simple algorithm
# for generating discrete prolate spheroidal sequences. IEEE
# Transactions on Signal Processing, 42(11), 3276-3278.
"""
    dpss(n, nw, ntapers=iceil(2*nw)-1)

The first `ntapers` discrete prolate spheroid sequences (Slepian
tapers) as an `n` × `ntapers` matrix. The signs of the tapers
follow the convention that the first element of the skew-symmetric
(odd) tapers is positive. The time-bandwidth product is given by
`nw`.
"""
function dpss(n::Int, nw::Real, ntapers::Int=ceil(Int, 2*nw)-1)
    0 < ntapers <= n || error("ntapers must be in interval (0, n]")
    0 <= nw < n/2 || error("nw must be in interval [0, n/2)")

    # Construct symmetric tridiagonal matrix
    v = cospi(2*nw/n)
    mat = SymTridiagonal([v*abs2((n - 1)/2 - i) for i=0:(n-1)],
                         [0.5.*(i*n - abs2(i)) for i=1:(n-1)])

    # Get tapers
    v = flipdim(eigfact!(mat, n-ntapers+1:n)[:vectors]::Matrix{Float64}, 2)

    # Slepian's convention; taper starts with a positive element
    sgn = ones(size(v, 2))
    for i = 2:2:size(v, 2)
        s = 0
        for j = 1:n
            s = sign(v[j, i])
            s != 0 && break
        end
        @assert s != 0
        sgn[i] = s
    end
    scale!(v, sgn)
end

# Eigenvalues of DPSS, following Percival & Walden p. 390, exercise 8.1
# See also implementation in MNE:
# https://github.com/mne-tools/mne-python/blob/d7082cf909ccab581667bc1f1ed3c23e6a24b567/mne/time_frequency/multitaper.py#L226
"""
    dpsseig(A, nw)

Eigenvalues of the DPSS matrix, representing the ratios of the
power within the main lobe to the total power (main and sidelobes).
`A` is the output of [`dpss`](@ref), and `nw` is the
time-bandwidth product provided to [`dpss`](@ref) as input.
"""
function dpsseig(A::Matrix{Float64}, nw::Real)
    0 <= nw < size(A, 1)/2 || error("nw must be in interval [0, n/2)")

    w = nw/size(A, 1)

    # Compute coefficients
    seq = Vector{Float64}(uninitialized, size(A, 1))
    seq[1] = 1.0
    for i = 1:size(A, 1)-1
        seq[i+1] = 2 * sinc(2w*i)
    end

    q = Vector{Float64}(uninitialized, size(A, 2))
    nfft = nextfastfft(2*size(A, 1)-1)

    tmp1 = Vector{Float64}(uninitialized, nfft)
    tmp2 = Vector{Complex{Float64}}(uninitialized, nfft >> 1 + 1)
    p1 = plan_rfft(tmp1)
    p2 = plan_brfft(tmp2, nfft)

    for i = 1:size(A, 2)
        fill!(tmp1, 0)
        copyto!(tmp1, 1, A, (i-1)*size(A, 1)+1, size(A, 1))
        A_mul_B!(tmp2, p1, tmp1)
        for j = 1:length(tmp2)
            @inbounds tmp2[j] = abs2(tmp2[j])
        end
        A_mul_B!(tmp1, p2, tmp2)

        eig = 0.0
        for j = 1:size(A, 1)
            eig += seq[j]*tmp1[j]
        end
        q[i] = 2w * eig / nfft
    end

    q
end

# tensor product window functions in 2-d, defined as w(x,y) = w(x) \otimes w(y)
for func in (:tukey, :gaussian, :kaiser)
    @eval begin
        $func(n::NTuple{2,Integer}, a::Real) = $func(n[1], a) * $func(n[2], a)'
    end
end
for func in (:rect, :hanning, :hamming, :cosine, :lanczos,
             :triang, :bartlett, :bartlett_hann, :blackman)
    @eval begin
        $func(n::NTuple{2,Integer}) = $func(n[1]) * $func(n[2])'
    end
end


end # end module definition
