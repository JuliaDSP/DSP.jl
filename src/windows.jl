module Windows
using ..DSP: @importffts, mul!, rmul!
using ..Util
import SpecialFunctions: besseli
import Compat
using Compat: copyto!, undef
if VERSION < v"0.7.0-DEV.5211"
    using Compat.LinearAlgebra: Diagonal, SymTridiagonal, eigfact!
else
    using LinearAlgebra: Diagonal, SymTridiagonal, eigen!
end
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

"""
    makewindow(winfunc::Function, len, padding, zerophase)

Generate a discrete window vector of the given length with `padding` zeros.
`winfunc` should be a function giving the window value in the range of
[-0.5, 0.5]. The function is assumed to be 0 outside of this range.

If `zerophase` is `true` the window is centered around index 1 (with the
negative half wrapped to the end of the vector), and is suitable for FFT
processing. These are usually even-length.
If `zerophase` is `false` the window is centered around `(n+1)/2`, which is
commonly used for FIR filter design. These are usually odd-length.
"""
function makewindow(winfunc::Function, len, padding, zerophase)
    win = zeros(len+padding)
    # TODO - handle odd-length zerophase case (e.g. check triang)
    if zerophase
        win[end-len÷2+1:end] .= winfunc.(linspace(-0.5, -1/len, len÷2))
        win[1:len÷2] .= winfunc.(linspace(0.0, 0.5-1/len, len÷2))
    else
        win[1:len] .= winfunc.(linspace(-0.5, 0.5, len))
    end

    win
end

#
# Window functions
#

"""
    rect(n; padding=0, zerophase=false)

Rectangular window function of length `n`, padded with `padding` zeros. This
window is 1 within the window, and 0 outside of it.

If `zerophase` is `true` the window is centered around index 1 (with the
negative half wrapped to the end of the vector), and is suitable for FFT
processing. These are usually even-length.
If `zerophase` is `false` the window is centered around `(n+1)/2`, which is
commonly used for FIR filter design. These are usually odd-length.
"""
function rect(n::Integer; padding=0, zerophase=false)
    makewindow(n, padding, zerophase) do x
        1.0
    end
end


"""
    hanning(n)

Hanning window of length `n`.
"""
function hanning(n::Integer; padding=0, zerophase=false)
    makewindow(n, padding, zerophase) do x
        0.5*(1+cos(2pi*x))
    end
end

"""
    hamming(n)

Hamming window of length `n`.
"""
function hamming(n::Integer; padding=0, zerophase=false)
    makewindow(n, padding, zerophase) do x
        0.54 + 0.46*cos(2*pi*x)
    end
end

"""
    tukey(n, alpha)

Tukey window of length `n`, parameterized by `alpha`. For
`alpha` = 0, the window is equivalent to a rectangular window.
For `alpha` = 1, the window is a Hann window.
"""
function tukey(n::Integer, alpha::Real; padding=0, zerophase=false)
    # check that alpha is reasonable
    !(0 <= alpha <= 1) && error("tukey window alpha parameter must be 0 <= alpha <= 1.")

    # if alpha is less than machine precision, call it zero and return the
    # rectangular window for this length.  if we don't short circuit this
    # here, it will blow up below.
    if abs(alpha) <= eps()
        rect(n, padding, zerophase)
    else
        m = alpha/2
        makewindow(n, padding, zerophase) do x
            # shift x so we define in terms of the range [0,1]
            x += 0.5
            if x <= m
                0.5*(1 + cos(pi*(x/m - 1)))
            elseif x <= 1-m
                1.0
            else
                0.5*(1 + cos(pi*(x/m - 2/alpha + 1)))
            end
        end
    end
end

"""
    cosine(n)

Cosine window of length `n`. Also called the sine window for
obvious reasons.
"""
function cosine(n::Integer; padding=0, zerophase=false)
    makewindow(n, padding, zerophase) do x
        cos(pi*x)
    end
end

"""
    lanczos(n)

Lanczos window of length `n`.
"""
function lanczos(n::Integer; padding=0, zerophase=false)
    makewindow(n, padding, zerophase) do x
        sinc(2x)
    end
end

"""
    triang(n)

Triangular window of length `n`.
"""
function triang(n::Integer; padding=0, zerophase=false)
    makewindow(n, padding, zerophase) do x
        1 - (n-1)/n*abs(2x)
    end
end

"""
    bartlett(n)

Bartlett window of length `n`. The Bartlett window is a triangular window that
reaches 0 at the edges.
"""
function bartlett(n::Integer; padding=0, zerophase=false)
    makewindow(n, padding, zerophase) do x
        1 - abs(2x)
    end
end

"""
    gaussian(n, σ)

Gives an n-sample gaussian window defined by sampling the function
\$w(x) = e^{-\\frac 1 2 \\left(\\frac x σ \\right)^2}\$ in the range
\$[-0.5,0.5]\$. This means that for \$σ=0.5\$ the endpoints of the window will
correspond to 1 standard deviation away from the center.
"""
function gaussian(n::Integer, σ::Real; padding=0, zerophase=false)
    σ > 0.0 || error("σ must be positive")
    makewindow(n, padding, zerophase) do x
        exp(-0.5*(x/σ)^2)
    end
end

"""
    bartlett_hann(n)

Bartlett-Hann window of length `n`.
"""
function bartlett_hann(n::Integer; padding=0, zerophase=false)
    a0, a1, a2 = 0.62, 0.48, 0.38
    makewindow(n, padding, zerophase) do x
        a0 - a1*abs(x) + a2*cos(2pi*x)
    end
end

"""
    blackman(n)

"Exact" Blackman window, alpha = 0.16.
"""
function blackman(n::Integer; padding=0, zerophase=false)
    a0, a1, a2 = 0.42, 0.5, 0.08
    makewindow(n, padding, zerophase) do x
        a0 + a1*cos(2pi*x) + a2*cos(4pi*x)
    end
end

"""
    kaiser(n, alpha)

Kaiser window of length `n` parameterized by `alpha`.
"""
function kaiser(n::Integer, alpha::Real; padding=0, zerophase=false)
    pf = 1.0/besseli(0,pi*alpha)
    makewindow(n, padding, zerophase) do x
        pf*besseli(0, pi*alpha*(sqrt(1 - (2x)^2)))
    end
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
function dpss(n::Integer, nw::Real, ntapers::Integer=ceil(Int, 2*nw)-1;
              padding=0, zerophase=false)
    # TODO there's probably a more efficient way to compute the zero-phase
    # version of this
    if zerophase
        n += 1
    end
    0 < ntapers <= n || error("ntapers must be in interval (0, n]")
    0 <= nw < n/2 || error("nw must be in interval [0, n/2)")

    # Construct symmetric tridiagonal matrix
    v = cospi(2*nw/n)
    dv = Vector{Float64}(undef, n)
    ev = Vector{Float64}(undef, n - 1)
    @inbounds dv[1] = v * abs2((n - 1) / 2)
    @inbounds @simd for i = 1:(n-1)
        dv[i + 1] = v * abs2((n - 1) / 2 - i)
        ev[i] = 0.5 * (i * n - i^2)
    end
    mat = SymTridiagonal(dv, ev)

    # Get tapers
    @static if VERSION < v"0.7.0-DEV.3159"
        eigvec = eigfact!(mat, n-ntapers+1:n)[:vectors]::Matrix{Float64}
    elseif VERSION < v"0.7.0-DEV.5211"
        eigvec = eigfact!(mat, n-ntapers+1:n).vectors
    else
        eigvec = eigen!(mat, n-ntapers+1:n).vectors
    end
    rv = Compat.reverse(eigvec, dims=2)::Matrix{Float64}

    # Slepian's convention; taper starts with a positive element
    sgn = ones(size(rv, 2))
    for i = 2:2:size(rv, 2)
        s = zero(eltype(rv))
        for j = 1:n
            s = sign(rv[j, i])
            s != 0 && break
        end
        @assert s != 0
        sgn[i] = s
    end
    rmul!(rv, Diagonal(sgn))

    if zerophase
        # remove the last element (should be equal to the first element)
        rv = rv[1:end-1, :]
    end

    if padding > 0
        rv = [rv; zeros(padding, ntapers)]
    end

    if zerophase
        # TODO make sure this is doing the right thing for odd-length windows
        rv = circshift(rv, (-n÷2, 0))
    end

    rv
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
    seq = Vector{Float64}(undef, size(A, 1))
    seq[1] = 1.0
    for i = 1:size(A, 1)-1
        seq[i+1] = 2 * sinc(2w*i)
    end

    q = Vector{Float64}(undef, size(A, 2))
    nfft = nextfastfft(2*size(A, 1)-1)

    tmp1 = Vector{Float64}(undef, nfft)
    tmp2 = Vector{Complex{Float64}}(undef, nfft >> 1 + 1)
    p1 = plan_rfft(tmp1)
    p2 = plan_brfft(tmp2, nfft)

    for i = 1:size(A, 2)
        fill!(tmp1, 0)
        copyto!(tmp1, 1, A, (i-1)*size(A, 1)+1, size(A, 1))
        mul!(tmp2, p1, tmp1)
        for j = 1:length(tmp2)
            @inbounds tmp2[j] = abs2(tmp2[j])
        end
        mul!(tmp1, p2, tmp2)

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
        $func(n::NTuple{2,Integer}, a::Real; padding=(0,0), kwargs...) =
            $func(n[1], a; padding=padding[1], kwargs...) * $func(n[2], a; padding=padding[2], kwargs...)'
    end
end
for func in (:rect, :hanning, :hamming, :cosine, :lanczos,
             :triang, :bartlett, :bartlett_hann, :blackman)
    @eval begin
        $func(n::NTuple{2,Integer}; padding=(0,0), kwargs...) =
            $func(n[1]; padding=padding[1], kwargs...) * $func(n[2]; padding=padding[2], kwargs...)'
    end
end


end # end module definition
