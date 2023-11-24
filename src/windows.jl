module Windows
using ..Util
import SpecialFunctions: besseli
using LinearAlgebra: Diagonal, SymTridiagonal, eigen!, mul!, rmul!
using FFTW

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

# this gets shared between the docstrings of all the window functions because
# they share these args.
const zerophase_docs = """
If `zerophase` is `false` (the default) the window is centered around index
`(n+1)/2`, which is commonly used for FIR filter design. These are often used
for FIR filter design, and usually odd-length. Note that for even-length windows
this will cause the window center to be located between samples.

If `zerophase` is `true` the window is centered around index 1 (with the
negative half wrapped to the end of the vector). Additionally this creates a
"periodic" window, which means that if there is no padding then the left and
right endpoints of the window wrap around to the same sample, so the window
length is the same as an `n+1`-length non-`zerophase` window. Alternatively you
can think of the continuous `zerophase` window being of width `n` and the
non-`zerophase` window being of length `n-1`. `zerophase` windows are often used
in FFT processing, and are usually even-length.
"""

function twoD_docs(arg=nothing)
    argstr = if arg !== nothing
        "`$arg`, "
    else
        ""
    end
    """
    Providing a `dims` `Tuple` rather than a single `n` constructs a 2D window.
    $argstr`padding` and `zerophase` can then be given either as a single value
    for both horizontal and vertical or a 2-tuple to specify them separately.
    """
end

"""
    padplot(plotstr)

Takes a multiline string and pre-pads it so that it shows up as
preformatted code when included in a docstring.
"""
function padplot(plotstr)
    minpad = typemax(Int)
    lines = split(plotstr, "\n")
    # count the minimum number of spaces preceeding any line
    for line in lines
        line == "" && continue
        minpad = min(minpad, length(match(r"^ *", line).match))
    end

    pad = " " ^ (4-minpad)
    join((string(pad, line) for line in lines), "\n")
end

include("winplots.jl")

"""
    makewindow(winfunc::Function, n::Integer, padding::Integer, zerophase::Bool)

Generate a discrete window vector of the given length with `padding` zeros.
`winfunc` should be a function giving the window value in the range of
[-0.5, 0.5]. The function is assumed to be 0 outside of this range.

$zerophase_docs

Example:

```julia
function hanning(n::Integer; padding::Integer=0, zerophase::Bool=false)
    makewindow(n, padding, zerophase) do x
        0.5*(1+cos(2pi*x))
    end
end
```
"""
function makewindow(winfunc::Function, n::Integer, padding::Integer, zerophase::Bool)
    if n < 0
        throw(ArgumentError("`n` must be nonnegative"))
    end
    if padding < 0
        throw(ArgumentError("`padding` must be nonnegative"))
    end
    win = zeros(n+padding)
    if n == 1
        win[1] = winfunc(0.0)
    elseif zerophase
        # note that the endpoint of the window gets set in both lines. In the
        # unpadded case this will set the same index (which shouldn't make a
        # difference if the window is symmetric), but it's necessary for when
        # there's padding, which ends up in the center of the vector length
        # n÷2+1
        win[1:n÷2+1] .= winfunc.(range(0.0, stop=(n÷2)/n, length=n÷2+1))
        # length n÷2
        win[end-n÷2+1:end] .= winfunc.(range(-(n÷2)/n, stop=-1/n, length=n÷2))
    else
        win[1:n] .= winfunc.(range(-0.5, stop=0.5, length=n))
    end

    win
end

#
# Window functions
#

# equations rendered into Unicode with https://arthursonzogni.com/Diagon/

"""
$rect_winplot

    rect(n::Integer; padding::Integer=0, zerophase::Bool=false)
    rect(dims; padding=0, zerophase=false)

Rectangular window of length `n`, padded with `padding` zeros. This window is 1
within the window, and 0 outside of it.

$(twoD_docs())

$zerophase_docs
"""
function rect(n::Integer; padding::Integer=0, zerophase::Bool=false)
    makewindow(n, padding, zerophase) do x
        1.0
    end
end


"""
$hanning_winplot

    hanning(n::Integer; padding::Integer=0, zerophase::Bool=false)
    hanning(dims; padding=0, zerophase=false)

Hanning window of length `n` with `padding` zeros. The Hanning (or Hann) window
is a raised-cosine window that reaches zero at the endpoints.

The window is defined by sampling the continuous function:

           1 + cos(2πx)
    w(x) = ──────────── = cos²(πx)
                2

in the range `[-0.5, 0.5]`

The `hanning` window satisfies the Constant Overlap-Add (COLA) property with an
hop of 0.5, which means that adding together a sequence of delayed windows with
50% overlap will result in a constant signal. This is useful when synthesizing
a signal from a number of overlapping frames (each with a roughly rectangular
window), to eliminate windowing amplitude modulation.

Note that the `hanning` window is the `cosine` window squared.

$(twoD_docs())

$zerophase_docs
"""
function hanning(n::Integer; padding::Integer=0, zerophase::Bool=false)
    makewindow(n, padding, zerophase) do x
        0.5*(1+cos(2pi*x))
    end
end

"""
$hamming_winplot

    hamming(n::Integer; padding::Integer=0, zerophase::Bool=false)
    hamming(dims; padding=0, zerophase=false)

Hamming window of length `n` with `padding` zeros. The Hamming window does not
reach zero at the endpoints and so has a shallower frequency roll-off when
compared to the Hanning window, but is designed to cancel the first side-lobe.

The window is defined by sampling the continuous function:

    w(x) = 0.54 + 0.46*cos(2pi*x)

in the range `[-0.5, 0.5]`

$(twoD_docs())

$zerophase_docs
"""
function hamming(n::Integer; padding::Integer=0, zerophase::Bool=false)
    makewindow(n, padding, zerophase) do x
        0.54 + 0.46*cos(2pi*x)
    end
end

"""
$tukey_winplot

    tukey(n::Integer, α::Real; padding::Integer=0, zerophase::Bool=false)
    tukey(dims, α; padding=0, zerophase=false)

Tukey window of length `n` with `padding` zeros. The Tukey window has a flat top
and reaches zero at the endpoints, with a sinusoidal transition area
parameterized by `α`. For `α == 0`, the window is equivalent to a rectangular
window. For `α == 1`, the window is a Hann window.

The window is defined by sampling the continuous function:

           ⎛              ⎛    ⎛    1 - α⎞⎞
           ⎜      1 + cos ⎜2πα ⎜x + ─────⎟⎟             1 - α
           ⎜              ⎝    ⎝      2  ⎠⎠         x ≤ ─────
           ⎜      ─────────────────────────               2
           ⎜                  2
           ⎜
    w(x) = ⎜      1                                 -α/2 < x ≤ α/2
           ⎜
           ⎜              ⎛    ⎛    1 - α⎞⎞
           ⎜      1 + cos ⎜2πα ⎜x - ─────⎟⎟             1 - α
           ⎜              ⎝    ⎝      2  ⎠⎠         x > ─────
           ⎜      ─────────────────────────               2
           ⎝                  2

in the range `[-0.5, 0.5]`

$(twoD_docs("α"))

$zerophase_docs
"""
function tukey(n::Integer, α::Real; padding::Integer=0, zerophase::Bool=false)
    # check that α is reasonable
    !(0 <= α <= 1) && error("α must be in the range 0 <= α <= 1.")

    # if α is less than machine precision, call it zero and return the
    # rectangular window for this length.  if we don't short circuit this
    # here, it will blow up below.
    abs(α) <= eps() && return rect(n; padding=padding, zerophase=zerophase)

    makewindow(n, padding, zerophase) do x
        if x <= -(1-α)/2
            0.5*(1 + cos(2pi/α*(x+(1-α)/2)))
        elseif x <= (1-α)/2
            1.0
        else
            0.5*(1 + cos(2pi/α*(x-(1-α)/2)))
        end
    end
end

"""
$cosine_winplot

    cosine(n::Integer; padding::Integer=0, zerophase::Bool=false)
    cosine(dims; padding=0, zerophase=false)

Cosine window of length `n` with `padding` zeros. The cosine window is the first
lobe of a cosine function (with the zero crossings at +/- π as endpoints). Also
called the sine window.

The window is defined by sampling the continuous function:

    w(x) = cos(πx)

in the range `[-0.5, 0.5]`

Note that the cosine window is the square root of the `hanning` window, so it is
sometimes used when you are applying the window twice, such as the analysis and
synthesis steps of an STFT.

$(twoD_docs())

$zerophase_docs
"""
function cosine(n::Integer; padding::Integer=0, zerophase::Bool=false)
    makewindow(n, padding, zerophase) do x
        cos(pi*x)
    end
end

"""
$lanczos_winplot

    lanczos(n::Integer; padding::Integer=0, zerophase::Bool=false)
    lanczos(dims; padding=0, zerophase=false)

Lanczos window of length `n` with `padding` zeros. The Lanczos window is the
main lobe of a `sinc` function.

The window is defined by sampling the continuous function:

                      sin(2πx)
    w(x) = sinc(2x) = ────────
                         2πx

in the range `[-0.5, 0.5]`

$(twoD_docs())

$zerophase_docs
"""
function lanczos(n::Integer; padding::Integer=0, zerophase::Bool=false)
    makewindow(n, padding, zerophase) do x
        sinc(2x)
    end
end

"""
$triang_winplot

    triang(n::Integer; padding::Integer=0, zerophase::Bool=false)
    triang(dims; padding=0, zerophase=false)

Triangular window of length `n` with `padding` zeros. The Triangular window does
not reach zero at the endpoints. For odd `n` the `triang` window is the center
`n` points of an `n+2`-point [`bartlett`](@ref) window (i.e. the samples just
outside the window would be zero). For even `n` the window slope is the same as
the `n-1` window but delayed by a half sample so the zero points would be 1/2
sample past the ends of the window.

The window is defined by sampling the continuous function:

            ⎛    2(n-1)
            ⎜1 - ────── abs(x)     n is even
            ⎜       n
    w(x) =  ⎜
            ⎜    2(n-1)
            ⎜1 - ────── abs(x)     n is odd
            ⎝     n+1

in the range `[-0.5, 0.5]`.

$(twoD_docs())

$zerophase_docs

When `zerophase` is `true` substitute `n+1` for `n` in the above window
expressions.
"""
function triang(n::Integer; padding::Integer=0, zerophase::Bool=false)
    # for the purpose of calculating the slope of the window, consider `n` to be
    # 1 larger to compensate for the fact that `zerophase` gives a periodic
    # window
    m = zerophase ? n+1 : n
    scale = iseven(m) ? 2(m-1)/m : 2(m-1)/(m+1)
    makewindow(n, padding, zerophase) do x
        1 - scale*abs(x)
    end
end

"""
$bartlett_winplot

    bartlett(n::Integer; padding::Integer=0, zerophase::Bool=false)
    bartlett(dims; padding=0, zerophase=false)

Bartlett window of length `n`. The Bartlett window is a triangular window that
reaches 0 at the endpoints. This is equivalent to convolving two rectangular
windows of length `(n-1)/2` and adding the zero endpoints. See `triang` for a
window that does not reach zero at the endpoints.

The window is defined by sampling the continuous function:

    1 - abs(2x)

in the range `[-0.5, 0.5]`

$(twoD_docs())

$zerophase_docs
"""
function bartlett(n::Integer; padding::Integer=0, zerophase::Bool=false)
    makewindow(n, padding, zerophase) do x
        1 - abs(2x)
    end
end

"""
$gaussian_winplot

    gaussian(n::Integer, σ::Real; padding::Integer=0, zerophase::Bool=false)
    gaussian(dims, σ; padding=0, zerophase=false)

Gives an n-sample gaussian window defined by sampling the function:

            ⎛        2⎞
            ⎜-1   ⎛x⎞ ⎟
            ⎜── ⋅ ⎜─⎟ ⎟
            ⎝ 2   ⎝σ⎠ ⎠
    w(x) = e

in the range `[-0.5,0.5]`. This means that for `σ=0.5` the endpoints of the
window will correspond to 1 standard deviation away from the center.

$(twoD_docs("σ"))

$zerophase_docs
"""
function gaussian(n::Integer, σ::Real; padding::Integer=0, zerophase::Bool=false)
    σ > 0.0 || error("σ must be positive")
    makewindow(n, padding, zerophase) do x
        exp(-0.5*(x/σ)^2)
    end
end

"""
$bartlett_hann_winplot

    bartlett_hann(n::Integer; padding::Integer=0, zerophase::Bool=false)
    bartlett_hann(dims; padding=0, zerophase=false)

Bartlett-Hann window of length `n` with `padding` zeros. The Bartlett-Hann
window is a weighted sum of the Bartlett and Hann windows.

The window is defined by sampling the continuous function:

    w(x) = 0.62 - 0.48*abs(x) + 0.38*cos(2π*x)

in the range `[-0.5, 0.5]`

$(twoD_docs())

$zerophase_docs
"""
function bartlett_hann(n::Integer; padding::Integer=0, zerophase::Bool=false)
    a0, a1, a2 = 0.62, 0.48, 0.38
    makewindow(n, padding, zerophase) do x
        a0 - a1*abs(x) + a2*cos(2pi*x)
    end
end

"""
$blackman_winplot

    blackman(n::Integer; padding::Integer=0, zerophase::Bool=false)
    blackman(dims; padding=0, zerophase=false)

Approximates the "Exact" Blackman window. This is the generalized Blackman
window with α = 0.16.

The window is defined by sampling the continuous function:

    w(x) = 0.42 + 0.5*cos(2π*x) + 0.08*cos(4π*x)

in the range `[-0.5, 0.5]`

$(twoD_docs())

$zerophase_docs
"""
function blackman(n::Integer; padding::Integer=0, zerophase::Bool=false)
    a0, a1, a2 = 0.42, 0.5, 0.08
    makewindow(n, padding, zerophase) do x
        a0 + a2*cospi(4*x) + a1*cospi(2*x)
    end
end

"""
$kaiser_winplot

    kaiser(n::Integer, α::Real; padding::Integer=0, zerophase::Bool=false)
    kaiser(dims, α; padding=0, zerophase=false)

Kaiser window of length `n` parameterized by `α`. The Kaiser window approximates
the DPSS window (given by `dpss`), using a simplified definition relying on a
Bessel function. Larger values for `α` give a wider main lobe but have lower
sidelobes. Typically `α` is set around 3.

The window is defined by sampling the continuous function:


             ⎛  ⎛   _________⎞⎞
             ⎜  ⎜  ╱        2⎟⎟
    w(x) = I₀⎝πα⎝╲╱ 1 - (2x) ⎠⎠
           ────────────────────
                   I₀(πα)

in the range `[-0.5, 0.5]`

Where I₀(⋅) is the zeroth-order modified Bessel function of the first kind.

$(twoD_docs("α"))

$zerophase_docs
"""
function kaiser(n::Integer, α::Real; padding::Integer=0, zerophase::Bool=false)
    pf = 1.0/besseli(0,pi*α)
    makewindow(n, padding, zerophase) do x
        pf*besseli(0, pi*α*(sqrt(1 - (2x)^2)))
    end
end

# Discrete prolate spheroid sequences (Slepian tapers)
#
# See Gruenbacher, D. M., & Hummels, D. R. (1994). A simple algorithm
# for generating discrete prolate spheroidal sequences. IEEE
# Transactions on Signal Processing, 42(11), 3276-3278.
"""
$dpss_winplot

    dpss(n::Integer, nw::Real, ntapers::Integer=iceil(2*nw)-1;
         padding::Integer=0, zerophase::Bool=false)

The first `ntapers` discrete prolate spheroid sequences (Slepian
tapers) as an `n` × `ntapers` matrix. The signs of the tapers
follow the convention that the first element of the skew-symmetric
(odd) tapers is positive. The time-bandwidth product is given by
`nw`.

The DPSS window maximizes the energy concentration in the main lobe.

$zerophase_docs
"""
function dpss(n::Integer, nw::Real, ntapers::Integer=ceil(Int, 2*nw)-1;
              padding::Integer=0, zerophase::Bool=false)
    if isodd(n) && zerophase
        throw(ArgumentError("`dpss` does not currently support odd-length zerophase windows"))
    end
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
    eigvec = eigen!(mat, n-ntapers+1:n).vectors
    rv = reverse(eigvec, dims=2)::Matrix{Float64}

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
        # TODO: figure out how to handle odd-length zerophase windows. See the
        # tests for one approach, but upsampling by 2 doesn't work here becasue
        # the overall amplitude varies with `n`
        rv = ifftshift(rv, 1)
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

# convert the given argument to a 2-Tuple so we can distribute it to the
# input windows when making a 2D window
argdup(arg::Tuple) = arg
argdup(arg::Real) = (arg, arg)
const IntegerOr2 = Union{Tuple{<:Integer, <:Integer}, Integer}
const RealOr2 = Union{Tuple{<:Real, <:Real}, Real}
const BoolOr2 = Union{Tuple{Bool, Bool}, Bool}

for func in (:rect, :hanning, :hamming, :cosine, :lanczos,
             :triang, :bartlett, :bartlett_hann, :blackman)
    @eval begin
        function $func(dims::Tuple; padding::IntegerOr2=0,
                                    zerophase::BoolOr2=false)
            length(dims) == 2 || throw(ArgumentError("`dims` must be length 2"))
            paddings = argdup(padding)
            zerophases = argdup(zerophase)
            w1 = $func(dims[1], padding=paddings[1], zerophase=zerophases[1])
            w2 = $func(dims[2], padding=paddings[2], zerophase=zerophases[2])
            w1 * w2'
        end
    end
end

for func in (:tukey, :gaussian, :kaiser)
    @eval begin
        function $func(dims::Tuple, arg::RealOr2;
                       padding::IntegerOr2=0, zerophase::BoolOr2=false)
            length(dims) == 2 || throw(ArgumentError("`dims` must be length 2"))
            args = argdup(arg)
            paddings = argdup(padding)
            zerophases = argdup(zerophase)
            w1 = $func(dims[1], args[1]; padding=paddings[1], zerophase=zerophases[1])
            w2 = $func(dims[2], args[2]; padding=paddings[2], zerophase=zerophases[2])
            w1 * w2'
        end
    end
end

end # end module definition
