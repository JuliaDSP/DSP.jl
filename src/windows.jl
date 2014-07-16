module Windows

export rect, hanning, hamming, tukey, cosine, lanczos, 
       triang, bartlett, gaussian, bartlett_hann, blackman, 
       kaiser, dpss
# 
# Window functions
#

# Rectangular window function of length N.
function rect(n::Integer)
    ones(n)
end

# Hanning window of length N.
function hanning(n::Integer)
    [0.5*(1 - cos(2*pi*k/(n-1))) for k=0:(n-1)]
end

# Hamming window of length N.
function hamming(n::Integer)
    [0.54 - 0.46*cos(2*pi*k/(n-1)) for k=0:(n-1)]
end

# Tukey window of length N, parameterized by alpha.  For
# alpha = 0, the window is equivalent to a rectangular window.
# For alpha = 1, the window is a hann window.
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

# Cosine window of length N.  Also called the sine window for obvious reasons.
function cosine(n::Integer)
    [sin(pi*k/(n-1)) for k=0:(n-1)]
end

# Lanczos window of length N.
function lanczos(n::Integer)
    [sinc(2*k/(n-1) - 1) for k=0:(n-1)]
end

# triangular window of length N.
function triang(n::Integer)
    [1 - abs((k - (n-1)/2))/(n/2) for k=0:(n-1)]
end

# bartlett window of length N.
function bartlett(n::Integer)
    [2/(n-1)*((n-1)/2 - abs(k - (n-1)/2)) for k=0:(n-1)]
end

# gaussian window of length N parameterized by the standard deviation
# sigma
function gaussian(n::Integer, sigma::Real)
    if !(0 < sigma <= 0.5)
        error("sigma must be greater than 0 and less than or equal to 0.5.")
    end
    [exp(-0.5*((k-(n-1)/2)/(sigma*(n-1/2)))^2) for k=0:(n-1)]
end

# bartlett-hann window of length n
function bartlett_hann(n::Integer)
    a0, a1, a2 = 0.62, 0.48, 0.38
    t = 2*pi/(n-1)
    [a0 - a1*abs(k/(n-1) - 0.5) - a2*cos(t*k) for k=0:(n-1)]
end

# "exact" blackman window, alpha=0.16
function blackman(n::Integer)
    a0, a1, a2 = 0.42, 0.5, 0.08
    t = 2*pi/(n-1)
    [a0 - a1*cos(t*k) + a2*cos(t*k*2) for k=0:(n-1)]
end

# kaiser window parameterized by alpha
function kaiser(n::Integer, alpha::Real)
    pf = 1.0/besseli(0,pi*alpha)
    [pf*besseli(0, pi*alpha*(sqrt(1 - (2*k/(n-1) - 1)^2))) for k=0:(n-1)]
end

# Discrete prolate spheroid sequences (Slepian tapers)
#
# See Gruenbacher, D. M., & Hummels, D. R. (1994). A simple algorithm
# for generating discrete prolate spheroidal sequences. IEEE
# Transactions on Signal Processing, 42(11), 3276-3278.
function dpss(n::Int, nw::Real, ntapers::Int=iceil(2*nw)-1)
    # Construct symmetric tridiagonal matrix
    mat = SymTridiagonal([cospi(2*nw/n)*abs2((n - 1)/2 - i) for i=0:(n-1)],
                         [0.5.*(i*n - abs2(i)) for i=1:(n-1)])

    # Get tapers
    v = fliplr(eigfact!(mat, n-ntapers+1:n)[:vectors]::Matrix{Float64})

    # Slepian's convention; taper starts with a positive element
    sgn = ones(size(v, 2))
    for i = 2:2:size(v, 2)
        sgn[i] = sign(v[1, i])
    end
    scale!(v, sgn)
end

end # end module definition
