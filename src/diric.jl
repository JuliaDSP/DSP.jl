export diric

"""
    kernel = diric(Ω::Real, n::Integer)

Dirichlet kernel, also known as periodic sinc function,
where `n` should be a positive integer.

Typically returns ``\\frac{\\sin(Ω \\cdot n/2)}{n * \\sin(Ω / 2)}``,
but returns ``±1`` when `Ω` is a multiple of 2π.

In the usual case where `n` is odd, the output is equivalent to
``\\frac{1}{n} \\sum_{k=-(n-1)/2}^{(n-1)/2} e^{i k Ω}``,
which is the discrete-time Fourier transform (DTFT)
of a `n`-point moving average filter.

When `n` is odd or even, the function is 2π or 4π periodic, respectively.
The formula for general `n` is
```math
\\mathrm{diric}(Ω,n) = \\frac{e^{-i (n-1) Ω/2}}{n} \\sum_{k=0}^{n-1} e^{i k Ω}
```
which is a real and symmetric function of `Ω`.

As of 2021-03-19, the Wikipedia definition has different factors.
The definition here is consistent with scipy and other software frameworks.

# Examples

```jldoctest
julia> round.(diric.((-2:0.5:2)*π, 5), digits=9)'
1×9 adjoint(::Vector{Float64}) with eltype Float64:
 1.0  -0.2  0.2  -0.2  1.0  -0.2  0.2  -0.2  1.0

julia> diric(0, 4)
1.0
```
"""
function diric(Ω::T, n::Integer) where T <: AbstractFloat
    n > 0 || throw(ArgumentError("n=$n not positive"))
    sign = one(T)
    if isodd(n)
        Ω = rem2pi(Ω, RoundNearest) # [-π,π)
    else # wrap to interval [-π,π) to improve precision near ±2π
        Ω = 2 * rem2pi(Ω/2, RoundNearest) # [-2π,2π)
        if Ω > π # [π,2π)
            sign = -one(T)
            Ω -= T(2π) # [-π,0)
        elseif Ω < -π # [-2π,-π)
            sign = -one(T)
            Ω += T(2π) # [0,π)
        end
    end

    denom = sin(Ω / 2)
    atol = eps(T)
    if abs(denom) ≤ atol # denom ≈ 0 ?
        return sign
    end

    return sign * sin(Ω * n/2) / (n * denom) # typical case
end

# handle non AbstractFloat types (e.g., Int, Rational, π)
diric(Ω::Real, n::Integer) = diric(float(Ω), n)
