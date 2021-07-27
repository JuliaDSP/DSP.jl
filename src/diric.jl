export diric

"""
    kernel = diric(Ω::Real, n::Int)

Dirichlet kernel, also known as periodic sinc function,
where `n` should be a positive integer.
Returns `sin(Ω * n/2) / (n * sin(Ω / 2))` typically,
but `±1` when `Ω` is a multiple of 2π.

In the usual case where 'n' is odd, the output is equivalent to
``1/n \\sum_{k=-(n-1)/2}^{(n-1)/2} e^{i k Ω}``,
which is the discrete-time Fourier transform (DTFT)
of a `n`-point moving average filter.

When `n` is odd or even, the function is 2π or 4π periodic, respectively.
The formula for general `n` is
`diric(Ω,n) = ``e^{-i (n-1) Ω/2}/n \\sum_{k=0}^{n-1} e^{i k Ω}``.

When `Ω` is an `AbstractFloat` (e.g., `Float32` or `BigFloat`),
the return type matches that of `Ω`.  Otherwise the return type is `Float64`.

As of 2021-03-19, the Wikipedia definition has different factors.
The definition here is consistent with scipy and other software frameworks.

This implementation treats inputs near multiples of 2π fairly carefully.
If the denominator `sin(Ω/2)` is close to 0, then this function uses
an accurate 2nd-order Taylor expansion rather than simple division. 

# Examples

```jldoctest
julia> round.(diric.((-2:0.5:2)*π, 5), digits=9)'
1×9 adjoint(::Vector{Float64}) with eltype Float64:
 1.0  -0.2  0.2  -0.2  1.0  -0.2  0.2  -0.2  1.0

julia> diric(0, 4)
1.0
```
"""
function diric(Ω::T, n::Int) where T <: AbstractFloat
    n > 0 || throw(ArgumentError("n=$n is non-positive"))
    sign = 1
    if isodd(n)
        Ω = rem2pi(Ω, RoundNearest) # [-π,π)
    else
        Ω = (VERSION < v"1.1") ?
            BigFloat(Ω, 2*precision(T)) :
            BigFloat(Ω ; precision = 2*precision(T))
        Ω = 2 * rem2pi(abs(Ω)/2, RoundNearest) # [-2π,2π)
        if Ω > π
            sign = -1
            Ω -= 2π # (-π,π]
        end
    end

    denom = sin(Ω / 2)
    atol = eps(T)
    if abs(denom) ≤ atol # denom ≈ 0 ?
#=
        if iseven(n) && abs(Ω) ≈ 2T(π)
            # 2nd-order Taylor expansion near ±2π for even n
            return abs2(abs(Ω)-2T(π)) * (n*n - 1) / 24 - one(T)
        end
=#
    #   return one(T) - abs2(Ω) * (n*n - 1) / 24 # 2nd-order Taylor near 0
        return sign * one(T)
    end

    return T(sign * sin(Ω * n/2) / (n * denom)) # typical case
end

# handle non AbstractFloat types (e.g., Int, Rational)
diric(Ω::Real, n::Int) = diric(Float64(Ω), n::Int)

# handle π exactly
function diric(::Irrational{:π}, n::Int)
    n > 0 || throw(ArgumentError("n=$n is non-positive"))
    return 1 // n
end
