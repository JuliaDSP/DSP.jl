export diric

"""
    kernel = diric(Ω::Real, n::Int)

Dirichlet kernel, also known as periodic sinc function,
where `n` should be a positive integer.
Returns `sin(Ω * n/2) / (n * sin(Ω / 2))` typically,
but `±1` when `Ω` is a multiple of 2π.

In the usual case where 'n' is odd, the output is equivalent to
`1/n \\sum_{k=-(n+1)/2}^{(n+1)/2} e^{i k Ω}`
which is the same as the discrete-time Fourier transform (DTFT)
of a `n`-point moving average filter.

When `n` is odd or even, the function is 2π or 4π periodic, respectively.

When `Ω` is an `AbstractFloat` (e.g., `Float32` or `BigFloat`),
the return type matches that of `Ω`.  Otherwise the return type is `Float64`.

As of 2021-03-19, the Wikipedia definition has different factors.
The definition here is consistent with scipy and other software frameworks.
This implementation treats inputs near multiples of 2π fairly carefully.
"""
function diric(Ω::T, n::Int) where T <: AbstractFloat
    n <= 0 && throw(ArgumentError("n=$n is non-positive"))
    denom = sin(Ω / 2)
    if isapprox(1 + denom, 1) # denom ≈ 0 ?
        if isodd(n)
            ω = mod2pi(Ω + π) - π # (-π,π)
        else
            π2 = T == BigFloat ? 2T(π) : T(2π) # for type stability
            π4 = T == BigFloat ? 4T(π) : T(4π)
            ω = mod(Ω + π2, π4) - π2 # (-2π,2π)
            if abs(ω) ≈ π2
                # 2nd-order Taylor expansion near ±2π
                return (abs(ω)-π2)^2 * (n^2 - 1) / 24 - one(T)
            end
        end
        return one(T) - ω^2 * (n^2 - 1) / 24 # 2nd-order Taylor near 0
    end

    return sin(Ω * n/2) / (n * denom) # typical case
end

# handle non AbstractFloat types (e.g., Int, Irrational)
diric(Ω::Real, n::Int) = diric(Float64(Ω), n::Int)
