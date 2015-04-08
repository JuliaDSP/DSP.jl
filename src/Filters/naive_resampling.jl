# Naive rational resampler
function naivefilt(h::Vector, x::Vector, resamplerate::Rational=1//1)

    upfactor     = num(resamplerate)
    downfactor   = den(resamplerate)
    xLen         = length(x)
    xZeroStuffed = zeros(eltype(x), length(x) * upfactor)

    for n in 0:length(x)-1
        xZeroStuffed[n*upfactor+1] = x[n+1]
    end

    y = Base.filt(h, one(eltype(x)), xZeroStuffed)
    y = [y[n] for n = 1:downfactor:length(y)]
end


# Naive arbitrary resampler
function naivefilt(h::Vector, x::Vector, resamplerate::FloatingPoint, numfilters::Integer=32)
    xLen          = length(x)
    xInterpolated = naivefilt(h, x, numfilters//1)
    xLen          = length(xInterpolated)
    yLen          = ceil(Int, xLen * resamplerate)
    y             = similar(x, yLen)
    yIdx          = 1
    xIdx          = 1
    α             = 0.0
    (δ, 𝜙Stride)  = modf(numfilters/resamplerate)
    𝜙Stride       = convert(Int, 𝜙Stride)

    while xIdx < xLen
        yLower  = xInterpolated[xIdx]
        yUpper  = xInterpolated[xIdx+1]
        y[yIdx] = yLower + α*(yUpper - yLower)
        yIdx   += 1
        α      += δ
        xIdx   += floor(Int, α) + 𝜙Stride
        α       = mod(α, 1.0)
    end

    resize!(y, yIdx-1)

    return y
end
