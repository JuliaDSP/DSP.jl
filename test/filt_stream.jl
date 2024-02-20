using DSP, Test

# Naive rational resampler
function naivefilt(h::Vector, x::Vector{T}, resamplerate::Union{Integer, Rational}=1) where T

    upfactor     = numerator(resamplerate)
    downfactor   = denominator(resamplerate)
    xLen         = length(x)
    xZeroStuffed = zeros(T, xLen * upfactor)

    for n in 0:xLen-1
        xZeroStuffed[n*upfactor+1] = x[n+1]
    end

    y = filt(h, one(T), xZeroStuffed)
    y = [y[n] for n = 1:downfactor:length(y)]
end


# Naive arbitrary resampler
function naivefilt(h::Vector, x::Vector, resamplerate::AbstractFloat, numfilters::Integer=32)
    xLen          = length(x)
    xInterpolated = naivefilt(h, x, numfilters//1)
    xLen          = length(xInterpolated)
    yLen          = ceil(Int, xLen * resamplerate)
    y             = similar(x, yLen)
    yIdx          = 1
    xIdx          = 1
    α             = 0.0
    (δ, ϕStride)  = modf(numfilters/resamplerate)
    ϕStride       = convert(Int, ϕStride)

    while xIdx < xLen
        yLower  = xInterpolated[xIdx]
        yUpper  = xInterpolated[xIdx+1]
        y[yIdx] = yLower + α*(yUpper - yLower)
        yIdx   += 1
        α      += δ
        xIdx   += floor(Int, α) + ϕStride
        α       = mod(α, 1.0)
    end

    resize!(y, yIdx-1)

    return y
end

# Disable time and printf macros when not running interactively ( for travis )
if isinteractive()
    macro timeifinteractive(ex)
        :(@time $(esc(ex)))
    end

    using Printf: @printf

    macro printfifinteractive(s, args...)
        :(@printf($s, $(map(esc, args)...)))
    end
else
    macro timeifinteractive(ex)
        esc(ex)
    end

    macro printfifinteractive(args...)
    end
end


#
# Single rate filtering
#

function test_singlerate(h, x)
    xLen       = length(x)
    hLen       = length(h)
    pivotPoint = min(rand(50:150), div(xLen, 4))
    x1         = x[1:pivotPoint]
    x2         = x[pivotPoint+1:end]

    @printfifinteractive( "\n\n" )
    @printfifinteractive( "____ _ _  _ ____ _    ____    ____ ____ ___ ____\n" )
    @printfifinteractive( "[__  | |\\ | | __ |    |___    |__/ |__|  |  |___\n" )
    @printfifinteractive( "___] | | \\| |__] |___ |___    |  \\ |  |  |  |___\n" )
    @printfifinteractive( "\nTesting single-rate filtering, h::%s, x::%s. xLen = %d, hLen = %d\n", typeof(h), typeof(x), xLen, hLen )

    @printfifinteractive( "\n\tfilt\n\t\t")
    @timeifinteractive naiveResult = filt(h, 1.0, x)

    @printfifinteractive( "\n\tfilt( h, x, 1//1 )\n\t\t" )
    @timeifinteractive statelessResult = filt( h, x )
    @test naiveResult ≈ statelessResult

    @printfifinteractive( "\n\tfilt. length(x1) = %d, length(x2) = %d\n\t\t", length(x1), length(x2) )
    myfilt = FIRFilter(h, 1//1)
    @timeifinteractive begin
        y1 = filt(myfilt, x1)
        y2 = filt(myfilt, x2)
    end
    statefulResult = [y1; y2]
    @test naiveResult ≈ statefulResult

    @printfifinteractive( "\n\tfilt filt. Piecewise for first %d inputs\n\t\t", length(x1) )
    reset!(myfilt)
    @timeifinteractive begin
        for i in 1:length(x1)
            y1[i] = filt(myfilt, x1[i:i])[1]
        end
        y2 = filt(myfilt, x2)
    end
    piecewiseResult = [y1; y2]
    @test naiveResult ≈ piecewiseResult

    reset!(myfilt)
    @test inputlength(myfilt, length(piecewiseResult)) == xLen
end


#
# Decimation
#

function test_decimation(h, x, decimation)
    xLen       = length(x)
    hLen       = length(h)
    pivotPoint = min(rand(50:150), div(xLen, 4))
    x1         = x[1:pivotPoint]
    x2         = x[pivotPoint+1:end]

    @printfifinteractive( "\n\n" )
    @printfifinteractive( "___  ____ ____ _ _  _ ____ ___ _ ____ _  _ \n" )
    @printfifinteractive( "|  \\ |___ |    | |\\/| |__|  |  | |  | |\\ | \n" )
    @printfifinteractive( "|__/ |___ |___ | |  | |  |  |  | |__| | \\| \n" )
    @printfifinteractive( "\nTesting decimation. h::%s, x::%s. xLen = %d, hLen = %d, decimation = %d\n", typeof(h), typeof(x), xLen, hLen, decimation )

    @printfifinteractive( "\n\tNaive decimation\n\t\t")
    @timeifinteractive naiveResult = naivefilt(h, x, 1//decimation)

    @printfifinteractive( "\n\tfilt( h, x, 1//%d)\n\t\t", decimation )
    @timeifinteractive statelessResult = filt(h, x, 1//decimation)
    @test naiveResult ≈ statelessResult

    @printfifinteractive( "\n\tfilt decimation. length(x1) = %d, length(x2) = %d\n\t\t", length(x1), length(x2) )
    myfilt = FIRFilter(h, 1//decimation)
    @timeifinteractive begin
        y1 = filt(myfilt, x1)
        y2 = filt(myfilt, x2)
    end
    statefulResult = [y1; y2]
    @test naiveResult ≈ statefulResult

    @printfifinteractive( "\n\tfilt decimation. Piecewise for first %d inputs.\n\t\t", length(x1) )
    reset!(myfilt)
    y1 = similar( x, 0 )
    @timeifinteractive begin
        for i in 1:length(x1)
            append!(y1, filt(myfilt, x1[i:i]))
        end
        y2 = filt(myfilt, x2)
    end
    piecewiseResult = [y1; y2]
    @test ≈(naiveResult, piecewiseResult, atol=sqrt(eps(real(one(eltype(x))))))

    reset!(myfilt)
    @test inputlength(myfilt, length(piecewiseResult)) == xLen
end


#
# Interpolation
#

function test_interpolation(h::AbstractVector{T}, x::AbstractVector{V}, interpolation) where {T,V}
    xLen = length(x)
    hLen       = length(h)
    pivotPoint = min(rand(50:150), div(xLen, 4))
    x1         = x[1:pivotPoint]
    x2         = x[pivotPoint+1:end]

    @printfifinteractive( "\n\n" )
    @printfifinteractive( "_ _  _ ___ ____ ____ ___  _    ____ ____ ___ _ ____ _  _ \n" )
    @printfifinteractive( "| |\\ |  |  |___ |__/ |__] |    |  | |__|  |  | |  | |\\ | \n" )
    @printfifinteractive( "| | \\|  |  |___ |  \\ |    |___ |__| |  |  |  | |__| | \\| \n" )
    @printfifinteractive( "\nTesting interpolation, h::%s, x::%s. xLen = %d, hLen = %d, interpolation = %d\n", typeof(h), typeof(x), xLen, hLen, interpolation )

    @printfifinteractive( "\n\tNaive interpolation with filt\n\t\t")
    @timeifinteractive begin
        xZeroStuffed = zeros(V, xLen * interpolation)
        for n = 0:xLen-1
            xZeroStuffed[n*interpolation+1] = x[n+1]
        end
        naiveResult = filt(h, one(T), xZeroStuffed)
    end

    @printfifinteractive( "\n\tfilt( h, x, %d//1 )\n\t\t", interpolation )
    @timeifinteractive statelessResult = filt( h, x, interpolation//1 )
    @test naiveResult ≈ statelessResult

    @printfifinteractive( "\n\tfilt interpolation. length(x1) = %d, length(x2) = %d\n\t\t", length(x1), length(x2) )
    myfilt = FIRFilter( h, interpolation//1 )
    @timeifinteractive begin
        y1 = filt(myfilt, x1)
        y2 = filt(myfilt, x2)
    end
    statefulResult = [y1; y2]
    @test naiveResult ≈ statefulResult

    reset!(myfilt)
    @test inputlength(myfilt, length(statefulResult)) == xLen

    @printfifinteractive( "\n\tfilt interpolation. Piecewise for first %d inputs\n\t\t", length(x1) )
    reset!(myfilt)
    y1 = similar(x, 0)
    @timeifinteractive begin
        for i in 1:length(x1)
            append!(y1, filt(myfilt, x1[i:i]))
        end
        y2 = filt(myfilt, x2)
    end
    piecewiseResult = [y1; y2]
    @test ≈(naiveResult, piecewiseResult, atol=sqrt(eps(real(one(V)))))
end


#
# Rational resampling
#

function test_rational(h, x, ratio)
    xLen       = length(x)
    hLen       = length(h)
    pivotPoint = min(rand(50:150), div(xLen, 4))
    x1         = x[1:pivotPoint]
    x2         = x[pivotPoint+1:end]
    upfactor   = numerator(ratio)
    downfactor = denominator(ratio)
    resultType = promote_type(eltype(h), eltype(x))

    @printfifinteractive( "\n\n" )
    @printfifinteractive( "      ____ ____ ___ _ ____ _  _ ____ _    \n" )
    @printfifinteractive( "      |__/ |__|  |  | |  | |\\ | |__| |    \n" )
    @printfifinteractive( "      |  \\ |  |  |  | |__| | \\| |  | |___ \n" )
    @printfifinteractive( "                                          \n" )
    @printfifinteractive( "____ ____ ____ ____ _  _ ___  _    _ _  _ ____\n" )
    @printfifinteractive( "|__/ |___ [__  |__| |\\/| |__] |    | |\\ | | __\n" )
    @printfifinteractive( "|  \\ |___ ___] |  | |  | |    |___ | | \\| |__]\n" )
    @printfifinteractive( "\n\nTesting rational resampling, h::%s, x::%s. xLen = %d, hLen = %d, ratio = %d//%d\n", typeof(h), typeof(x), xLen, hLen, upfactor, downfactor )

    @printfifinteractive( "\n\tNaive rational resampling\n\t\t")
    @timeifinteractive naiveResult = naivefilt(h, x, ratio)

    @printfifinteractive( "\n\tfilt( h, x, %d//%d )\n\t\t", upfactor, downfactor )
    @timeifinteractive statelessResult = filt(h, x, ratio)
    @test naiveResult ≈ statelessResult

    @printfifinteractive( "\n\tfilt rational resampling. length(x1) = %d, length(x2) = %d\n\t\t", length(x1), length(x2) )
    myfilt = FIRFilter(h, ratio)
    @timeifinteractive begin
        s1 = filt(myfilt, x1)
        s2 = filt(myfilt, x2)
    end
    statefulResult = [s1; s2]
    @test naiveResult ≈ statefulResult

    @printfifinteractive( "\n\tfilt rational. Piecewise for all %d inputs\n\t\t", length( x ) )
    reset!(myfilt)
    y1 = similar(x, 0)
    @timeifinteractive begin
        for i in 1:length(x)
            append!(y1, filt(myfilt, x[i:i]))
        end
    end
    piecewiseResult = y1
    @test naiveResult ≈ piecewiseResult

    reset!(myfilt)
    @test inputlength(myfilt, length(piecewiseResult)) == xLen
end


#
# Arbitrary resampling
#

function test_arbitrary(Th, x, resampleRate, numFilters)
    cutoffFreq      = 0.45
    transitionWidth = 0.05
    h               = digitalfilter(Lowpass(cutoffFreq, fs=numFilters), FIRWindow(transitionwidth=transitionWidth/numFilters)) .* numFilters
    h               = convert(Vector{Th}, h)
    myfilt          = FIRFilter(h, resampleRate, numFilters)
    xLen            = length(x)

    @printfifinteractive("\n\n")
    @printfifinteractive( "____ ____ ___      ____ ____ ____ ____ _  _ ___  _    _ _  _ ____\n" )
    @printfifinteractive( "|__| |__/ |__]     |__/ |___ [__  |__| |\\/| |__] |    | |\\ | | __\n" )
    @printfifinteractive( "|  | |  \\ |__] .   |  \\ |___ ___] |  | |  | |    |___ | | \\| |__]\n" )
    @printfifinteractive( "\nh::%s, x::%s, rate = %f, Nϕ = %d, xLen = %d\n", typeof(h), typeof(x), resampleRate, numFilters, xLen )

    @printfifinteractive( "\n\tNaive arbitrary resampling\n\t\t" )
    @timeifinteractive naiveResult = naivefilt(h, x, resampleRate, numFilters)

    @printfifinteractive( "\n\tStateless arbitrary resampling\n\t\t" )
    @timeifinteractive statelessResult = filt(h, x, resampleRate, numFilters)

    @printfifinteractive( "\n\tStateful arbitrary resampling\n\t\t" )
    @timeifinteractive statefulResult = filt(myfilt, x)

    # reset!(myfilt)
    # TODO: figure out why this fails
    # @test inputlength(myfilt, length(statefulResult)) == xLen

    @printfifinteractive( "\n\tPiecewise arbitrary resampling\n\t\t" )
    reset!(myfilt)
    piecewiseResult = similar(x, 0)
    sizehint!(piecewiseResult, ceil(Int, length(x)*resampleRate))
    @timeifinteractive for i in 1:length(x)
        thisY = filt(myfilt, x[i:i])
        append!(piecewiseResult, thisY)
    end

    commonLen = minimum(length.((naiveResult, statelessResult, statefulResult, piecewiseResult)))
    resize!(naiveResult, commonLen)
    resize!(statelessResult, commonLen)
    resize!(statefulResult, commonLen)
    resize!(piecewiseResult, commonLen)

    @test naiveResult ≈ statelessResult
    @test naiveResult ≈ statefulResult
    @test naiveResult ≈ piecewiseResult
end

#
# Run the tests
#

@testset "interp.=$interpolation, dec.=$decimation, Th=$Th, Tx=$Tx" for interpolation in [1, 5, 14, 23],
            decimation in [1, 9, 17, 21],
                Th in [Float32, Float64],
                    Tx in [Float32, Float64, ComplexF32, ComplexF64]

    h     = rand(Th, rand(16:128))
    xLen  = rand(200:300)
    xLen  = xLen-mod(xLen, decimation)
    x     = rand(Tx, xLen)
    ratio = interpolation//decimation
    if ratio == 1
        test_singlerate(h, x)
    end
    if decimation != 1
        test_decimation(h, x, decimation)
    end
    if interpolation != 1
        test_interpolation(h, x, interpolation)
    end
    if numerator(ratio) == interpolation && denominator(ratio) == decimation && numerator(ratio) != 1 && denominator(ratio) != 1
        test_rational(h, x, ratio)
        if Tx in [Float32, ComplexF32]
            test_arbitrary(Th, x, convert(Float64, ratio)+rand(), 32)
        end
    end
    if decimation == 1
        test_rational(h, x, interpolation)
    end
end
