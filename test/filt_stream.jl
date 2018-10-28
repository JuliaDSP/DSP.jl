using DSP, Compat, Compat.Test

# Naive rational resampler
function naivefilt(h::Vector, x::Vector, resamplerate::Union{Integer, Rational}=1)

    upfactor     = numerator(resamplerate)
    downfactor   = denominator(resamplerate)
    xLen         = length(x)
    xZeroStuffed = zeros(eltype(x), length(x) * upfactor)

    for n in 0:length(x)-1
        xZeroStuffed[n*upfactor+1] = x[n+1]
    end

    y = filt(h, one(eltype(x)), xZeroStuffed)
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

# Disable time and printf macros when not running interactivly ( for travis )
if isinteractive()
    macro timeifinteractive(ex)
        :(@time $(esc(ex)))
    end

    macro printfifinteractive(args...)
        :(@printf($(map(esc, args)...)))
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
    @printfifinteractive( "\nTesting single-rate fitering, h is %s, x is %s. xLen = %d, hLen = %d", string(eltype(h)), string(eltype(x)), xLen, hLen )

    @printfifinteractive( "\n\tBase.filt\n\t\t")
    @timeifinteractive naiveResult = filt(h, 1.0, x)

    @printfifinteractive( "\n\tDSP.filt( h, x, 1//1 )\n\t\t" )
    @timeifinteractive statelesResult = DSP.filt( h, x )
    @test naiveResult ≈ statelesResult

    @printfifinteractive( "\n\tDSP.filt. length(x1) = %d, length(x2) = %d\n\t\t", length(x1), length(x2) )
    myfilt = DSP.FIRFilter(h, 1//1)
    @timeifinteractive begin
        y1 = DSP.filt(myfilt, x1)
        y2 = DSP.filt(myfilt, x2)
    end
    statefulResult = [y1; y2]
    @test naiveResult ≈ statefulResult

    @printfifinteractive( "\n\tDSP.filt filt. Piecewise for first %d inputs\n\t\t", length(x1) )
    DSP.reset!(myfilt)
    @timeifinteractive begin
        for i in 1:length(x1)
            y1[i] = DSP.filt(myfilt, x1[i:i])[1]
        end
        y2 = DSP.filt(myfilt, x2)
    end
    piecewiseResult = [y1; y2]
    @test naiveResult ≈ piecewiseResult

    DSP.reset!(myfilt)
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
    @printfifinteractive( "\nTesting decimation. h::%s, x::%s. xLen = %d, hLen = %d, decimation = %d", string(typeof(h)), string(typeof(h)), xLen, hLen, decimation )

    @printfifinteractive( "\n\tNaive decimation\n\t\t")
    @timeifinteractive naiveResult = naivefilt(h, x, 1//decimation)

    @printfifinteractive( "\n\tDSP.filt( h, x, 1//%d)\n\t\t", decimation )
    @timeifinteractive statelesResult = DSP.filt(h, x, 1//decimation)
    @test naiveResult ≈ statelesResult

    @printfifinteractive( "\n\tDSP.filt decimation. length(x1) = %d, length(x2) = %d\n\t\t", length(x1), length(x2) )
    myfilt = DSP.FIRFilter(h, 1//decimation)
    @timeifinteractive begin
        y1 = DSP.filt(myfilt, x1)
        y2 = DSP.filt(myfilt, x2)
    end
    statefulResult = [y1; y2]
    @test naiveResult ≈ statefulResult

    @printfifinteractive( "\n\tDSP.filt decimation. Piecewise for first %d inputs.\n\t\t", length(x1) )
    DSP.reset!(myfilt)
    y1 = similar( x, 0 )
    @timeifinteractive begin
        for i in 1:length(x1)
            append!(y1, DSP.filt(myfilt, x1[i:i]))
        end
        y2 = DSP.filt(myfilt, x2)
    end
    piecewiseResult = [y1; y2]
    @test ≈(naiveResult, piecewiseResult, atol=sqrt(eps(real(one(eltype(x))))))

    DSP.reset!(myfilt)
    @test inputlength(myfilt, length(piecewiseResult)) == xLen
end


#
# Interpolation
#

function test_interpolation(h, x, interpolation)
    xLen       = length(x)
    hLen       = length(h)
    pivotPoint = min(rand(50:150), div(xLen, 4))
    x1         = x[1:pivotPoint]
    x2         = x[pivotPoint+1:end]

    @printfifinteractive( "\n\n" )
    @printfifinteractive( "_ _  _ ___ ____ ____ ___  _    ____ ____ ___ _ ____ _  _ \n" )
    @printfifinteractive( "| |\\ |  |  |___ |__/ |__] |    |  | |__|  |  | |  | |\\ | \n" )
    @printfifinteractive( "| | \\|  |  |___ |  \\ |    |___ |__| |  |  |  | |__| | \\| \n" )
    @printfifinteractive( "\nTesting interpolation, h::%s, x::%s. xLen = %d, hLen = %d, interpolation = %d", typeof(h), typeof(x), xLen, hLen, interpolation )

    @printfifinteractive( "\n\tNaive interpolation with Base.filt\n\t\t")
    @timeifinteractive begin
        xZeroStuffed = zeros(eltype(x), xLen * interpolation)
        for n = 0:xLen-1;
            xZeroStuffed[n*interpolation+1] = x[n+1]
        end
        naiveResult = filt(h, one(eltype(h)), xZeroStuffed)
    end

    @printfifinteractive( "\n\tDSP.filt( h, x, %d//1 )\n\t\t", interpolation )
    @timeifinteractive statelesResult = DSP.filt( h, x, interpolation//1 )
    @test naiveResult ≈ statelesResult

    @printfifinteractive( "\n\tDSP.filt interpolation. length(x1) = %d, length(x2) = %d\n\t\t", length(x1), length(x2) )
    myfilt = DSP.FIRFilter( h, interpolation//1 )
    @timeifinteractive begin
        y1 = DSP.filt(myfilt, x1)
        y2 = DSP.filt(myfilt, x2)
    end
    statefulResult = [y1; y2]
    @test naiveResult ≈ statefulResult

    DSP.reset!(myfilt)
    @test inputlength(myfilt, length(statefulResult)) == xLen

    @printfifinteractive( "\n\tDSP.filt interpolation. Piecewise for first %d inputs\n\t\t", length(x1) )
    DSP.reset!(myfilt)
    y1 = similar(x, 0)
    @timeifinteractive begin
        for i in 1:length(x1)
            append!(y1, DSP.filt(myfilt, x1[i:i]))
        end
        y2 = DSP.filt(myfilt, x2)
    end
    piecewiseResult = [y1; y2]
    @test ≈(naiveResult, piecewiseResult, atol=sqrt(eps(real(one(eltype(x))))))
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
    @printfifinteractive( "\n\nTesting rational resampling, h::%s, x::%s. xLen = %d, hLen = %d, ratio = %d//%d", string(typeof(h)), string(typeof(x)), xLen, hLen, upfactor, downfactor )

    @printfifinteractive( "\n\tNaive rational resampling\n\t\t")
    @timeifinteractive naiveResult = naivefilt(h, x, ratio)

    @printfifinteractive( "\n\tDSP.filt( h, x, %d//%d )\n\t\t", upfactor, downfactor )
    @timeifinteractive statelesResult = DSP.filt(h, x, ratio)
    @test naiveResult ≈ statelesResult

    @printfifinteractive( "\n\tDSP.filt rational resampling. length(x1) = %d, length(x2) = %d\n\t\t", length(x1), length(x2) )
    myfilt = DSP.FIRFilter(h, ratio)
    @timeifinteractive begin
        s1 = DSP.filt(myfilt, x1)
        s2 = DSP.filt(myfilt, x2)
    end
    statefulResult = [s1; s2]
    @test naiveResult ≈ statefulResult

    @printfifinteractive( "\n\tDSP.filt rational. Piecewise for all %d inputs\n\t\t", length( x ) )
    reset!(myfilt)
    y1 = similar(x, 0)
    @timeifinteractive begin
        for i in 1:length(x)
            append!(y1, DSP.filt(myfilt, x[i:i]))
        end
    end
    piecewiseResult = y1
    @test naiveResult ≈ piecewiseResult

    DSP.reset!(myfilt)
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
    myfilt          = DSP.FIRFilter(h, resampleRate, numFilters)
    xLen            = length(x)

    @printfifinteractive( "____ ____ ___      ____ ____ ____ ____ _  _ ___  _    _ _  _ ____\n" )
    @printfifinteractive( "|__| |__/ |__]     |__/ |___ [__  |__| |\\/| |__] |    | |\\ | | __\n" )
    @printfifinteractive( "|  | |  \\ |__] .   |  \\ |___ ___] |  | |  | |    |___ | | \\| |__]\n" )
    @printfifinteractive( "\n\nh::%s, x::%s, rate = %f, Nϕ = %d, xLen = %d, ", string(typeof(h)), string(typeof(x)), resampleRate, numFilters, length(x) )

    @printfifinteractive( "\n\tNaive arbitrary resampling\n\t\t" )
    @timeifinteractive naiveResult = naivefilt(h, x, resampleRate, numFilters)

    @printfifinteractive( "\n\tStateless arbitrary resampling\n\t\t" )
    @timeifinteractive statelessResult = DSP.filt(h, x, resampleRate, numFilters)

    @printfifinteractive( "\n\tStateful arbitrary resampling\n\t\t" )
    @timeifinteractive statefulResult = DSP.filt(myfilt, x)

    # DSP.reset!(myfilt)
    # TODO: figure out why this fails
    # @test inputlength(myfilt, length(statefulResult)) == xLen

    @printfifinteractive( "\n\tPiecewise arbitrary resampling\n\t\t" )
    reset!(myfilt)
    piecwiseResult = eltype(x)[]
    sizehint!(piecwiseResult, ceil(Int, length(x)*resampleRate))
    @timeifinteractive for i in 1:length(x)
        thisY = filt(myfilt, x[i:i])
        append!(piecwiseResult, thisY)
    end

    commonLen = min(length(naiveResult), length(statelessResult), length(statefulResult), length(piecwiseResult))
    resize!(naiveResult, commonLen)
    resize!(statelessResult, commonLen)
    resize!(statefulResult, commonLen)
    resize!(piecwiseResult, commonLen)

    @test naiveResult ≈ statelessResult
    @test naiveResult ≈ statefulResult
    @test naiveResult ≈ piecwiseResult
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
