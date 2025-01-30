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
    return y[1:downfactor:length(y)]
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

function test_singlerate(h::AbstractVector{T}, x::AbstractVector) where T
    xLen       = length(x)
    hLen       = length(h)
    pivotPoint = min(rand(50:150), xLen ÷ 4)
    x1         = x[1:pivotPoint]
    x2         = x[pivotPoint+1:end]

    @printfifinteractive( """\n\n
        ____ _ _  _ ____ _    ____    ____ ____ ___ ____
        [__  | |\\ | | __ |    |___    |__/ |__|  |  |___
        ___] | | \\| |__] |___ |___    |  \\ |  |  |  |___\n\n""" )
    @printfifinteractive( "Testing single-rate filtering, h::%s, x::%s. xLen = %d, hLen = %d\n", typeof(h), typeof(x), xLen, hLen )

    @printfifinteractive( "\n\tfilt\n\t\t")
    @timeifinteractive naiveResult = filt(h, one(T), x)

    @printfifinteractive( "\n\tfilt( h, x, 1//1 )\n\t\t" )
    @timeifinteractive statelessResult = filt(h, x)
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
        one_buffer = similar(x1, 1)
        for (i, p) in enumerate(x1)
            y1[i] = filt(myfilt, setindex!(one_buffer, p, 1))[1]
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

function test_decimation(h, x, decimation::Integer)
    xLen       = length(x)
    hLen       = length(h)
    pivotPoint = min(rand(50:150), xLen ÷ 4)
    x1         = x[1:pivotPoint]
    x2         = x[pivotPoint+1:end]

    @printfifinteractive( """\n\n
        ___  ____ ____ _ _  _ ____ ___ _ ____ _  _
        |  \\ |___ |    | |\\/| |__|  |  | |  | |\\ |
        |__/ |___ |___ | |  | |  |  |  | |__| | \\|\n\n""" )
    @printfifinteractive( "Testing decimation. h::%s, x::%s. xLen = %d, hLen = %d, decimation = %d\n", typeof(h), typeof(x), xLen, hLen, decimation )

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
    y1 = similar(x, 0)
    @timeifinteractive begin
        one_buffer = similar(x1, 1)
        for p in x1
            append!(y1, filt(myfilt, setindex!(one_buffer, p, 1)))
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

function test_interpolation(h::AbstractVector{T}, x::AbstractVector{V}, interpolation::Integer) where {T,V}
    xLen = length(x)
    hLen       = length(h)
    pivotPoint = min(rand(50:150), xLen ÷ 4)
    x1         = x[1:pivotPoint]
    x2         = x[pivotPoint+1:end]

    @printfifinteractive( """\n\n
        _ _  _ ___ ____ ____ ___  ____ _   ____ ___ _ ____ _  _
        | |\\ |  |  |___ |__/ |__] |  | |   |__|  |  | |  | |\\ |
        | | \\|  |  |___ |  \\ |    |__| |__ |  |  |  | |__| | \\|\n\n""")
    @printfifinteractive( "Testing interpolation, h::%s, x::%s. xLen = %d, hLen = %d, interpolation = %d\n", typeof(h), typeof(x), xLen, hLen, interpolation )

    @printfifinteractive( "\n\tNaive interpolation with filt\n\t\t")
    @timeifinteractive begin
        xZeroStuffed = zeros(V, xLen * interpolation)
        for n = 0:xLen-1
            xZeroStuffed[n*interpolation+1] = x[n+1]
        end
        naiveResult = filt(h, one(T), xZeroStuffed)
    end

    @printfifinteractive( "\n\tfilt( h, x, %d//1 )\n\t\t", interpolation )
    @timeifinteractive statelessResult = filt(h, x, interpolation//1)
    @test naiveResult ≈ statelessResult

    @printfifinteractive( "\n\tfilt interpolation. length(x1) = %d, length(x2) = %d\n\t\t", length(x1), length(x2) )
    myfilt = FIRFilter(h, interpolation//1)
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
        one_buffer = similar(x1, 1)
        for p in x1
            append!(y1, filt(myfilt, setindex!(one_buffer, p, 1)))
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
    pivotPoint = min(rand(50:150), xLen ÷ 4)
    x1         = x[1:pivotPoint]
    x2         = x[pivotPoint+1:end]
    upfactor   = numerator(ratio)
    downfactor = denominator(ratio)
    resultType = promote_type(eltype(h), eltype(x))

    @printfifinteractive( """\n\n
              ____ ____ ___ _ ____ _  _ ____ _
              |__/ |__|  |  | |  | |\\ | |__| |
              |  \\ |  |  |  | |__| | \\| |  | |___

        ____ ____ ____ ____ _  _ ___  _    _ _  _ ____
        |__/ |___ [__  |__| |\\/| |__] |    | |\\ | | __
        |  \\ |___ ___] |  | |  | |    |___ | | \\| |__]\n\n""")
    @printfifinteractive( "Testing rational resampling, h::%s, x::%s. xLen = %d, hLen = %d, ratio = %d//%d\n", typeof(h), typeof(x), xLen, hLen, upfactor, downfactor )

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
        one_buffer = similar(x1, 1)
        for p in x
            append!(y1, filt(myfilt, setindex!(one_buffer, p, 1)))
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

function test_arbitrary(Th, x, resampleRate, numFilters::Integer)
    cutoffFreq      = 0.45
    transitionWidth = 0.05
    h               = digitalfilter(Lowpass(cutoffFreq), FIRWindow(transitionwidth=transitionWidth/numFilters); fs=numFilters) .* numFilters
    h               = convert(Vector{Th}, h)
    myfilt          = FIRFilter(h, resampleRate, numFilters)
    xLen            = length(x)

    @printfifinteractive("""\n\n
        ____ ____ ___      ____ ____ ____ ____ _  _ ___  _    _ _  _ ____
        |__| |__/ |__]     |__/ |___ [__  |__| |\\/| |__] |    | |\\ | | __
        |  | |  \\ |__] .   |  \\ |___ ___] |  | |  | |    |___ | | \\| |__]\n\n""")
    @printfifinteractive( "h::%s, x::%s, rate = %f, Nϕ = %d, xLen = %d\n", typeof(h), typeof(x), resampleRate, numFilters, xLen )

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
    one_buffer = similar(x, 1)
    @timeifinteractive for p in x
        thisY = filt(myfilt, setindex!(one_buffer, p, 1))
        append!(piecewiseResult, thisY)
    end

    results = (naiveResult, statelessResult, statefulResult, piecewiseResult)
    commonLen = minimum(length, results)
    foreach(x -> resize!(x, commonLen), results)

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
    if isone(ratio)
        test_singlerate(h, x)
    elseif numerator(ratio) == interpolation
        test_rational(h, x, ratio)
        if Tx in (Float32, ComplexF32)
            test_arbitrary(Th, x, convert(Float64, ratio) + rand(), 32)
        end
    end
    if decimation != 1
        test_decimation(h, x, decimation)
    else
        test_rational(h, x, interpolation)
    end
    if interpolation != 1
        test_interpolation(h, x, interpolation)
    end
end

@test resample(1:2, 3, [zeros(2); 1; zeros(3)]) == [1, 0, 0, 2, 0, 0]
@test resample(1:2, 3//2, [zeros(2); 1; zeros(3)]) == [1, 0, 0]
let H = FIRFilter(2.22)
    setphase!(H, 0.99)
    @test length(filt(H, 1:2)) == 3
end
let H = FIRFilter(122.2)
    setphase!(H, 0.99)
    @test length(filt(H, 1:2)) == 124
end
