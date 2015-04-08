using DSP
using Base.Test

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
    @timeifinteractive naiveResult = Base.filt(h, 1.0, x)

    if method_exists(DSP.firfilt, (typeof(h), typeof(x)))
        @printfifinteractive( "\n\tDSP.firfilt\n\t\t")
        @timeifinteractive dspResult = DSP.firfilt( h, x )
    end

    @printfifinteractive( "\n\tDSP.filt( h, x, 1//1 )\n\t\t" )
    @timeifinteractive statelesResult = DSP.filt( h, x )
    @test_approx_eq naiveResult statelesResult

    @printfifinteractive( "\n\tDSP.filt. length(x1) = %d, length(x2) = %d\n\t\t", length(x1), length(x2) )
    self = DSP.FIRFilter(h, 1//1)
    @timeifinteractive begin
        y1 = DSP.filt(self, x1)
        y2 = DSP.filt(self, x2)
    end
    statefulResult = [y1; y2]
    @test_approx_eq naiveResult statefulResult

    @printfifinteractive( "\n\tDSP.filt filt. Piecewise for first %d inputs\n\t\t", length(x1) )
    DSP.reset(self)
    @timeifinteractive begin
        for i in 1:length(x1)
            y1[i] = DSP.filt(self, x1[i:i])[1]
        end
        y2 = DSP.filt(self, x2)
    end
    piecewiseResult = [y1; y2]
    @test_approx_eq naiveResult piecewiseResult
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
    @timeifinteractive naiveResult = DSP.naivefilt(h, x, 1//decimation)

    @printfifinteractive( "\n\tDSP.filt( h, x, 1//%d)\n\t\t", decimation )
    @timeifinteractive statelesResult = DSP.filt(h, x, 1//decimation)
    @test_approx_eq naiveResult statelesResult

    @printfifinteractive( "\n\tDSP.filt decimation. length(x1) = %d, length(x2) = %d\n\t\t", length(x1), length(x2) )
    self = DSP.FIRFilter(h, 1//decimation)
    @timeifinteractive begin
        y1 = DSP.filt(self, x1)
        y2 = DSP.filt(self, x2)
    end
    statefulResult = [y1; y2]
    @test_approx_eq naiveResult statefulResult

    @printfifinteractive( "\n\tDSP.filt decimation. Piecewise for first %d inputs.\n\t\t", length(x1) )
    DSP.reset( self )
    y1 = similar( x, 0 )
    @timeifinteractive begin
        for i in 1:length(x1)
            append!(y1, DSP.filt(self, x1[i:i]))
        end
        y2 = DSP.filt(self, x2)
    end
    piecewiseResult = [y1; y2]
    # @test_approx_eq naiveResult, piecewiseResult
    @test all(map(isapprox, naiveResult, piecewiseResult))
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
        naiveResult = Base.filt(h, one(eltype(h)), xZeroStuffed)
    end

    @printfifinteractive( "\n\tDSP.filt( h, x, %d//1 )\n\t\t", interpolation )
    @timeifinteractive statelesResult = DSP.filt( h, x, interpolation//1 )
    @test_approx_eq naiveResult statelesResult

    @printfifinteractive( "\n\tDSP.filt interpolation. length(x1) = %d, length(x2) = %d\n\t\t", length(x1), length(x2) )
    self = DSP.FIRFilter( h, interpolation//1 )
    @timeifinteractive begin
        y1 = DSP.filt(self, x1)
        y2 = DSP.filt(self, x2)
    end
    statefulResult = [y1; y2]
    @test_approx_eq naiveResult statefulResult

    @printfifinteractive( "\n\tDSP.filt interpolation. Piecewise for first %d inputs\n\t\t", length(x1) )
    DSP.reset(self)
    y1 = similar(x, 0)
    @timeifinteractive begin
        for i in 1:length(x1)
            append!(y1, DSP.filt(self, x1[i:i]))
        end
        y2 = DSP.filt(self, x2)
    end
    piecewiseResult = [y1; y2]
    # @test_approx_eq naiveResult piecewiseResult
    @test all(map(isapprox, naiveResult, piecewiseResult))
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
    upfactor   = num(ratio)
    downfactor = den(ratio)
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
    @timeifinteractive naiveResult = DSP.naivefilt(h, x, ratio)

    @printfifinteractive( "\n\tDSP.filt( h, x, %d//%d )\n\t\t", upfactor, downfactor )
    @timeifinteractive statelesResult = DSP.filt(h, x, ratio)
    @test_approx_eq naiveResult statelesResult

    @printfifinteractive( "\n\tDSP.filt rational resampling. length(x1) = %d, length(x2) = %d\n\t\t", length(x1), length(x2) )
    self = DSP.FIRFilter(h, ratio)
    @timeifinteractive begin
        s1 = DSP.filt(self, x1)
        s2 = DSP.filt(self, x2)
    end
    statefulResult = [s1; s2]
    @test_approx_eq naiveResult statefulResult

    @printfifinteractive( "\n\tDSP.filt rational. Piecewise for all %d inputs\n\t\t", length( x ) )
    self = DSP.FIRFilter(h, ratio)
    y1 = similar(x, 0)
    @timeifinteractive begin
        for i in 1:length(x)
            append!(y1, DSP.filt(self, x[i:i]))
        end
    end
    piecewiseResult = y1
    @test_approx_eq naiveResult piecewiseResult

    return false
end


#
# Arbitrary resampling
#

function test_arbitrary(Th, x, resampleRate, numFilters)
    cutoffFreq      = 0.45
    transitionWidth = 0.05
    h               = digitalfilter(Lowpass(cutoffFreq, fs=numFilters), WindowFIR(transitionwidth=transitionWidth/numFilters)) .* numFilters
    h               = convert(Vector{Th}, h)

    @printfifinteractive( "____ ____ ___      ____ ____ ____ ____ _  _ ___  _    _ _  _ ____\n" )
    @printfifinteractive( "|__| |__/ |__]     |__/ |___ [__  |__| |\\/| |__] |    | |\\ | | __\n" )
    @printfifinteractive( "|  | |  \\ |__] .   |  \\ |___ ___] |  | |  | |    |___ | | \\| |__]\n" )
    @printfifinteractive( "\n\nh::%s, x::%s, rate = %f, N𝜙 = %d, xLen = %d, ", string(typeof(h)), string(typeof(x)), resampleRate, numFilters, length(x) )

    @printfifinteractive( "\n\tNaive arbitrary resampling\n\t\t" )
    @timeifinteractive naiveResult = DSP.naivefilt(h, x, resampleRate, numFilters)

    @printfifinteractive( "\n\tStateless arbitrary resampling\n\t\t" )
    @timeifinteractive statelessResult = DSP.filt(h, x, resampleRate, numFilters)

    @printfifinteractive( "\n\tPiecewise arbitrary resampling\n\t\t" )
    self           = DSP.FIRFilter(h, resampleRate, numFilters)
    piecwiseResult = eltype(x)[]
    sizehint!(piecwiseResult, ceil(Int, length(x)*resampleRate))
    @timeifinteractive for i in 1:length(x)
        thisY = filt(self, x[i:i])
        append!(piecwiseResult, thisY)
    end

    commonLen = min(length(naiveResult), length(statelessResult), length(piecwiseResult))

    resize!(naiveResult, commonLen)
    resize!(statelessResult, commonLen)
    resize!(piecwiseResult, commonLen)

    @test_approx_eq naiveResult statelessResult
    @test_approx_eq naiveResult piecwiseResult
end

#
# Run the tests
#

function test_all()
    for interpolation in sort([1; unique(rand(2:32,8))]),
            decimation in sort([1; unique(rand(2:32,8))]),
                Th in [Float32, Float64],
                    Tx in [Float32, Float64, Complex64, Complex128]

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


        if num(ratio) == interpolation && den(ratio) == decimation && num(ratio) != 1 && den(ratio) != 1
            test_rational(h, x, ratio)
            if Tx in [Float32, Complex64]
                test_arbitrary(Th, x, convert(Float64, ratio)+rand(), 32)
            end
        end
    end
end

function test_nextphase()
    for interpolation in 1:8
        for decimation in 1:8
            ratio           = interpolation//decimation
            interpolation   = num(ratio)
            decimation      = den(ratio)
            x               = repmat([1:interpolation;], decimation)
            reference = [x[n] for n = 1:decimation:length(x)]
            result = [1]
            for i in 2:interpolation
                append!(result, [DSP.Filters.nextphase(result[end], ratio)])
            end
            @test_approx_eq reference result
        end
    end
end

test_nextphase()
test_all()
