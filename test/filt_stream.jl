using DSP
using Base.Test

# Disable time and printf macros when not running interactivly ( for travis )
if !isinteractive()
    macro time( ex )
        quote
            $(esc(ex))
        end
    end

    macro printf( args... )
    end
end

function Base.isapprox( x1::Vector, x2::Vector )
    Nx1 = length( x1 )
    Nx2 = length( x2 )

    if Nx1 != Nx2
        @printf("x1 & x2 are different lengths vectors\n")
        return false
    end

    for i = 1:Nx1
        if !isapprox( x1[i], x2[i] )
            @printf( "Something went wrong at index %d\n", i )
            return false
        end
    end

    return true
end




#==============================================================================#
#               ____ _ _  _ ____ _    ____    ____ ____ ___ ____               #
#               [__  | |\ | | __ |    |___    |__/ |__|  |  |___               #
#               ___] | | \| |__] |___ |___    |  \ |  |  |  |___               #
#==============================================================================#

function test_singlerate( h, x )
    xLen       = length( x )
    hLen       = length( h )
    pivotPoint = min( rand(50:150, 1)[1], ifloor( xLen/4 ))
    x1         = x[ 1 : pivotPoint ]
    x2         = x[ pivotPoint+1 : end ]

    @printf( "\n\n" )
    @printf( "____ _ _  _ ____ _    ____    ____ ____ ___ ____\n" )
    @printf( "[__  | |\\ | | __ |    |___    |__/ |__|  |  |___\n" )
    @printf( "___] | | \\| |__] |___ |___    |  \\ |  |  |  |___\n" )
    @printf( "\nTesting single-rate fitering, h is %s, x is %s. xLen = %d, hLen = %d", string(eltype(h)), string(eltype(x)), xLen, hLen )

    @printf( "\n\tBase.filt\n\t\t")
    @time baseResult = Base.filt( h, 1.0, x )

    if method_exists( DSP.firfilt, ( typeof(h), typeof(x) ))
        @printf( "\n\tDSP.firfilt\n\t\t")
        @time dspResult = DSP.firfilt( h, x )
    end

    @printf( "\n\tDSP.filt( h, x, 1//1 )\n\t\t" )
    @time statelesResult = DSP.filt( h, x )

    @printf( "\n\tDSP.filt. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = DSP.FIRFilter( h, 1//1 )
    @time begin
        y1 = DSP.filt( self, x1 )
        y2 = DSP.filt( self, x2 )
    end
    statefulResult = [ y1, y2 ]

    @printf( "\n\tDSP.filt filt. Piecewise for first %d inputs\n\t\t", length( x1 ) )
    DSP.reset( self )
    @time begin
        for i in 1:length(x1)
            y1[i] = DSP.filt( self, x1[i:i] )[1]
        end
        y2 = DSP.filt( self, x2 )
    end
    piecewiseResult = [ y1, y2 ]


     if isapprox( baseResult, statelesResult ) && isapprox( baseResult, statefulResult ) && isapprox( baseResult, piecewiseResult )
         return true
     end

     display( [ baseResult statelesResult statefulResult piecewiseResult ] )

     return false

end




#==============================================================================#
#                      ___  ____ ____ _ _  _ ____ ___ ____                     #
#                      |  \ |___ |    | |\/| |__|  |  |___                     #
#                      |__/ |___ |___ | |  | |  |  |  |___                     #
#==============================================================================#

function test_decimation( h, x, decimation )
    xLen       = length( x )
    hLen       = length( h )
    pivotPoint = min( rand(50:150, 1)[1], ifloor( xLen/4 ))
    x1         = x[ 1 : pivotPoint ]
    x2         = x[ pivotPoint+1 : end ]

    @printf( "\n\n" )
    @printf( "___  ____ ____ _ _  _ ____ ___ _ ____ _  _ \n" )
    @printf( "|  \\ |___ |    | |\\/| |__|  |  | |  | |\\ | \n" )
    @printf( "|__/ |___ |___ | |  | |  |  |  | |__| | \\| \n" )
    @printf( "\nTesting decimation. h::%s, x::%s. xLen = %d, hLen = %d, decimation = %d", string(typeof(h)), string(typeof(h)), xLen, hLen, decimation )

    @printf( "\n\tNaive decimation\n\t\t")
    @time begin naiveResult = DSP.naivefilt( h, x, 1//decimation )
    end

    @printf( "\n\tDSP.filt( h, x, 1//%d)\n\t\t", decimation )
    @time statelesResult = DSP.filt( h, x, 1//decimation )

    @printf( "\n\tDSP.filt decimation. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = DSP.FIRFilter( h, 1//decimation )
    @time begin
        y1 = DSP.filt( self, x1 )
        y2 = DSP.filt( self, x2 )
    end
    statefulResult = [ y1, y2 ]

    @printf( "\n\tDSP.filt decimation. Piecewise for first %d inputs.\n\t\t", length( x1 ) )
    DSP.reset( self )
    y1 = similar( x, 0 )
    @time begin
        for i in 1:length(x1)
            append!( y1, DSP.filt( self, x1[i:i] ) )
        end
        y2 = DSP.filt( self, x2 )
    end
    piecewiseResult = [ y1, y2 ]

    if isapprox( naiveResult, statelesResult ) && isapprox( naiveResult, statefulResult ) && isapprox( naiveResult, piecewiseResult )
        return true
    end

    display( [ naiveResult statefulResult piecewiseResult ] )
    return false

end




#==============================================================================#
#               _ _  _ ___ ____ ____ ___  _    ____ ____ ___ ____              #
#               | |\ |  |  |___ |__/ |__] |    |  | |__|  |  |___              #
#               | | \|  |  |___ |  \ |    |___ |__| |  |  |  |___              #
#==============================================================================#

function test_interpolation( h, x, interpolation )
    xLen       = length( x )
    hLen       = length( h )
    pivotPoint = min( rand(50:150, 1)[1], ifloor( xLen/4 ))
    x1         = x[ 1 : pivotPoint ]
    x2         = x[ pivotPoint+1 : end ]

    @printf( "\n\n" )
    @printf( "_ _  _ ___ ____ ____ ___  _    ____ ____ ___ _ ____ _  _ \n" )
    @printf( "| |\\ |  |  |___ |__/ |__] |    |  | |__|  |  | |  | |\\ | \n" )
    @printf( "| | \\|  |  |___ |  \\ |    |___ |__| |  |  |  | |__| | \\| \n" )
    @printf( "\nTesting interpolation, h::%s, x::%s. xLen = %d, hLen = %d, interpolation = %d", typeof(h), typeof(x), xLen, hLen, interpolation )

    @printf( "\n\tNaive interpolation with Base.filt\n\t\t")
    @time begin
        xZeroStuffed = zeros( eltype(x), xLen * interpolation )
        for n = 0:xLen-1;
            xZeroStuffed[ n*interpolation+1 ] = x[ n+1 ]
        end
        baseResult = Base.filt( h, one(eltype(h)), xZeroStuffed )
    end

    if method_exists( DSP.firfilt, ( typeof(h), typeof(x) ))
        @printf( "\n\tNaive interpolation with DSP.firfilt\n\t\t")
        @time begin
            xZeroStuffed = zeros( eltype(x), xLen * interpolation )
            for n = 0:xLen-1;
                xZeroStuffed[ n*interpolation+1 ] = x[ n+1 ]
            end
            dspResult = DSP.firfilt( h, xZeroStuffed )
        end
    end

    @printf( "\n\tDSP.filt( h, x, %d//1 )\n\t\t", interpolation )
    @time statelesResult = DSP.filt( h, x, interpolation//1 )

    @printf( "\n\tDSP.filt interpolation. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = DSP.FIRFilter( h, interpolation//1 )
    @time begin
        y1 = DSP.filt( self, x1 )
        y2 = DSP.filt( self, x2 )
    end
    statefulResult = [ y1, y2 ]

    @printf( "\n\tDSP.filt interpolation. Piecewise for first %d inputs\n\t\t", length( x1 ) )
    DSP.reset( self )
    y1 = similar( x, 0 )
    @time begin
        for i in 1:length(x1)
            append!( y1, DSP.filt( self, x1[i:i] ) )
        end
        y2 = DSP.filt( self, x2 )
    end
    piecewiseResult = [ y1, y2 ]

    if isapprox( baseResult, statelesResult ) && isapprox( baseResult, statefulResult ) && isapprox( piecewiseResult, baseResult )
        return true
    end

    display( [ [1:length(baseResult)] baseResult statefulResult piecewiseResult ] )

    return false
end




#==============================================================================#
#           ____ ____ ___     ____ ____ ____ ____ _  _ ___  _    ____          #
#           |__/ |__|  |      |__/ |___ [__  |__| |\/| |__] |    |___          #
#           |  \ |  |  |  .   |  \ |___ ___] |  | |  | |    |___ |___          #
#==============================================================================#

function test_rational( h, x, ratio )
    xLen       = length( x )
    hLen       = length( h )
    pivotPoint = min( rand(50:150, 1)[1], ifloor( xLen/4 ))
    x1         = x[ 1 : pivotPoint ]
    x2         = x[ pivotPoint+1 : end ]
    upfactor   = num( ratio )
    downfactor = den( ratio )
    resultType = promote_type( eltype(h), eltype(x) )

    @printf( "\n\n" )
    @printf( "      ____ ____ ___ _ ____ _  _ ____ _    \n" )
    @printf( "      |__/ |__|  |  | |  | |\\ | |__| |    \n" )
    @printf( "      |  \\ |  |  |  | |__| | \\| |  | |___ \n" )
    @printf( "                                          \n" )
    @printf( "____ ____ ____ ____ _  _ ___  _    _ _  _ ____\n" )
    @printf( "|__/ |___ [__  |__| |\\/| |__] |    | |\\ | | __\n" )
    @printf( "|  \\ |___ ___] |  | |  | |    |___ | | \\| |__]\n" )
    @printf( "\n\nTesting rational resampling, h::%s, x::%s. xLen = %d, hLen = %d, ratio = %d//%d", string(typeof(h)), string(typeof(x)), xLen, hLen, upfactor, downfactor )

    @printf( "\n\tNaive rational resampling\n\t\t")
    @time naiveResult = DSP.naivefilt( h, x, ratio )

    @printf( "\n\tDSP.filt( h, x, %d//%d )\n\t\t", upfactor, downfactor )
    @time statelesResult = DSP.filt( h, x, ratio )

    @printf( "\n\tDSP.filt rational resampling. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = DSP.FIRFilter( h, ratio )
    @time begin
        s1 = DSP.filt( self, x1 )
        s2 = DSP.filt( self, x2 )
    end
    statefulResult = [ s1, s2 ]

    @printf( "\n\tDSP.filt rational. Piecewise for all %d inputs\n\t\t", length( x ) )
    self = DSP.FIRFilter( h, ratio )
    y1 = similar( x, 0 )
    @time begin
        for i in 1:length(x)
            append!( y1, DSP.filt( self, x[i:i] ) )
        end
    end
    piecewiseResult = y1

    if isapprox( naiveResult, statelesResult ) && isapprox( naiveResult, statefulResult ) && isapprox( naiveResult, piecewiseResult )
        return true
    end

    display( [  naiveResult statefulResult piecewiseResult ] )

    return false
end




#==============================================================================#
#        ____ ____ ___      ____ ____ ____ ____ _  _ ___  _    _ _  _ ____     #
#        |__| |__/ |__]     |__/ |___ [__  |__| |\/| |__] |    | |\ | | __     #
#        |  | |  \ |__] .   |  \ |___ ___] |  | |  | |    |___ | | \| |__]     #
#==============================================================================#

function test_arbitrary( Th, x, resampleRate, numFilters )
    cutoffFreq      = 0.45
    transitionWidth = 0.05
    h               = DSP.firdes( Lowpass(cutoffFreq, fs=numFilters), transitionWidth/numFilters ) .* numFilters
    h               = convert( Vector{Th}, h )

    @printf( "____ ____ ___      ____ ____ ____ ____ _  _ ___  _    _ _  _ ____\n" )
    @printf( "|__| |__/ |__]     |__/ |___ [__  |__| |\\/| |__] |    | |\\ | | __\n" )
    @printf( "|  | |  \\ |__] .   |  \\ |___ ___] |  | |  | |    |___ | | \\| |__]\n" )
    @printf( "\n\nh::%s, x::%s, rate = %f, N𝜙 = %d, xLen = %d, ", string(typeof(h)), string(typeof(x)), resampleRate, numFilters, length(x) )

    @printf( "\n\tNaive arbitrary resampling\n\t\t" )
    @time naiveResult = DSP.naivefilt( h, x, resampleRate, numFilters )

    @printf( "\n\tStateless arbitrary resampling\n\t\t" )
    @time statelessResult = DSP.filt( h, x, resampleRate, numFilters )

    @printf( "\n\tPiecewise arbitrary resampling\n\t\t" )
    self           = DSP.FIRFilter( h, resampleRate, numFilters )
    piecwiseResult = eltype(x)[]
    sizehint( piecwiseResult, iceil( length(x)*resampleRate ) )
    @time for i in 1:length( x )
        thisY = filt( self, x[i:i] )
        append!( piecwiseResult, thisY )
    end

    commonLen = min( length(naiveResult), length(statelessResult), length(piecwiseResult) )

    resize!( naiveResult, commonLen )
    resize!( statelessResult, commonLen )
    resize!( piecwiseResult, commonLen )

    if isapprox( naiveResult, statelessResult ) && isapprox( naiveResult, piecwiseResult )
        return true
    end

    display( [  [1:commonLen] naiveResult statelessResult  piecwiseResult abs(naiveResult.-statelessResult) abs(naiveResult.-piecwiseResult) ] )


    display( [ [1:commonLen-1] diff(naiveResult) diff(statelessResult) diff(piecwiseResult) ] )
    return false
end




#==============================================================================#
#                  ____ _  _ _  _    ___ ____ ____ ___ ____                    #
#                  |__/ |  | |\ |     |  |___ [__   |  [__                     #
#                  |  \ |__| | \|     |  |___ ___]  |  ___]                    #
#==============================================================================#

function test_all()
    for interpolation in sort([1, unique(rand(2:32,8))] ),
            decimation in sort([1, unique(rand(2:32,8))] ),
                Th in [Float32, Float64],
                    Tx in [Float32, Float64, Complex64, Complex128]

        h     = rand(Th, rand(16:128,1)[1] )
        xLen  = int(rand( 200:300, 1 )[1])
        xLen  = xLen-mod( xLen, decimation )
        x     = rand( Tx, xLen )
        ratio = interpolation//decimation

        if ratio == 1
            @test test_singlerate( h, x )
        end

        if decimation != 1
            @test test_decimation( h, x, decimation )
        end

        if interpolation != 1
            @test test_interpolation( h, x, interpolation )
        end


        if num(ratio) == interpolation && den(ratio) == decimation && num(ratio) != 1 && den(ratio) != 1
            @test test_rational( h, x, ratio )
            if Tx in [ Float32, Complex64 ]
                @test test_arbitrary( Th, x, float64(ratio)+rand(), 32 )
            end
        end
    end
end

# function test_nextphase()
#     for interpolation in 1:8
#         for decimation in 1:8
#             ratio           = interpolation//decimation
#             interpolation   = num(ratio)
#             decimation      = den(ratio)
#             x               = repmat( [1:interpolation], decimation )
#             reference = [ x[n] for n = 1:decimation:length( x ) ]
#             result = [ 1 ]
#             for i in 2:interpolation
#                 append!( result, [ DSP.nextphase( result[end], ratio ) ] )
#             end
#             @test isapprox( reference, result )
#         end
#     end
# end
#
# test_nextphase()
test_all()
