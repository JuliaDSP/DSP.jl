# This file was formerly a part of Julia. License is MIT: https://julialang.org/license
# TODO: parameterize conv tests
using Test, DSP, OffsetArrays
import DSP: filt, filt!, deconv, conv, xcorr
using DSP: optimalfftfiltlength, unsafe_conv_kern_os!, _conv_kern_fft!, _conv_similar,
    nextfastfft



@testset "filt" begin
    # Filter
    b = [1., 2., 3., 4.]
    x = [1., 1., 0., 1., 1., 0., 0., 0.]
    @test filt(b, 1., x)  == [1., 3., 5., 8., 7., 5., 7., 4.]
    @test filt(b, [1., -0.5], x)  == [1., 3.5, 6.75, 11.375, 12.6875, 11.34375, 12.671875, 10.3359375]
    # With ranges
    @test filt(b, 1., 1.0:10.0) == [1., 4., 10., 20., 30., 40., 50., 60., 70., 80.]
    @test filt(1.:4., 1., 1.0:10.0) == [1., 4., 10., 20., 30., 40., 50., 60., 70., 80.]
    # Across an array is the same as channel-by-channel
    @test filt(b, 1., [x 1.0:8.0]) == [filt(b, 1., x) filt(b, 1., 1.0:8.0)]
    @test filt(b, [1., -0.5], [x 1.0:8.0]) == [filt(b, [1., -0.5], x) filt(b, [1., -0.5], 1.0:8.0)]
    si = zeros(3)
    @test filt(b, 1., [x 1.0:8.0], si) == [filt(b, 1., x, si) filt(b, 1., 1.0:8.0, si)]
    @test si == zeros(3) # Will likely fail if/when arrayviews are implemented
    si = [zeros(3) ones(3)]
    @test filt(b, 1., [x 1.0:8.0], si) == [filt(b, 1., x, zeros(3)) filt(b, 1., 1.0:8.0, ones(3))]
    # With initial conditions: a lowpass 5-pole butterworth filter with W_n = 0.25,
    # and a stable initial filter condition matched to the initial value
    b = [0.003279216306360201,0.016396081531801006,0.03279216306360201,0.03279216306360201,0.016396081531801006,0.003279216306360201]
    a = [1.0,-2.4744161749781606,2.8110063119115782,-1.703772240915465,0.5444326948885326,-0.07231566910295834]
    si = [0.9967207836936347,-1.4940914728163142,1.2841226760316475,-0.4524417279474106,0.07559488540931815]
    @test filt(b, a, ones(10), si) ≈ ones(10) # Shouldn't affect DC offset

    @test_throws ArgumentError filt!([1, 2], [1], [1], [1])
end

@testset "conv" begin

    @test optimalfftfiltlength(1, 3) == 1 # Should be at least the first input

    @testset "conv-1D" begin
        # Convolution
        a = [1, 2, 1, 2]
        b = [1, 2, 3]
        expectation = [1, 4, 8, 10, 7, 6]
        im_expectation = [1, 3, 6, 6, 5, 3]
        a32 = convert(Array{Int32}, a)
        @test conv(a, b) == expectation
        @test conv(a32, b) == expectation
        fa = convert(Array{Float64}, a)
        f32a = convert(Array{Float32}, a)
        fb = convert(Array{Float64}, b)
        fexp = convert(Array{Float64}, expectation)
        im_fexp = convert(Array{Float64}, im_expectation)
        @test conv(fa, fb) ≈ fexp
        @test conv(complex.(fa, 1.), complex.(fb)) ≈ complex.(fexp, im_fexp)

        @test conv(fa, b) ≈ fexp
        @test conv(f32a, b) ≈ fexp
        @test conv(fb, a) ≈ fexp

        offset_arr = OffsetArray{Int}(undef, -1:2)
        offset_arr[:] = a
        @test conv(offset_arr, 1:3) == OffsetVector(expectation, 0:5)
        # Issue #352
        @test conv([1//2, 1//3, 1//4], [1, 2]) ≈ [1//2, 4//3, 11//12, 1//2]
        # Non-numerical arrays should not be convolved
        @test_throws MethodError conv([sin], [cos])
    end


    @testset "conv-2D" begin
        a =[1 2 1;
            2 3 1;
            1 2 1]

        b = [3 2;
             0 1]

        expectation = [3 8 7 2;
                       6 14 11 3;
                       3 10 10 3;
                       0 1 2 1]
        im_expectation = [3 5 5 2;
                          3 6 6 3;
                          3 6 6 3;
                          0 1 1 1]

        # Integers
        # Real Integers
        a32 = convert(Array{Int32}, a)
        @test conv(a, b) == expectation
        @test conv(a32, b) == expectation
        # Floats
        fa = convert(Array{Float64}, a)
        f32a = convert(Array{Float32}, a)
        fb = convert(Array{Float64}, b)
        fexp = convert(Array{Float64}, expectation)
        im_fexp = convert(Array{Float64}, im_expectation)
        # Real
        @test conv(fa, fb) == fexp
        # Complex
        @test conv(complex.(fa, 1), complex.(fb)) == complex.(fexp,
                                                               im_fexp)
        @test conv(fa, b) ≈ fexp
        @test conv(f32a, b) ≈ fexp
        @test conv(fb, a) ≈ fexp

        offset_arr = OffsetArray{Int}(undef, -1:1, -1:1)
        offset_arr[:] = a
        @test conv(offset_arr, b) == OffsetArray(expectation, 0:3, 0:3)
    end

    @testset "seperable conv" begin
        u = [1, 2, 3, 2, 1]
        v = [6, 7, 3, 2]
        A = [1 2 3 4 5 6 7;
             8 9 10 11 12 13 14;
             15 16 17 18 19 20 21;
             22 23 24 25 26 27 28]
        exp = [6 19 35 53 71 89 107 77 33 14;
               60 148 217 285 339 393 447 315 134 56;
               204 478 658 822 930 1038 1146 798 338 140;
               468 1062 1400 1684 1828 1972 2116 1456 614 252;
               636 1426 1848 2188 2332 2476 2620 1792 754 308;
               624 1388 1778 2082 2190 2298 2406 1638 688 280;
               354 785 1001 1167 1221 1275 1329 903 379 154;
               132 292 371 431 449 467 485 329 138 56]
        @test_broken conv(u, v, A) == exp

        fu = convert(Array{Float64}, u)
        fv = convert(Array{Float64}, v)
        fA = convert(Array{Float64}, A)
        fexp = convert(Array{Float64}, exp)
        @test conv(fu, fv, fA) ≈ fexp

    end

    @testset "conv-ND" begin
        # is it safe to assume that if conv works for
        # int/float/complex in 1 and 2 D, it does in ND?
        a = convert(Array, reshape(1:27, (3, 3, 3)))
        b = ones(Int64, 2, 2, 2)
        exp = convert(Array, reshape([1, 3, 5, 3,
                                      5, 12, 16, 9,
                                      11, 24, 28, 15,
                                      7, 15, 17, 9,

                                      11, 24, 28, 15,
                                      28, 60, 68, 36,
                                      40, 84, 92, 48,
                                      23, 48, 52, 27,
                                      29, 60, 64, 33,

                                      64, 132, 140, 72,
                                      76, 156, 164, 84,
                                      41, 84, 88, 45,
                                      19, 39, 41, 21,

                                      41, 84, 88, 45,
                                      47, 96,  100, 51,
                                      25, 51, 53, 27], (4, 4, 4)))
        @test conv(a, b) == exp


        #6D, trivial, just to see if it works
        a = ones(2, 2, 2, 2, 2, 2)
        b = ones(1, 1, 1, 1, 1, 1)
        @test conv(a, b) == a

        # promote dims to largest
        a = cat([fill(n, 3, 3) for n in 1:6]..., dims=3)
        b = ones(Int64, 2, 2)
        expf1 = conv(a[:, :, 1], b)
        exp = cat([expf1 * n for n in 1:6]..., dims=3)
        @test conv(a, b) == exp
        fb = convert(Array{Float64}, b)
        @test conv(a, fb) ≈ convert(Array{Float64}, exp)
    end

    @testset "Overlap-Save" begin
        to_sizetuple(diml, N) = ntuple(_ -> diml, N)
        function os_test_data(::Type{T}, diml, N) where T
            s = to_sizetuple(diml, N)
            arr = rand(T, s)
            s, arr
        end
        function test_os(::Type{T}, nu, nv, N, nfft) where T
            nffts = to_sizetuple(nfft, N)
            su, u = os_test_data(T, nu, N)
            sv, v = os_test_data(T, nv, N)
            sout = su .+ sv .- 1
            out = _conv_similar(u, sout, axes(u), axes(v))
            unsafe_conv_kern_os!(out, u, v, su, sv, sout, nffts)
            os_out = copy(out)
            fft_nfft = nextfastfft(sout)
            _conv_kern_fft!(out, u, v, su, sv, sout, fft_nfft)
            @test out ≈ os_out
        end
        Ns = [1, 2, 3]
        eltypes = [Float32, Float64, Complex{Float64}]
        nlarge = 128

        regular_nsmall = [12, 128]
        for numdim in Ns, elt in eltypes, nsmall in regular_nsmall
            nfft = optimalfftfiltlength(nsmall, nlarge)
            test_os(elt, nlarge, nsmall, Val{numdim}(), nfft)
        end

        # small = 12, fft = 256 exercises output being smaller than the normal
        # valid region of OS convolution block
        #
        # small = 13, fft = 32 makes sout evenly divisble by nfft - sv + 1
        # small = 12, fft = 32 does not
        adversarial_nsmall_nfft = [(12, 256), (13, 32), (12, 32)]
        for (nsmall, nfft) in adversarial_nsmall_nfft
            test_os(Float64, nlarge, nsmall, Val{1}(), nfft)
        end

        # Should not bug:
        conv(zeros(4, 7, 1), zeros(3, 3, 3))

        # three blocks need to be padded in the following case:
        test_os(Float64, 25, 4, Val{1}(), 16)
    end
end

@testset "xcorr" begin
    a = [1, 2, 3]
    b = [4, 5]
    exp = [5, 14, 23, 12]
    @test xcorr([1, 2], [3, 4]) == [4, 11, 6]
    @test xcorr(a, b) == exp
    @test xcorr(a, b, padmode = :longest) == [0, 5, 14, 23, 12]
    @test xcorr(a, b, padmode = :none) == exp
    @test xcorr([1, 2], [3, 4, 5]) == [5, 14, 11, 6]
    @test xcorr([1, 2], [3, 4, 5], padmode = :longest) == [5, 14, 11, 6, 0]
    @test xcorr([1, 2], [3, 4, 5], padmode = :none) == [5, 14, 11, 6]
    @test xcorr([1.0im], [1.0im]) == [1]
    @test xcorr([1, 2, 3]*1.0im, ComplexF64[4, 5]) ≈ [5, 14, 23, 12]*im
    @test xcorr([1, 2]*1.0im, ComplexF64[3, 4, 5]) ≈ [5, 14, 11, 6]*im
    @test xcorr(ComplexF64[1, 2, 3], [4, 5]*1.0im) ≈ -[5, 14, 23, 12]*im
    @test xcorr(ComplexF64[1, 2], [3, 4, 5]*1.0im) ≈ -[5, 14, 11, 6]*im
    @test xcorr([1, 2, 3]*1.0im, [4, 5]*1.0im) ≈ [5, 14, 23, 12]
    @test xcorr([1, 2]*1.0im, [3, 4, 5]*1.0im) ≈ [5, 14, 11, 6]

    # Base Julia issue #17351
    let
        x = rand(10)
        u = rand(10, 3)
        su = view(u, :, 1)
        @test size(@inferred(xcorr(su, x))) == (19,)
        @test size(@inferred(xcorr(x, su))) == (19,)
    end

    # xcorr only supports 1d inputs at the moment
    @test_throws MethodError xcorr(rand(2, 2), rand(2, 2))

    off_a = OffsetVector(a, -1:1)
    off_b = OffsetVector(b, 0:1)
    @test xcorr(off_a, off_b, padmode = :none) == OffsetVector(exp, -2:1)

    # Issue #288
    @test xcorr(off_a, off_b, padmode = :longest) == OffsetVector(vcat(0, exp), -3:1)
end

@testset "deconv" begin
    # Test for issue #188: deconv mutates inputs
    let b = [4.0, 2.0, 1.0]; a = [2.0, 1.0]
        bb = b[:]
        aa = a[:]
        c = deconv(b,a)
        @test c == [2.0, 0.0]
        @test a == aa
        @test b == bb
    end
end
