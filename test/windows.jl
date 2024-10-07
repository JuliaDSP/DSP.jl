using FFTW, DSP, Test, DelimitedFiles

# not exported, but used for some internal testing
using DSP.Windows: makewindow
@testset "makewindow" begin
    # make sure we're evaluating the given function at the correct points for
    # all combinations of arguments
    @test makewindow(identity, 6, 0, false) ≈ [-0.5, -0.3, -0.1, 0.1, 0.3, 0.5]
    win = makewindow(identity, 6, 0, true)
    # doesn't matter whether the middle value is evaluated at +/- 0.5
    @test win ≈ [0.0, 1/6, 1/3, 1/2, -1/3, -1/6] || win ≈ [0.0, 1/6, 1/3, -1/2, -1/3, -1/6]
    # we actually only end up with one real zero here, because one of the
    # "padding" values is used to split the -1/2 from +1/2. For windows that go
    # to zero padding will actually add `padding` zeros, because evaluting at
    # +/- 0.5 will be zero.
    @test makewindow(identity, 6, 2, false) ≈ [-0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.0, 0.0]
    @test makewindow(identity, 6, 2, true) ≈ [0.0, 1/6, 1/3, 1/2, 0.0, -1/2, -1/3, -1/6]
    @test makewindow(identity, 5, 0, false) ≈ [-0.5, -0.25, 0.0, 0.25, 0.5]
    @test makewindow(identity, 5, 0, true) ≈ [0.0, 0.2, 0.4, -0.4, -0.2]
    @test makewindow(identity, 5, 2, false) ≈ [-0.5, -0.25, 0.0, 0.25, 0.5, 0.0, 0.0]
    @test makewindow(identity, 5, 2, true) ≈ [0.0, 0.2, 0.4, 0.0, 0.0, -0.4, -0.2]

    @test makewindow(x->42.0, 1, 0, false) ≈ [42.0]
    @test makewindow(x->42.0, 1, 0, true) ≈ [42.0]
    @test makewindow(x->42.0, 1, 2, false) ≈ [42.0, 0.0, 0.0]
    @test makewindow(x->42.0, 1, 2, true) ≈ [42.0, 0.0, 0.0]
end

@testset "dspss" begin
    # Test dpss against dpss computed with MATLAB
    d1 = dpss(128, 4)
    d2 = readdlm(joinpath(dirname(@__FILE__), "data", "dpss128,4.txt"), '\t')
    @test d1 ≈ d2

    # Test dpsseig against dpss from MATLAB
    lambda = [0.9999999997159923,0.9999999731146645,0.9999988168667646,
              0.9999680890685374,0.9994167543397652,0.9925560207018469,
              0.9368556668429153]
    @test dpsseig(d1, 4) ≈ lambda
end

@testset "common windows" begin
    # Checking Hanning, Hamming, Triangular, Bartlett, Bartlett-Hann, Tukey,
    # and Blackman windows against values computed with MATLAB.
    # Lanczos and cosine are checked against values generated with DSP.jl v0.4.0
    # to test for regressions, as there's no reference MATLAB implementation
    # Gaussian is compared against DSP.jl from commit da1b195, when the
    # implementation was corrected (see GH issue #204)
    @test rect(128) == ones(128)

    hanning_jl = hanning(128)
    hann_jl = hann(128)
    hanning_ml = readdlm(joinpath(dirname(@__FILE__), "data", "hanning128.txt"), '\t')
    @test hanning_jl ≈ hanning_ml
    @test hann_jl ≈ hanning_ml

    hamming_jl = hamming(128)
    hamming_ml = readdlm(joinpath(dirname(@__FILE__), "data", "hamming128.txt"), '\t')
    @test hamming_jl ≈ hamming_ml

    triang_jl = triang(128)
    triang_ml = readdlm(joinpath(dirname(@__FILE__), "data", "triang128.txt"), '\t')
    @test triang_jl ≈ triang_ml

    # for odd `n` the `triang` window should be the middle n-2 samples of the
    # bartlett
    @test triang(5) ≈ bartlett(7)[2:6]

    bartlett_jl = bartlett(128)
    bartlett_ml = readdlm(joinpath(dirname(@__FILE__), "data", "bartlett128.txt"), '\t')
    @test bartlett_jl ≈ bartlett_ml

    barthann_jl = bartlett_hann(128)
    barthann_ml = readdlm(joinpath(dirname(@__FILE__), "data", "bartlett_hann128.txt"), '\t')
    @test bartlett_jl ≈ bartlett_ml

    blackman_jl = blackman(128)
    blackman_ml = readdlm(joinpath(dirname(@__FILE__), "data", "blackman128.txt"), '\t')
    @test blackman_jl ≈ blackman_ml
    @test minimum(blackman_jl) == 0.0

    kaiser_jl = kaiser(128, 0.4/π)
    kaiser_ml = readdlm(joinpath(dirname(@__FILE__), "data", "kaiser128,0.4.txt"), '\t')
    @test kaiser_jl ≈ kaiser_ml

    gaussian_jl = gaussian(128, 0.2)
    gaussian_ref = readdlm(joinpath(dirname(@__FILE__), "data", "gaussian128,0.2.txt"), '\t')
    @test gaussian_jl ≈ gaussian_ref

    tukey_jl = tukey(128, 0.4)
    tukey_ml = readdlm(joinpath(dirname(@__FILE__), "data", "tukey128,0.4.txt"), '\t')
    @test tukey_jl ≈ tukey_ml

    @test tukey(128, 0) == rect(128)

    lanczos_jl = lanczos(128)
    lanczos_ref = readdlm(joinpath(dirname(@__FILE__), "data", "lanczos128.txt"), '\t')
    @test lanczos_jl ≈ lanczos_ref

    cosine_jl = cosine(128)
    cosine_ref = readdlm(joinpath(dirname(@__FILE__), "data", "cosine128.txt"), '\t')
    @test cosine_jl ≈ cosine_ref
end

zeroarg_wins = [rect, hanning, hamming, cosine, lanczos,
                bartlett, bartlett_hann, blackman, triang]
onearg_wins = [gaussian, kaiser, tukey]
@testset "zero-phase windows" begin
    for winf in zeroarg_wins
        if winf == triang
            # triang needs to be special-cased here because it has different
            # definitions for odd and even `n` (and also the underlying
            # continuous function changes with `n`), so it doesn't match the
            # zerophase assumption below
            @test triang(6, zerophase=true) ≈ [1.0, 0.75, 0.5, 0.25, 0.5, 0.75]
            @test triang(7, zerophase=true) ≈ [1.0, 0.75, 0.5, 0.25, 0.25, 0.5, 0.75]
        else
            @test winf(8, zerophase=true) ≈ ifftshift(winf(9)[1:8])
            # the `zerophase=false` version doesn't hit the center point, but the
            # zerophase one does, so this ends up introducing a 1/2-sample shift.
            @test winf(9, zerophase=true) ≈ ifftshift(winf(19)[2:2:end])
        end
    end

    # test the window functions that need extra args
    for winf in onearg_wins
        @test winf(8, 0.5, zerophase=true) == ifftshift(winf(9, 0.5)[1:8])
        @test winf(9, 0.5, zerophase=true) == ifftshift(winf(19, 0.5)[2:2:end])
    end

    @test dpss(8, 2, 1, zerophase=true)[:] == ifftshift(dpss(9, 2, 1)[1:8])
    # odd-length zerophase dpss windows not currently supported
    @test_throws ArgumentError dpss(9, 2, 1, zerophase=true)
end

@testset "window return types" begin
    # make sure return types are correct
    n = 12
    ft = typeof(1.0)
    for W in(:rect, :hanning, :hamming, :cosine, :lanczos, :triang, :bartlett, :bartlett_hann, :blackman)
        @eval @test Array{$ft,1} == typeof($W($n)) && length($W($n)) == $n
    end

    @test Array{ft,1} == typeof(tukey(n, 0.4)) && length(tukey(n, 0.4)) == n
    @test Array{ft,1} == typeof(gaussian(n, 0.4)) && length(gaussian(n, 0.4)) == n
    @test Array{ft,1} == typeof(kaiser(n, 0.4)) && length(kaiser(n, 0.4)) == n
    @test Array{ft,2} == typeof(dpss(n, 1.5)) && size(dpss(n, 1.5),1) == n  # size(,2) depends on the parameters
end

@testset "tensor product windows" begin
    # test all combinations of arguments. Each arg and kwarg can be not present,
    # a single value, or a 2-tuple
    function push_args!((w1e, w2e, w3e), arg, ::Type{T}, f) where {T}
        if arg isa T
            push!(w1e.args, f(arg))
            push!(w2e.args, f(arg))
            push!(w3e.args, f(arg))
        elseif arg isa Tuple
            push!(w1e.args, f(arg[1]))
            push!(w2e.args, f(arg[2]))
            push!(w3e.args, f(arg))
        end
    end

    for winf in [zeroarg_wins; onearg_wins],
         arg in (nothing, 0.4, (0.4, 0.5))
        # skip invalid combinations
        winf in zeroarg_wins && arg !== nothing && continue
        winf in onearg_wins && arg === nothing && continue
        for padding in (nothing, 4, (4,5)),
          zerophase in (nothing, true, (true,false))
            w1_expr = :($winf(15))
            w2_expr = :($winf(20))
            w3_expr = :($winf((15,20)))
            w_all = (w1_expr, w2_expr, w3_expr)

            push_args!(w_all, arg, Real, identity)
            push_args!(w_all, padding, Integer, s -> Expr(:kw, :padding, s))
            push_args!(w_all, zerophase, Bool, s -> Expr(:kw, :zerophase, s))

            w1 = eval(w1_expr)
            w2 = eval(w2_expr)
            w3 = eval(w3_expr)
            @test w3 ≈ w1 * w2'
        end
    end
end
