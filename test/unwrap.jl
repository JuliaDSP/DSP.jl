using DSP, Test
using Random: MersenneTwister, default_rng, seed!
using StableRNGs
using Statistics: mean

@testset "Unwrap 1D" begin
    unwrapped = [0.1, 0.2, 0.3, 0.4]
    @test unwrap([0.1, 0.2, 0.3, 0.4]) ≈ unwrapped
    @test unwrap([0.1, 0.2 + 2pi, 0.3, 0.4]) ≈ unwrapped
    @test unwrap([0.1, 0.2 - 2pi, 0.3, 0.4]) ≈ unwrapped
    @test unwrap([0.1, 0.2 - 2pi, 0.3 - 2pi, 0.4]) ≈ unwrapped
    @test unwrap([0.1 + 2pi, 0.2, 0.3, 0.4]) ≈ unwrapped .+ 2pi
    @test unwrap([0.1, 0.2 + 6pi, 0.3, 0.4]) ≈ unwrapped

    test_v = [0.1, 0.2, 0.3 + 2pi, 0.4]
    res_v = unwrap(test_v)
    @test res_v ≈ unwrapped
    @test test_v == [0.1, 0.2, 0.3 + 2pi, 0.4]

    res_v .= 0
    unwrap!(res_v, test_v)
    @test res_v ≈ unwrapped
    @test test_v == [0.1, 0.2, 0.3 + 2pi, 0.4]

    @test unwrap!(test_v) === test_v
    @test test_v ≈ unwrapped

    # test unwrapping within multi-dimensional array
    wrapped = repeat([0.1, 0.2 + 2pi, 0.3, 0.4], 1, 2)
    unwrapped = hcat(unwrapped, unwrapped)
    @test unwrap(wrapped, dims=2) ≈ wrapped
    @test unwrap(wrapped, dims=1) ≈ unwrapped
    @test unwrap!(copy(wrapped), dims=2) ≈ wrapped
    @test unwrap!(copy(wrapped), dims=1) ≈ unwrapped

    # this should eventually default to the multi-dimensional case
    @test_throws ArgumentError unwrap!(similar(wrapped), wrapped)

    # test unwrapping with other ranges
    unwrapped = [1.0:100;]
    wrapped = Float64[i % 10 for i in unwrapped]
    @test unwrap(wrapped, range=10) ≈ unwrapped

    # test generically typed unwrapping
    types = (Float32, Float64, BigFloat)
    for T in types
        A_unwrapped = collect(range(0, stop=4convert(T, π), length=10))
        A_wrapped = A_unwrapped .% (2convert(T, π))

        @test unwrap(A_wrapped) ≈ A_unwrapped
        unwrap!(A_wrapped)
        @test A_wrapped ≈ A_unwrapped

        A_unwrapped_range = collect(range(0, stop=4, length=10))
        test_range = convert(T, 2)
        A_wrapped_range = A_unwrapped_range .% test_range
        @test unwrap(A_wrapped_range; range=test_range) ≈ A_unwrapped_range
    end
end

# tests for multi-dimensional unwrapping
@testset "Unwrap 2D" begin
    types = (Float32, Float64, BigFloat)
    for T in types
        v_unwrapped = collect(range(0, stop=4convert(T, π), length=7))
        A_unwrapped = v_unwrapped .+ v_unwrapped'
        A_wrapped = A_unwrapped .% (2convert(T, π))

        test_unwrapped = unwrap(A_wrapped, dims=1:2)
        d = first(A_unwrapped) - first(test_unwrapped)
        @test (test_unwrapped .+ d) ≈ A_unwrapped
        unwrap!(A_wrapped, dims=1:2)
        d = first(A_unwrapped) - first(A_wrapped)
        @test (A_wrapped .+ d) ≈ A_unwrapped

        v_unwrapped_range = collect(range(0, stop=4, length=7))
        A_unwrapped_range = v_unwrapped_range .+ v_unwrapped_range'
        test_range = convert(T, 2)
        A_wrapped_range = A_unwrapped_range .% test_range

        test_unwrapped_range = unwrap(A_wrapped_range, dims=1:2; range=test_range)
        d = first(A_unwrapped_range) - first(test_unwrapped_range)
        @test (test_unwrapped_range .+ d) ≈ A_unwrapped_range

        # Test circular_dims
        # after unwrapping, pixels at borders should be equal to corresponding pixels
        # on other side
        circular_dims = (true, true)
        wa_vec = range(0, stop=4convert(T, π), length=10)
        wa_uw = wa_vec .+ zeros(10)'
        # make periodic
        wa_uw[end, :] = wa_uw[1, :]
        wa_w = wa_uw .% (2π)
        wa_test = unwrap(wa_w, dims=1:2, circular_dims=circular_dims, rng=MersenneTwister(0))
        # with wrap-around, the borders should be equal, but for this problem the
        # image may not be recovered exactly
        @test wa_test[:, 1] ≈ wa_test[:, end]
        @test wa_test[end, :] ≈ wa_test[1, :]
        # In this case, calling unwrap w/o circular_dims does not recover the borders
        wa_test_nowa = unwrap(wa_w, dims=1:2)
        @test !(wa_test_nowa[end, :] ≈ wa_test_nowa[1, :])

    end
end

@testset "Unwrap 3D" begin
    types = (Float32, Float64, BigFloat)

    f_nowrap(x, y, z) = 0.1x^2 - 2y + 2z
    f_wraparound2(x, y, z) = 5*sin(x) + 2*cos(y) + z
    f_wraparound3(x, y, z) = 5*sin(x) + 2*cos(y) - 4*cos(z)

    for T in types
        grid = range(zero(T), stop=2convert(T, π), length=11)

        for (func, circular_dims) in (
            (f_nowrap,      (false, false, false)),
            (f_wraparound2, (true, true, false)),
            (f_wraparound3, (true, true, true))
        )
            f_uw = func.(grid, grid', reshape(grid, 1, 1, :))
            f_wr = f_uw .% (2convert(T, π))
            uw_test = unwrap(f_wr; dims=1:3, circular_dims)
            offset = first(f_uw) - first(uw_test)
            @test (uw_test .+ offset) ≈ f_uw rtol = eps(T)  # oop

            # test in-place version
            unwrap!(f_wr; dims=1:3, circular_dims)
            offset = first(f_uw) - first(f_wr)
            @test (f_wr .+ offset) ≈ f_uw rtol = eps(T)     # ip
        end
    end
end

@testset "reproducible unwrap" begin
    function measure(s, rng, seed=rand(UInt))
        u1 = unwrap(s; dims=1:ndims(s), rng=seed!(rng, seed))
        u2 = unwrap(s; dims=1:ndims(s), rng=seed!(rng, seed))
        return sqrt(mean(abs2, u1 - u2))
    end

    x1 = 1234 * rand(10_000)
    x2 = 1357 * rand(200, 200)
    x3 = 1248 * rand(30, 30, 30)

    for x in (x1, x2, x3)
        @test measure(x, default_rng()) == measure(x, MersenneTwister()) ==
              measure(x, StableRNG(10)) == 0.0
    end
end
