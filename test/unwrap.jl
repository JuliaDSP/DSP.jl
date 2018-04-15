using DSP, Compat, Compat.Test

@testset "Unwrap 1D" begin
    @test unwrap([0.1, 0.2, 0.3, 0.4]) ≈ [0.1, 0.2, 0.3, 0.4]
    @test unwrap([0.1, 0.2 + 2pi, 0.3, 0.4]) ≈ [0.1, 0.2, 0.3, 0.4]
    @test unwrap([0.1, 0.2 - 2pi, 0.3, 0.4]) ≈ [0.1, 0.2, 0.3, 0.4]
    @test unwrap([0.1, 0.2 - 2pi, 0.3 - 2pi, 0.4]) ≈ [0.1, 0.2, 0.3, 0.4]
    @test unwrap([0.1 + 2pi, 0.2, 0.3, 0.4]) ≈ [0.1 + 2pi, 0.2 + 2pi, 0.3 + 2pi, 0.4 + 2pi]

    test_v = [0.1, 0.2, 0.3 + 2pi, 0.4]
    res_v = unwrap(test_v)
    @test test_v ≈ [0.1, 0.2, 0.3 + 2pi, 0.4]
    unwrap!(test_v)
    @test test_v ≈ [0.1, 0.2, 0.3, 0.4]

    # test unwrapping with other ranges
    unwrapped = [1.0:100;]
    wrapped = Float64[i % 10 for i in unwrapped]
    @test unwrap(wrapped, range=10) ≈ unwrapped

    types = (Float32, Float64, BigFloat)
    for T in types
        srand(1234)
        A_unwrapped = collect(linspace(0, 4convert(T, π), 10))
        A_wrapped = A_unwrapped .% (2convert(T, π))

        @test unwrap(A_wrapped) ≈ A_unwrapped
        unwrap!(A_wrapped)
        @test A_wrapped ≈ A_unwrapped

        A_unwrapped_range = collect(linspace(0, 4, 10))
        test_range = convert(T, 2)
        A_wrapped_range = A_unwrapped_range .% test_range
        @test unwrap(A_wrapped_range; range=test_range) ≈ A_unwrapped_range
    end
end

@testset "Unwrap 2D" begin
    types = (Float32, Float64, BigFloat)
    for T in types
        srand(1234)
        v_unwrapped = collect(linspace(0, 4convert(T, π), 7))
        A_unwrapped = v_unwrapped .+ v_unwrapped'
        A_wrapped = A_unwrapped .% (2convert(T, π))

        test_unwrapped = unwrap(A_wrapped)
        d = first(A_unwrapped) - first(test_unwrapped)
        @test (test_unwrapped + d) ≈ A_unwrapped
        unwrap!(A_wrapped)
        d = first(A_unwrapped) - first(A_wrapped)
        @test (A_wrapped + d) ≈ A_unwrapped

        v_unwrapped_range = collect(linspace(0, 4, 7))
        A_unwrapped_range = v_unwrapped_range .+ v_unwrapped_range'
        test_range = convert(T, 2)
        A_wrapped_range = A_unwrapped_range .% test_range

        test_unwrapped_range = unwrap(A_wrapped_range; range=test_range)
        d = first(A_unwrapped_range) - first(test_unwrapped_range)
        @test (test_unwrapped_range + d) ≈ A_unwrapped_range

        # Test wrap_around
        # after unwrapping, pixels at borders should be equal to corresponding pixels
        # on other side
        wrap_around = (true, true)
        wa_vec = linspace(0, 4convert(T, π), 10)
        wa_uw = wa_vec .+ zeros(10)'
        # make periodic
        wa_uw[end, :] = wa_uw[1, :]
        wa_w = wa_uw .% (2π)
        wa_test = unwrap(wa_w, wrap_around=wrap_around, seed=0)
        # with wrap-around, the borders should be equal, but for this problem the
        # image may not be recovered exactly
        @test wa_test[:, 1] ≈ wa_test[:, end]
        @test wa_test[end, :] ≈ wa_test[1, :]
        # In this case, calling unwrap w/o wrap_around does not recover the borders
        wa_test_nowa = unwrap(wa_w)
        @test !(wa_test_nowa[end, :] ≈ wa_test_nowa[1, :])

    end
end

@testset "Unwrap 3D" begin
    types = (Float32, Float64, BigFloat)
    f(x, y, z) = 0.1x^2 - 2y + 2z
    f_wraparound2(x, y, z) = 5*sin(x) + 2*cos(y) + z
    f_wraparound3(x, y, z) = 5*sin(x) + 2*cos(y) - 4*cos(z)
    for T in types
        grid = linspace(zero(T), convert(T, 2π), 11)
        f_uw = f.(grid, grid', reshape(grid, 1, 1, :))
        f_wr = f_uw .% (2π)
        uw_test = unwrap(f_wr)
        offset = first(f_uw) - first(uw_test)
        @test isapprox(uw_test + offset, f_uw, atol=1e-8)
        # test in-place version
        unwrap!(f_wr)
        offset = first(f_uw) - first(f_wr)
        @test isapprox(f_wr + offset, f_uw, atol=1e-8)

        f_uw = f_wraparound2.(grid, grid', reshape(grid, 1, 1, :))
        f_wr = f_uw .% (2π)
        uw_test = unwrap(f_wr, wrap_around=(true, true, false))
        offset = first(f_uw) - first(uw_test)
        @test isapprox(uw_test + offset, f_uw, atol=1e-8)
        # test in-place version
        unwrap!(f_wr, wrap_around=(true, true, false))
        offset = first(f_uw) - first(f_wr)
        @test isapprox(f_wr + offset, f_uw, atol=1e-8)

        f_uw = f_wraparound3.(grid, grid', reshape(grid, 1, 1, :))
        f_wr = f_uw .% (2π)
        uw_test = unwrap(f_wr, wrap_around=(true, true, true))
        offset = first(f_uw) - first(uw_test)
        @test isapprox(uw_test + offset, f_uw, atol=1e-8)
        # test in-place version
        unwrap!(f_wr, wrap_around=(true, true, true))
        offset = first(f_uw) - first(f_wr)
        @test isapprox(f_wr + offset, f_uw, atol=1e-8)
    end
end
