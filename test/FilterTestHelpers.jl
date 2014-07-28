module FilterTestHelpers
using DSP, Base.Test

export tffilter_eq, zpkfilter_eq, tffilter_accuracy, zpkfilter_accuracy

function lt(a, b)
    if abs(real(a) - real(b)) > 1e-10
        isless(real(a), real(b))
    else
        isless(imag(a), imag(b))
    end
end

function tffilter_eq(f1, f2)
    b1, a1 = (coefb(f1), coefa(f1))
    b2, a2 = (coefb(f2), coefa(f2))
    @test_approx_eq float64(b1) float64(b2)
    @test_approx_eq float64(a1) float64(a2)
end

function zpkfilter_eq(f1, f2)
    if !isempty(f1.z) || !isempty(f2.z)
        @test_approx_eq complex128(sort(f1.z, lt=lt)) complex128(sort(f2.z, lt=lt))
    end
    @test_approx_eq complex128(sort(f1.p, lt=lt)) complex128(sort(f2.p, lt=lt))
    @test_approx_eq float64(f1.k) float64(f2.k)
end

function zpkfilter_eq(f1, f2, eps)
    if !isempty(f1.z) || !isempty(f2.z)
        @test_approx_eq_eps complex128(sort(f1.z, lt=lt)) complex128(sort(f2.z, lt=lt)) eps
    end
    @test_approx_eq_eps complex128(sort(f1.p, lt=lt)) complex128(sort(f2.p, lt=lt)) eps
    @test_approx_eq_eps float64(f1.k) float64(f2.k) eps
end

loss(x::Real, y::Real) = abs(float(x) - float(y))/eps(float(x))
loss(x::Union(Real,Complex), y::Union(Real,Complex)) = loss(real(x), real(y)) + loss(imag(x), imag(y))
loss(x::AbstractVector, y::AbstractVector) = sum(map(loss, x, y))

function accuracy_check(err1, err2, part, relerr=1)
    try
        @test err1 <= relerr*err2
    catch e
        println("Filter 1 $part error (ULP): ", float64(err1))
        println("Filter 2 $part error (ULP): ", float64(err2))
        println("Ratio: ", float64(err1/err2))
        rethrow(e)
    end
end

function tffilter_accuracy(f1, f2, accurate_f)
    b1, a1 = (coefb(f1), coefa(f1))
    b2, a2 = (coefb(f2), coefa(f2))
    accurate_b, accurate_a = (coefb(accurate_f), coefa(accurate_f))
    accuracy_check(loss(b1, accurate_b), loss(b2, accurate_b), "b")
    accuracy_check(loss(a1, accurate_a), loss(a2, accurate_a), "a")
end

function zpkfilter_accuracy(f1, f2, accurate_f; relerr=1, eps=nothing)
    z1, p1 = sort(f1.z, lt=lt), sort(f1.p, lt=lt)
    z2, p2 = sort(f2.z, lt=lt), sort(f2.p, lt=lt)
    accurate_z, accurate_p = sort(accurate_f.z, lt=lt), sort(accurate_f.p, lt=lt)
    if !isempty(z1) || !isempty(z2) || !isempty(accurate_z)
        if eps != nothing
            @test_approx_eq_eps z1 accurate_z eps
            @test_approx_eq_eps z2 accurate_z eps
        else
            @test_approx_eq z1 accurate_z
            @test_approx_eq z2 accurate_z
        end
        accuracy_check(loss(z1, accurate_z), loss(z2, accurate_z), "z", relerr)
    end
    if eps != nothing
        @test_approx_eq_eps p1 accurate_p eps
        @test_approx_eq_eps p2 accurate_p eps
        @test_approx_eq_eps f1.k accurate_f.k eps
        @test_approx_eq_eps f2.k accurate_f.k eps
    else
        @test_approx_eq p1 accurate_p
        @test_approx_eq p2 accurate_p
        @test_approx_eq f1.k accurate_f.k
        @test_approx_eq f2.k accurate_f.k
    end
    accuracy_check(loss(p1, accurate_p), loss(p2, accurate_p), "p", relerr)
    accuracy_check(loss(f1.k, accurate_f.k), loss(f2.k, accurate_f.k), "k", relerr)
end
end