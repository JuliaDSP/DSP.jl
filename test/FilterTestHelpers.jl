module FilterTestHelpers
using DSP, Test

export tffilter_eq, zpkfilter_eq, tffilter_accuracy, zpkfilter_accuracy,
       matrix_to_sosfilter, sosfilter_to_matrix

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
    @test map(Float64, b1) ≈ map(Float64, b2)
    @test map(Float64, a1) ≈ map(Float64, a2)
end

function zpkfilter_eq(f1, f2)
    if !isempty(f1.z) || !isempty(f2.z)
        @test map(ComplexF64, sort(f1.z, lt=lt)) ≈ map(ComplexF64, sort(f2.z, lt=lt))
    end
    @test map(ComplexF64, sort(f1.p, lt=lt)) ≈ map(ComplexF64, sort(f2.p, lt=lt))
    @test map(Float64, f1.k) ≈ map(Float64, f2.k)
end

function zpkfilter_eq(f1, f2, eps)
    if !isempty(f1.z) || !isempty(f2.z)
        @test ≈(map(ComplexF64, sort(f1.z, lt=lt)), map(ComplexF64, sort(f2.z, lt=lt)), atol=eps)
    end
    @test ≈(map(ComplexF64, sort(f1.p, lt=lt)), map(ComplexF64, sort(f2.p, lt=lt)), atol=eps)
    @test ≈(map(Float64, f1.k), map(Float64, f2.k), atol=eps)
end

loss(x::Real, y::Real) = abs(float(x) - float(y))/eps(float(x))
loss(x::Union{Real,Complex}, y::Union{Real,Complex}) = loss(real(x), real(y)) + loss(imag(x), imag(y))
loss(x::AbstractVector, y::AbstractVector) = sum(map(loss, x, y))

function accuracy_check(err1, err2, part, relerr=1)
    try
        @test err1 <= relerr*err2
    catch e
        println("Filter 1 $part error (ULP): ", map(Float64, err1))
        println("Filter 2 $part error (ULP): ", map(Float64, err2))
        println("Ratio: ", map(Float64, err1/err2))
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

function zpkfilter_accuracy(f1, f2, accurate_f; relerr=1, compare_gain_at=nothing, eps=nothing)
    z1, p1 = sort(f1.z, lt=lt), sort(f1.p, lt=lt)
    z2, p2 = sort(f2.z, lt=lt), sort(f2.p, lt=lt)
    accurate_z, accurate_p = sort(accurate_f.z, lt=lt), sort(accurate_f.p, lt=lt)
    if !isempty(z1) || !isempty(z2) || !isempty(accurate_z)
        if eps !== nothing
            @test ≈(z1, accurate_z, atol=eps)
            @test ≈(z2, accurate_z, atol=eps)
        else
            @test z1 ≈ accurate_z
            @test z2 ≈ accurate_z
        end
        accuracy_check(loss(z1, accurate_z), loss(z2, accurate_z), "z", relerr)
    end
    if eps !== nothing
        @test ≈(p1, accurate_p, atol=eps)
        @test ≈(p2, accurate_p, atol=eps)
        @test ≈(f1.k, accurate_f.k, atol=eps)
        @test ≈(f2.k, accurate_f.k, atol=eps)
    else
        @test p1 ≈ accurate_p
        @test p2 ≈ accurate_p
        @test f1.k ≈ accurate_f.k
        @test f2.k ≈ accurate_f.k
    end
    accuracy_check(loss(p1, accurate_p), loss(p2, accurate_p), "p", relerr)
    if compare_gain_at === nothing
        accuracy_check(loss(f1.k, accurate_f.k), loss(f2.k, accurate_f.k), "k", relerr)
    else
        jω = compare_gain_at*im
        gain1 = Float64(abs(f1.k * prod(Complex{BigFloat}.(z1) .- jω) / prod(Complex{BigFloat}.(p1) .- jω)))
        gain2 = Float64(abs(f2.k * prod(Complex{BigFloat}.(z2) .- jω) / prod(Complex{BigFloat}.(p2) .- jω)))
        accurate_gain = abs(accurate_f.k * prod(Complex{BigFloat}.(accurate_z) .- jω) / prod(Complex{BigFloat}.(accurate_p) .- jω))
        accuracy_check(loss(gain1, accurate_gain), loss(gain2, accurate_gain), "k", relerr)
    end
end

# Convert an SOS matrix, as returned by MATLAB, to an SecondOrderSections
matrix_to_sosfilter(sos::Matrix, g::Number) =
    SecondOrderSections([Biquad(sos[i, 1], sos[i, 2], sos[i, 3],
                            sos[i, 4], sos[i, 5], sos[i, 6])
               for i = 1:size(sos, 1)], g)

# Convert an SecondOrderSections to an SOS matrix as returned by MATLAB
function sosfilter_to_matrix(sos::SecondOrderSections)
    A = ones(length(sos.biquads), 6)
    for (i, biquad) in enumerate(sos.biquads)
        A[i, 1] = biquad.b0
        A[i, 2] = biquad.b1
        A[i, 3] = biquad.b2
        A[i, 5] = biquad.a1
        A[i, 6] = biquad.a2
    end
    A
end

# Show the poles and zeros in each biquad
# This is not currently used for testing, but is useful for debugging
function sosfilter_poles_zeros(sos::SecondOrderSections)
    z = fill(complex(nan(T)), 2, length(sos.biquads))
    p = fill(complex(nan(T)), 2, length(sos.biquads))
    for (i, biquad) in enumerate(sos.biquads)
        zpk = convert(ZeroPoleGain, biquad)
        z[1:length(zpk.z), i] = zpk.z
        p[1:length(zpk.p), i] = zpk.p
    end
    (z, p)
end
end
