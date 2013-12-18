using DSP, Base.Test

@test_approx_eq(unwrap([0.1, 0.2, 0.3, 0.4]), [0.1, 0.2, 0.3, 0.4])
@test_approx_eq(unwrap([0.1, 0.2 + 2pi, 0.3, 0.4]), [0.1, 0.2, 0.3, 0.4])
@test_approx_eq(unwrap([0.1, 0.2 - 2pi, 0.3, 0.4]), [0.1, 0.2, 0.3, 0.4])
@test_approx_eq(unwrap([0.1, 0.2 - 2pi, 0.3 - 2pi, 0.4]), [0.1, 0.2, 0.3, 0.4])
@test_approx_eq(unwrap([0.1 + 2pi, 0.2, 0.3, 0.4]),
                [0.1 + 2pi, 0.2 + 2pi, 0.3 + 2pi, 0.4 + 2pi])
@test_approx_eq(unwrap([0.1, 0.2 + 6pi, 0.3, 0.4]), [0.1, 0.2, 0.3, 0.4])

test_v = [0.1, 0.2, 0.3 + 2pi, 0.4]
res_v = unwrap(test_v)
@test_approx_eq(test_v, [0.1, 0.2, 0.3 + 2pi, 0.4])
unwrap!(test_v)
@test_approx_eq(test_v, [0.1, 0.2, 0.3, 0.4])
