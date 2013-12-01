using DSP, Base.Test

# Test dpss against dpss computed with MATLAB
d1 = dpss(128, 4)
d2 = readdlm(joinpath(dirname(@__FILE__), "data", "dpss128,4.txt"), '\t')
@test_approx_eq d1 d2
