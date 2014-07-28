require(joinpath(dirname(@__FILE__), "FilterTestHelpers.jl"))
using DSP, Base.Test, FilterTestHelpers

# Test that a numerically challenging filter (high order, clustered
# roots) has acceptable errors in its coefficients after conversion to
# SOS
f = ZPKFilter(ones(100), 0.99*ones(100), 1)
g = convert(SOSFilter, f)
tffilter_eq(convert(TFFilter, f), convert(TFFilter, g))

for f in (digitalfilter(Lowpass(0.5), Butterworth(1)), digitalfilter(Lowpass(0.5), Butterworth(2)),
          digitalfilter(Bandpass(0.25, 0.75), Butterworth(1)))
    for ftype1 in (ZPKFilter, TFFilter, BiquadFilter, SOSFilter)
        f2 = convert(ftype1, f)
        for ftype2 in (ZPKFilter, TFFilter, BiquadFilter, SOSFilter)
            f3 = convert(ftype2, f)
            try
                zpkfilter_eq(convert(ZPKFilter, f), convert(ZPKFilter, f3), sqrt(eps()))
            catch e
                println("Conversion from $ftype1 to $ftype2 failed:")
                rethrow(e)
            end
        end
    end
end

for proto in (Butterworth(3), Chebyshev1(3, 1), Chebyshev2(3, 1))
    f = digitalfilter(Lowpass(0.5), proto)
    for ftype1 in (ZPKFilter, TFFilter, SOSFilter)
        f2 = convert(ftype1, f)
        for ftype2 in (ZPKFilter, TFFilter, SOSFilter)
            f3 = convert(ftype2, f2)
            zpkfilter_eq(convert(ZPKFilter, f), convert(ZPKFilter, f3), 2e-5)
        end
    end
end
