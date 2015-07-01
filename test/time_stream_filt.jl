using DSP

f = FIRFilter(3.1415926535897)
y = Array(Float64, 400)
x = rand(Float64, 100)

filt!(y, f, x)

function run_test()
    x  = rand(Float64, 10_000_000)
    y  = Array(Float64, 40_000_000)
    ff = FIRFilter(3.1415926535897)
    t  = @timed filt!(y, ff, x)
    return t
end

t = run_test()
print(t)

# Master
# (31415927,1.598210385,16,0.0,Base.GC_Diff(16,0,0,0,1,0,0,0,0,0,0))

# new
# # (31415927,3.510104474,16,0.0,Base.GC_Diff(16,0,0,0,1,0,0,0,0,0,0))
