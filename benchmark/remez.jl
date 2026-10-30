# julia --startup-file=no benchmark/remez.jl [baseline git revision]
# No non-stdlib dependencies; both implementations run in the same process.
using Printf

module BaselineRemez
    revision = isempty(ARGS) ? "7c798756cca39251da14bb5d38146feb2cdd1717" : ARGS[1]
    source = read(`git -C $(@__DIR__) show $revision:src/Filters/remez_fir.jl`, String)
    include_string(@__MODULE__, source, "baseline_remez.jl")
end

module CurrentRemez
    include("../src/Filters/remez_fir.jl")
end

const REMEZ_CASES = [
    ("small LP", 15, [(0., .2) => 1., (.3, .5) => 0.], (;)),
    ("LP odd", 151, [(0., .475) => 1., (.5, 1.) => 0.], (; Hz=2.)),
    ("LP even weighted", 152, [(0., .475) => (1., 1.), (.5, 1.) => (0., 2.)], (; Hz=2.)),
    ("HP", 51, [(0., .75) => 0., (.8, 1.) => 1.], (; Hz=2.)),
    ("BP", 180, [(0., .375) => 0., (.4, .5) => 1., (.525, 1.) => 0.], (; Hz=2., maxiter=30)),
    ("Hilbert even", 20, [(.1, .95) => 1.], (; Hz=2., neg=true)),
    ("Hilbert odd", 21, [(.1, .95) => 1.], (; Hz=2., neg=true)),
    ("diff even", 200, [(.01, .99) => (f -> f / 2, f -> 1 / f)], (; Hz=2., neg=true)),
    ("diff odd", 201, [(.05, .95) => (f -> f / 2, f -> 1 / f)], (; Hz=2., neg=true)),
    ("partial band", 63, [(.05, .2) => 1., (.25, .45) => 0.], (;)),
    ("dense grid", 63, [(0., .2) => 1., (.23, .5) => 0.], (; grid_density=64)),
    ("long LP", 511, [(0., .2) => 1., (.205, .5) => 0.], (; maxiter=100)),
    ("inverse sinc", 201, [(0., 2880.) => (f -> inv(sinc(f / 4800)), 1.),
                          (10000., 153600.) => (0., 100.)], (; Hz=307200.)),
]

function elapsed_batch(f, args, kwargs, repetitions)
    start = time_ns()
    for _ in 1:repetitions
        f(args...; kwargs...)
    end
    return (time_ns() - start) / repetitions
end

function compare_case(name, args, kwargs; newkwargs=kwargs, samples=31, repetitions=5)
    old, new = BaselineRemez.remez, CurrentRemez.remez
    reference = old(args...; kwargs...)
    @assert isequal(reference, new(args...; newkwargs...)) name
    # Compile both the filters and the timing harness before collecting samples.
    elapsed_batch(old, args, kwargs, 2)
    elapsed_batch(new, args, newkwargs, 2)
    oldbytes = @allocated old(args...; kwargs...)
    newbytes = @allocated new(args...; newkwargs...)
    oldtimes, newtimes = Float64[], Float64[]
    GC.gc()
    for i in 1:samples
        # Alternate order to reduce systematic thermal / frequency bias.
        if isodd(i)
            push!(oldtimes, elapsed_batch(old, args, kwargs, repetitions))
            push!(newtimes, elapsed_batch(new, args, newkwargs, repetitions))
        else
            push!(newtimes, elapsed_batch(new, args, newkwargs, repetitions))
            push!(oldtimes, elapsed_batch(old, args, kwargs, repetitions))
        end
    end
    oldtime = sort!(oldtimes)[cld(samples, 2)] / 1000
    newtime = sort!(newtimes)[cld(samples, 2)] / 1000
    @printf("%-24s %10.2f %10.2f %7.3f %10d %10d\n",
            name, oldtime, newtime, newtime / oldtime, oldbytes, newbytes)
    return (; oldtime, newtime, oldbytes, newbytes)
end

function main()
    println("Julia ", VERSION, " / ", Sys.CPU_NAME, " / threads=", Threads.nthreads())
    println("Times: median microseconds per call; bytes: allocations per call, after warmup.")
    @printf("%-24s %10s %10s %7s %10s %10s\n",
            "case", "old μs", "new μs", "new/old", "old bytes", "new bytes")
    for (name, taps, bands, kwargs) in REMEZ_CASES
        compare_case(name, (taps, bands), kwargs)
    end
    compare_case("three-argument weighted", (152, [0., .475, .5, 1.], [1., 0.]),
                 (; weight=[1., 2.], Hz=2.))
    for n in (20, 21)
        compare_case("three-arg Hilbert $n", (n, [.1, .95], [1.]),
                     (; Hz=2., filter_type=BaselineRemez.filter_type_hilbert);
                     newkwargs=(; Hz=2., filter_type=CurrentRemez.filter_type_hilbert))
    end
    for (n, edges) in ((200, [.01, .99]), (201, [.05, .95]))
        compare_case("three-arg diff $n", (n, edges, [1.]),
                     (; Hz=2., filter_type=BaselineRemez.filter_type_differentiator);
                     newkwargs=(; Hz=2., filter_type=CurrentRemez.filter_type_differentiator))
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
