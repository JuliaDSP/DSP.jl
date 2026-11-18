# julia --startup-file=no benchmark/remez_compatibility.jl [baseline git revision]
# Keep this differential check outside the package tests: it requires git history.
using Test, Random, Logging
include("remez.jl")

function outcome(f, taps, bands, kwargs)
    logger = Test.TestLogger()
    result = with_logger(logger) do
        try
            f(taps, bands; kwargs...)
        catch e
            (typeof(e), sprint(showerror, e))
        end
    end
    return result, [(entry.level, entry.message) for entry in logger.logs]
end

@testset "Remez differential compatibility" begin
    rng = MersenneTwister(9138)
    for trial in 1:600
        taps = rand(rng, 6:201)
        nbands = rand(rng, 1:4)
        edges = sort!(rand(rng, 2nbands) ./ 2)
        if isodd(trial)
            edges[1], edges[end] = 0., .5
        end
        bands = [(edges[2i-1], edges[2i]) => (rand(rng), 10.0^rand(rng, -1:1))
                 for i in 1:nbands]
        kwargs = (; neg=rand(rng, Bool), grid_density=rand(rng, (4, 16, 32)),
                   maxiter=rand(rng, (1, 2, 10, 25, 40)))
        @testset "case $trial: taps=$taps, $kwargs" begin
            @test isequal(outcome(BaselineRemez.remez, taps, bands, kwargs),
                          outcome(CurrentRemez.remez, taps, bands, kwargs))
        end
    end
end
