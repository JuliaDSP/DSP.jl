using DSP, BenchmarkTools, PyPlot, Statistics
using DSP: _conv_kern_fft!, _conv_kern_os!, _conv_similar, optimalfftfiltlength, nextfastfft

struct BenchTestResult
    usize::Vector{Int}
    vsize::Vector{Int}
    fft_time::Float64
    fft_memory::Float64
    os_time::Float64
    os_memory::Float64
end

struct BenchTestSet
    ndim::Integer
    results::Vector{BenchTestResult}
end

volume_pow2s = 2:20

function prepare_test_objects(N, edge_pow_u, edge_pow_v)
    su = ntuple(i -> 2 ^ edge_pow_u, N)
    sv = ntuple(i -> 2 ^ edge_pow_v, N)
    sout = su .+ sv .- 1
    u = rand(Float64, su)
    v = rand(Float64, sv)
    out = _conv_similar(u, v, sout)
    os_nffts = map((nu, nv)-> optimalfftfiltlength(nu, nv), su, sv)
    out, u, v, su, sv, sout, os_nffts
end

dimresults = Vector{BenchTestSet}(undef, 4)
for N = 1:4
    valid_edge_pows = unique(max.(2, div.(volume_pow2s, N)))
    results = BenchTestResult[]
    for u_edge_pow2 in valid_edge_pows
        last_idx = searchsortedlast(valid_edge_pows, u_edge_pow2)
        for v_edge_pow2 in valid_edge_pows[1:last_idx]
            @show N, u_edge_pow2, v_edge_pow2
            out, u, v, su, sv, sout, os_nffts = prepare_test_objects(N, u_edge_pow2, v_edge_pow2)
            nffts = nextfastfft(sout)
            t_fft = @benchmark(_conv_kern_fft!($out, $u, $v, $su, $sv, $sout, $nffts))
            t_os = @benchmark(_conv_kern_os!($out, $u, $v, $su, $sv, $sout, $os_nffts))
            push!(
                results,
                BenchTestResult(
                    collect(su), collect(sv), median(t_fft.times),
                    median(t_os.times), t_fft.memory, t_os.memory
                )
            )
        end
    end
    dimresults[N] = BenchTestSet(N, results)
end

per_arr = Matrix{Float64}(undef, 19, 19)

fig, axs = subplots(nrows = 2, ncols = 2, figsize = (6.4, 5.5))

I = nothing
for N = 1:4

    fill!(per_arr, NaN)
    pos = N == 1 ? (1, 1) :
        N == 2 ? (2, 1) :
        N == 3 ? (1, 2) :
        N == 4 ? (2, 2) :
        error()

    idx = CartesianIndex(pos)

    for res in dimresults[N].results
        upos = round(Int, log2(prod(res.usize))) - 1
        vpos = round(Int, log2(prod(res.vsize))) - 1
        per_arr[vpos, upos] = 100 * (res.os_memory / res.fft_memory)
    end
    global I = axs[idx].imshow(per_arr, cmap = "RdBu_r", vmin = 0, vmax = 200, origin = "lower", extent = [1.5, 20.5, 1.5, 20.5])
    axs[idx].set_title("N = $N")
end
fig.tight_layout(rect = [0.02, 0.06, .98, .94])
fig.colorbar(I, ax = axs[:])
hth = fig.text(0.5, 0.05, "\$log_2(length(u))\$", ha = "center")
htv = fig.text(0.05, 0.5, "\$log_2(length(v))\$", va = "center", rotation = "vertical")
htt = fig.suptitle("Overlap-Save Memory Usage (% of FFT)")
fig.savefig("os_grid.png")
