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

volume_pow2s = 2:28

function prepare_test_objects(N, edge_pow_u, edge_pow_v)
    su = ntuple(i -> edge_pow_u ^ 2, N)
    sv = ntuple(i -> edge_pow_v ^ 2, N)
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
