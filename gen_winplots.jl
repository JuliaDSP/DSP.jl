# this script is used to generate the ascii plots used in the window docstrings

using UnicodePlots
using DSP

open(joinpath(@__DIR__, "src", "winplots.jl"), "w") do io

    println(io, """
    #####################################################
    # THIS FILE IS AUTO-GENERATED
    # to re-generate use the `gen_winplots.jl` script
    #####################################################
    """)

    for (winfunc, args) in [
        (rect, ()),
        (hanning, ()),
        (hamming, ()),
        (tukey, (0.4)),
        (cosine, ()),
        (lanczos, ()),
        (triang, ()),
        (bartlett, ()),
        (gaussian, (0.2)),
        (bartlett_hann, ()),
        (blackman, ()),
        (kaiser, (3)),
        (dpss, (2, 1)),
        ]

        ymax = winfunc == dpss ? 0.2 : 1.0
        n = winfunc == triang ? 7 : 69

        fname = split(string(winfunc), ".")[end]
        println(io, "const $(fname)_winplot = padplot(\"\"\"")
        # convert Nx1 matrices to vectors with [:] - necessary for dpss
        print(io, lineplot(winfunc(n, args...)[:], ylim=[0,ymax], xlim=[1,n], width=70, canvas=BlockCanvas))
        println(io, "\"\"\")\n")
    end
end
