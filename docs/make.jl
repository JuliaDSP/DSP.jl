using Documenter, DSP

makedocs(
    modules = [DSP],
    format = :html,
    sitename = "DSP.jl",
    pages = Any[
        "Contents" => "contents.md",
        "periodograms.md",
        "estimation.md",
        "windows.md",
        "filters.md",
        "util.md",
        "convolutions.md",
        "lpc.md",
        "index.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaDSP/DSP.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
)
