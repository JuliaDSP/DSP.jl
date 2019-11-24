using Documenter, DSP

DocMeta.setdocmeta!(DSP, :DocTestSetup, :(using DSP); recursive=true)

makedocs(
    modules = [DSP],
    format = Documenter.HTML(),
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
