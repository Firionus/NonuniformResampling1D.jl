using Pkg

# build example plot
cd("../examples")
Pkg.activate(".")
Pkg.instantiate()

include("../examples/explanation_diagram.jl")
include("../examples/plot_windows.jl")

# build docs
cd("../docs")
Pkg.activate(".")
Pkg.instantiate()

push!(LOAD_PATH,"../src/")

using Documenter, DocumenterMarkdown, NonuniformResampling1D

DocMeta.setdocmeta!(NonuniformResampling1D, :DocTestSetup, 
    :(using NonuniformResampling1D); recursive=true
)

makedocs(
    format=Markdown(),
)

cp("build/README.md", "../README.md", force=true)