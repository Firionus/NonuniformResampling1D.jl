using Pkg

Pkg.activate(".")
Pkg.instantiate()

push!(LOAD_PATH,"../src/")

using Documenter, DocumenterMarkdown, NonuniformResampling1D

DocMeta.setdocmeta!(NonuniformResampling1D, :DocTestSetup, 
    :(using NonuniformResampling1D); recursive=true
)

makedocs(
    format=Markdown(),
    modules=[NonuniformResampling1D]
)

cp("build/README.md", "../README.md", force=true)