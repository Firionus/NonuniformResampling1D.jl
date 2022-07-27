import Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

using Documenter, DocumenterMarkdown, NonuniformResampling1D

DocMeta.setdocmeta!(NonuniformResampling1D, :DocTestSetup, 
    :(using NonuniformResampling1D); recursive=true
)

makedocs(
    format=Markdown(),
    modules=[NonuniformResampling1D],
)

cp(joinpath(@__DIR__, "build", "README.md"), joinpath(@__DIR__, "..", "README.md"), force=true)