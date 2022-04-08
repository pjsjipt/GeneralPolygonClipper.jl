using GeneralPolygonClipper
using Documenter

DocMeta.setdocmeta!(GeneralPolygonClipper, :DocTestSetup, :(using GeneralPolygonClipper); recursive=true)

makedocs(;
    modules=[GeneralPolygonClipper],
    authors="Paulo Jabardo <pjabardo@ipt.br>",
    repo="https://github.com/pjsjipt/GeneralPolygonClipper.jl/blob/{commit}{path}#{line}",
    sitename="GeneralPolygonClipper.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
