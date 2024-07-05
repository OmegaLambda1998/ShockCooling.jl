using Documenter
push!(LOAD_PATH, "../src/")
using ShockCooling

DocMeta.setdocmeta!(ShockCooling, :DocTestSetup, :(using ShockCooling); recursive = true)

makedocs(
    sitename = "ShockCooling Documentation",
    modules = [ShockCooling],
    pages = ["ShockCooling" => "index.md", "API" => "api.md"],
    format = Documenter.HTML(assets = ["assets/favicon.ico"]),
)

deploydocs(repo = "github.com/OmegaLambda1998/ShockCooling.jl.git")
