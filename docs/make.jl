using Documenter
using HyperCloud

makedocs(;
    modules=[HyperCloud],
    authors="Pascal Jung", # Update with your desired author name
    sitename="HyperCloud.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://blhackslash.github.io/HyperCloud.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "PDEs" => [
            "Architecture" => "PDEs/main.md",
            "Linear Advection" => "PDEs/advection.md",
            "Burgers" => "PDEs/burgers.md",
            "Euler" => "PDEs/euler.md",
        ],
        "Source Terms" => [
            "Kinetic Relaxation" => "SourceTerms/kinetic.md",
        ],
        "Initial Conditions" => "ICs.md",
        "Simulation Engine" => "simulation.md",
        "Time Integration" => "time_integration.md",
        "Execution Macros" => "macros.md",
        "Convenience Wrappers" => "wrappers.md",
    ],
)

deploydocs(;
    repo="github.com/blhackslash/HyperCloud.jl",
    devbranch="main",
)