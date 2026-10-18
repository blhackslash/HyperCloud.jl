using Documenter
using DocumenterInterLinks
using HyperCloud

# Cross-references into the numerical core, e.g. [`reserve!`](@extref HyperCloudCore).
# HYPERCLOUDCORE_INVENTORY overrides the inventory, e.g. with a local HyperCloudCore.jl/docs/build/objects.inv
# to check references to names that are not deployed yet.
links = InterLinks(
    "HyperCloudCore" => (
        "https://blhackslash.github.io/HyperCloudCore.jl/dev/",
        get(ENV, "HYPERCLOUDCORE_INVENTORY", "https://blhackslash.github.io/HyperCloudCore.jl/dev/objects.inv"),
    ),
)

makedocs(;
    modules=[HyperCloud],
    # Explicit, so the build does not depend on the `origin` remote (Forgejo checkouts).
    repo=Remotes.GitHub("blhackslash", "HyperCloud.jl"),
    plugins=[links],
    authors="Pascal Jung",
    sitename="HyperCloud.jl",
    checkdocs = :exports,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://blhackslash.github.io/HyperCloud.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Configuration Reference" => "configuration.md",
        "Execution Macros" => "macros.md",
        "Simulation Engine" => "simulation.md",
        "Time Integration" => "time_integration.md",
        "Initial Conditions" => "ICs.md",
        "PDEs" => [
            "Architecture" => "PDEs/main.md",
            "Linear Advection" => "PDEs/advection.md",
            "Burgers" => "PDEs/burgers.md",
            "Euler" => "PDEs/euler.md",
        ],
        "Source Terms" => [
            "Kinetic Relaxation" => "SourceTerms/kinetic.md",
        ],
        "Convenience Wrappers" => "wrappers.md",
        "Builder API" => "builder.md",
    ],
)

deploydocs(;
    repo="github.com/blhackslash/HyperCloud.jl",
    devbranch="main",
    push_preview = true,
)
