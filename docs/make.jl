using Documenter, SpinEnsembles

makedocs(
    sitename = "SpinEnsembles",
    authors = "Yuning Zhang",
    pages = [
        "Home" => "index.md",
        ])

deploydocs(repo="https://github.com/Neuromancer43/SpinEnsembles.jl.git")
