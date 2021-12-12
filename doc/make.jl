using Documenter, SpinEnsembles

Manual = [
    "Function" => "manual/functions.md",
    "Test" => "manual/manual.md",
]

makedocs(
    sitename = "SpinEnsembles",
    authors = "Yuning Zhang",
    pages = [
        "Home" => "index.md",
        "Manual" => Manual,
])

deploydocs(repo="https://github.com/Neuromancer43/SpinEnsembles.jl.git",
target = "build",
push_preview = true)