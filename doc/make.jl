using Documenter, SpinControl

Manual = ["Function" => "manual/functions.md", "Test" => "manual/manual.md"]

makedocs(
    sitename = "SpinControl",
    authors = "Yuning Zhang",
    pages = ["Home" => "index.md", "Manual" => Manual],
)

deploydocs(
    repo = "https://github.com/Neuromancer43/SpinControl.jl.git",
    target = "build",
    push_preview = true,
)
