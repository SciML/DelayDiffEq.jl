using Documenter, DelayDiffEq
include("pages.jl")

makedocs(
    sitename = "DelayDiffEq.jl",
    authors = "Chris Rackauckas et al.",
    clean = true,
    doctest = false,
    modules = [DelayDiffEq],
    warnonly = [:docs_block, :missing_docs, :eval_block],
    format = Documenter.HTML(
        analytics = "UA-90474609-3",
        canonical = "https://delaydiffeq.sciml.ai/stable/"
    ),
    pages = pages
)

deploydocs(
    repo = "github.com/SciML/DelayDiffEq.jl";
    push_preview = true
)
