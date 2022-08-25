# DelayDiffEq: Generating Delay Differential Equation Solvers via Recursive Embedding of Ordinary Differential Equation Solvers

[![Build Status](https://github.com/SciML/DelayDiffEq.jl/workflows/CI/badge.svg?branch=paper)](https://github.com/SciML/DelayDiffEq.jl/actions?query=workflow%3ACI+branch%3Apaper)
[![Codecov](https://codecov.io/gh/SciML/DelayDiffEq.jl/branch/paper/graph/badge.svg)](https://codecov.io/gh/SciML/DelayDiffEq.jl)
[![Coveralls](https://coveralls.io/repos/github/SciML/DelayDiffEq.jl/badge.svg?branch=paper)](https://coveralls.io/github/SciML/DelayDiffEq.jl?branch=paper)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

This branch contains the paper about DelayDiffEq and its Julia code.
The LaTeX source code of the paper is contained in the folder `/paper`.
The directory `/src` contains the Julia code, with accompanying tests in `/test`.

## Build document

To generate the PDF version of the paper, the following dependencies need to be installed:
- [Ruby](https://www.ruby-lang.org/en/)
- [latexmk](https://www.ctan.org/pkg/latexmk/)

Navigate to the folder `/paper` and run:

```shell
latexmk -bibtex -pdf paper.tex
```

This will re-generate the PDF file `paper.pdf` in the same folder.

## Reproduce results

All results in the paper were computed with Julia 1.7.3:

```julia
julia> versioninfo()
Julia Version 1.7.3
Commit 742b9abb4d (2022-05-06 12:58 UTC)
Platform Info:
  OS: Linux (x86_64-pc-linux-gnu)
  CPU: 11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-12.0.1 (ORCJIT, tigerlake)
Environment:
  JULIA_NUM_THREADS = 3
```

Install Julia 1.7.3 with [`juliaup`](https://github.com/JuliaLang/juliaup).
Open a terminal in the root folder of this repository and start Julia with

```shell
$ julia +1.7.3 --project=. --startup-file=no
```

Install the required Julia packages (their versions are saved in `Manifest.toml`):

```julia
julia> using Pkg; Pkg.instantiate()
```

Load the Julia package for this paper:

```julia
julia> using DelayDiffEqPaper
```

### Figures

The data for Figures 1 and 2 in the folder `/data` (`/data/mackey_glass.csv` and `/data/waltman.csv`) can be reproduced exactly modulo hardware-dependent floating point differences.
The benchmark results in `/data/benchmark.csv`, shown in Figure 3, are system-dependent and hence can be reproduced only qualitatively but not quantitatively.

In the open Julia REPL run:

```julia
julia> mackey_glass("data/mackey_glass.csv"); # Figure 1

julia> waltman("data/waltman.csv"); # Figure 2

julia> benchmark("data/benchmark.csv"); # Figure 3 (results are system-dependent!)
```

Afterwards you can recompile the PDF version of the paper as described above.

### Anderson acceleration

You can also inspect the effect of Anderson acceleration mentioned in Section 4.1.

In the open Julia REPL run:

```julia
julia> (; default, anderson) = compare_anderson();

julia> default.destats
DelayDiffEq.DDEStats
Number of function 1 evaluations:                            5517
Number of function 2 evaluations:                            0
Number of W matrix evaluations:                              0
Number of linear solves:                                     0
Number of Jacobians created:                                 0
Number of nonlinear solver iterations:                       0
Number of nonlinear solver convergence failures:             0
Number of fixed-point solver iterations:                     720
Number of fixed-point solver convergence failures:           58
Number of rootfind condition calls:                          0
Number of accepted steps:                                    129
Number of rejected steps:                                    8

julia> anderson.destats
DelayDiffEq.DDEStats
Number of function 1 evaluations:                            2829
Number of function 2 evaluations:                            0
Number of W matrix evaluations:                              0
Number of linear solves:                                     124
Number of Jacobians created:                                 0
Number of nonlinear solver iterations:                       0
Number of nonlinear solver convergence failures:             0
Number of fixed-point solver iterations:                     308
Number of fixed-point solver convergence failures:           30
Number of rootfind condition calls:                          0
Number of accepted steps:                                    121
Number of rejected steps:                                    8
```

### Sensitivity example

The sensitivity example in Section 5.2 can be reproduced as well.
You can inspect the Jacobians computed by [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) and [Zygote](https://github.com/FluxML/Zygote.jl).
In the open Julia REPL execute:

```julia
julia> (; forwarddiff, zygote) = jacobian_lv();

julia>  forwarddiff ≈ zygote
true
```

### Tests

The folder `/test` contains tests of the Julia code.
You can run them in the open Julia REPL with:

```julia
julia> using Pkg; Pkg.test()
...
```
