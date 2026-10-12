# Installation

RibosomeAnalysis supports Julia 1.10 and later. Install the released source
directly from GitHub while initial Julia General registration is pending:

```julia
using Pkg
Pkg.add(url="https://github.com/eganeganegan/RibosomeAnalysis.jl", rev="main")
```

For development, clone the repository and instantiate the package environment:

```bash
git clone https://github.com/eganeganegan/RibosomeAnalysis.jl.git
cd RibosomeAnalysis.jl
```

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

As is customary for a library package, the root `Manifest.toml` is not
committed. Julia resolves dependency versions compatible with `Project.toml`;
developers can keep a local manifest for reproducible work.

## Optional plotting

The core package does not install or load Makie. Add and load CairoMakie to
activate deterministic static rendering for PNG, SVG, and PDF:

```julia
using Pkg
Pkg.add("CairoMakie")
using CairoMakie, RibosomeAnalysis
```

## Optional trajectory formats

Static structure analysis and array-based trajectory analysis require no
Python. XTC, TRR, DCD, and other formats supported by MDAnalysis require both
PythonCall.jl and a Python environment containing MDAnalysis:

```julia
using Pkg
Pkg.add("PythonCall")
using PythonCall
CondaPkg = pyimport("conda") # use your established PythonCall environment policy
```

Install MDAnalysis according to PythonCall/CondaPkg policy at your institution.
The package never silently installs it. Confirm topology/trajectory atom counts,
selection ordering, time units, and periodic-boundary handling for each format.

GLMakie is not required.
