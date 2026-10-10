# Installation

RibosomeAnalysis supports Julia 1.10 and later. Install the released source
directly from GitHub while initial Julia General registration is pending:

```julia
using Pkg
Pkg.add(url="https://github.com/eganeganegan/RibosomeAnalysis.jl", rev="v0.1.0")
```

For development, clone the repository and instantiate its checked-in
environment:

```bash
git clone https://github.com/eganeganegan/RibosomeAnalysis.jl.git
cd RibosomeAnalysis.jl
```

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

The checked-in `Manifest.toml` records the environment used for this release.
Julia may resolve platform-specific binary artifacts on another operating
system.

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

GLMakie is not required. CairoMakie supplies deterministic static rendering for
PNG, SVG, and PDF.
