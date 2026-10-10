# Performance considerations

Binding-site and contact searches use three-dimensional KD trees, avoiding a
full ribosome-by-ligand distance matrix. Atom positions use fixed-size vectors;
bulk numerical operations use contiguous `3 × N` matrices. Contact-persistence
analysis builds one tree per frame and can work on a sliced trajectory.

The MDAnalysis reader accepts `start`, `stop`, `step`, and `selection` so users
can limit data before allocating a Julia array. `TrajectoryData` is in memory;
`foreach_trajectory_frame` instead keeps only the current selected frame in
Julia memory. For production ribosome trajectories, select only scientifically
needed atoms, process replicas independently, and test convergence before
increasing sampling density.

Nonperiodic contact persistence uses a KD tree per frame. Periodic contact
persistence currently checks selected atom pairs directly because the KD-tree
backend has no triclinic periodic metric. Keep periodic groups focused and
benchmark representative selections before scaling to whole ribosomes.

Run controlled benchmarks with:

```bash
julia --project=. benchmarks/benchmarks.jl
```

The bundled fixture has only 14 atoms. Its timings are regression smoke tests,
not evidence of performance on full ribosomes. Benchmark real workloads locally
and record hardware, Julia version, package manifest, atom/frame counts, and
selection definitions. Profile allocations before changing algorithms.

## Recorded smoke benchmark

The following minimum trial times were measured on 2026-10-09 with Julia 1.13.1
on macOS arm64 using the committed manifest and the 14-atom controlled fixture:

| Operation | Minimum trial time |
| --- | ---: |
| Structure loading | 61.417 μs |
| 5 Å binding-site search | 6.583 μs |
| 4 Å atomic contacts | 2.875 μs |
| Four-atom Kabsch fit | 2.667 μs |
| RMSD over 100 four-atom frames | 311.375 μs |

These are actual development-run measurements, not portable performance
claims. The public-data validation additionally completed its 4V7U assertions
for a 284,525-atom structure in 8.0 seconds after package startup on the same
machine; that test combines parsing, selection, binding-site, and contact work
and is not a microbenchmark.
