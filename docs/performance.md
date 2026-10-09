# Performance considerations

Binding-site and contact searches use three-dimensional KD trees, avoiding a
full ribosome-by-ligand distance matrix. Atom positions use fixed-size vectors;
bulk numerical operations use contiguous `3 × N` matrices. Contact-persistence
analysis builds one tree per frame and can work on a sliced trajectory.

The MDAnalysis reader accepts `start`, `stop`, `step`, and `selection` so users
can limit data before allocating a Julia array. The present core
`TrajectoryData` is in memory; true one-pass streaming statistics remain future
work. For production ribosome trajectories, select only scientifically needed
atoms, process replicas independently, and test convergence before increasing
sampling density.

Run controlled benchmarks with:

```bash
julia --project=. benchmarks/benchmarks.jl
```

The bundled fixture has only 14 atoms. Its timings are regression smoke tests,
not evidence of performance on full ribosomes. Benchmark real workloads locally
and record hardware, Julia version, package manifest, atom/frame counts, and
selection definitions. Profile allocations before changing algorithms.
