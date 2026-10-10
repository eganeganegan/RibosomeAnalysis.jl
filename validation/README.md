# Public-data and independent validation

Core tests remain offline. This directory contains separately runnable checks
that use public RCSB structures and an independently generated MDAnalysis
reference.

## Public structures

- [1J7T](https://www.rcsb.org/structure/1J7T): paromomycin bound to an
  *E. coli* 16S rRNA A-site model, X-ray diffraction, 2.5 Å resolution.
- [4V7U](https://www.rcsb.org/structure/4V7U): *E. coli* 70S ribosome bound to
  erythromycin, X-ray diffraction, 3.1 Å resolution.

Download and validate exact RCSB files:

```bash
julia --project=. validation/download_public_structures.jl validation/data
julia --project=. validation/validate_public_structures.jl validation/data
```

`public_structures.toml` records the source URLs and SHA-256 checksums. The
downloaded data directory is excluded from version control.

## Independent MDAnalysis reference

The checked-in reference was generated with MDAnalysis 2.10.0 and NumPy 2.4.6.
It independently evaluates a known rigid-body/deformation trajectory and the
5 Å paromomycin contact set in the official PDB-format 1J7T entry.

To reproduce it in a separate Python environment:

```bash
python -m venv .validation-venv
.validation-venv/bin/pip install 'MDAnalysis==2.10.0' 'numpy==2.4.6'
gzip -dk validation/data/1J7T.pdb.gz
.validation-venv/bin/python validation/generate_mdanalysis_reference.py \
  validation/data/1J7T.pdb
```

The reference and Julia calculations agree within `5e-8` Å for RMSD, RMSF,
and the selected distance series. For 1J7T at a 5 Å cutoff, both implementations
report 356 atom-pair contacts, 88 contacting non-ligand atoms, 18 contacting
residues, and a minimum distance of approximately 2.51381 Å.

## Full-ribosome result

The 4V7U integration run loads 284,525 atoms in 104 author chains and identifies
ERY at author chain BA, residue 3136. At 5 Å it reports 283 atom-pair contacts,
96 contacting atoms, and 12 contacting residues. These are geometric results,
not affinity or energetic claims.
