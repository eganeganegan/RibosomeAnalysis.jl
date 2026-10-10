# Scientific methodology

## Distances and contacts

For ligand atoms `L`, cutoff `r`, and residue `i`, the binding site is

`B(L,r) = { i : min(distance(atom in i, atom in L)) <= r }`.

The minimum uses heavy atoms by default. The software supports 3, 4, 5, and 8 Å
as useful reporting thresholds, but none is a universal chemical definition.
Contact persistence is the fraction of sampled frames satisfying the chosen
cutoff. It is not a probability of binding.

## Interaction geometry

Atomic contacts need coordinates only. Hydrogen bonds require explicit hydrogen
coordinates and an explicit covalent bond list, then apply donor–acceptor
distance and D–H–A angle thresholds. Protonation, tautomer state, water
competition, and interaction energy remain unresolved. Ionic and aromatic
analysis accepts explicit chemical groups rather than guessing them from atom
names. Metal-mediated results are candidate geometries, not coordination bonds.

## Alignment and RMSD

Needleman–Wunsch sequence alignment establishes candidate residue
correspondence within user-confirmed chain pairs. Matching backbone atom names
are fitted with a proper-rotation Kabsch solution. The reported RMSD is

`sqrt(sum(||R x_i + t - y_i||^2) / N)`.

Coverage and unmatched residues accompany each fit. Equal residue numbers are
never the sole basis for correspondence. Cross-species homology needs external
biological validation.

## Geometry

Centroid, radius of gyration, axis-aligned extent, pairwise-distance summaries,
and local contact density describe selected coordinates. These are not an
experimentally validated pocket volume, solvent accessibility, affinity, or
free energy.

## Dynamics and statistics

RMSD frames are aligned before measurement. RMSF uses the aligned ensemble
mean. PCA is singular-value decomposition of mean-centered Cartesian features;
components are mathematical variance directions and receive no automatic
physical interpretation.

Periodic distances convert a Cartesian displacement to lattice coordinates and
search the 27 translations surrounding the nearest lattice-coordinate image.
This matches minimum-image distances for orthorhombic and reduced triclinic MD
cells. Periodic contact persistence applies that distance independently to each
atom pair and frame. Periodic handling is never enabled implicitly. Molecular
unwrapping is topology dependent and remains the caller's responsibility before
alignment-based observables.

Ordinary bootstrap intervals assume independent observations. Moving-block
bootstrap intervals retain within-block time correlation. Effective sample size
uses the positive initial autocorrelation sequence. Replica comparisons first
reduce each independent replica to a summary, preventing frames from being
misrepresented as independent simulations.

## Physical units

All stored structure and trajectory distances are Angstrom. MDAnalysis reports
positions in Angstrom under its standard conversion behavior. Array input
declared as nanometers is multiplied by 10. Times are picoseconds. Every figure
labels units explicitly. Stored unit-cell vectors receive the same length-unit
conversion as coordinates.
