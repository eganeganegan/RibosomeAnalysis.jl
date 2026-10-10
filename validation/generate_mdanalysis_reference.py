import argparse
import json

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align, rms
from MDAnalysis.lib.distances import distance_array


def trajectory_reference():
    reference = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    translated = reference @ rotation.T + np.array([2.0, 3.0, 4.0])
    deformed = translated.copy()
    deformed[3] += np.array([0.1, 0.0, 0.0])
    coordinates = np.stack([reference, translated, deformed])
    rmsd_values = [rms.rmsd(reference, frame, center=True, superposition=True) for frame in coordinates]

    universe = mda.Universe.empty(4, trajectory=True)
    universe.load_new(coordinates.copy(), order="fac")
    align.AlignTraj(universe, universe, select="all", in_memory=True, ref_frame=0).run()
    fitted = np.stack(
        [universe.trajectory[index].positions.copy() for index in range(len(coordinates))]
    )
    rmsf_values = rms.RMSF(universe.atoms).run().results.rmsf
    distances = [float(distance_array(frame[:1], frame[1:2])[0, 0]) for frame in fitted]
    return {
        "rmsd_angstrom": rmsd_values,
        "rmsf_angstrom": rmsf_values.tolist(),
        "distance_1_2_angstrom": distances,
    }


def structure_reference(path):
    universe = mda.Universe(path)
    ligand = universe.select_atoms("resname PAR and chainID A and resid 45")
    other = universe.atoms.difference(ligand)
    distances = distance_array(other.positions, ligand.positions)
    contacting_atoms = (distances <= 5.0).any(axis=1)
    return {
        "format": "PDB",
        "ligand_selection": "resname PAR and chainID A and resid 45",
        "cutoff_angstrom": 5.0,
        "atoms": len(universe.atoms),
        "ligand_atoms": len(ligand),
        "contact_pairs": int((distances <= 5.0).sum()),
        "contacting_atoms": int(contacting_atoms.sum()),
        "contacting_residues": len(other[contacting_atoms].residues),
        "minimum_distance_angstrom": float(distances.min()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="decompressed official RCSB 1J7T.pdb")
    arguments = parser.parse_args()
    print(
        json.dumps(
            {
                "mdanalysis_version": mda.__version__,
                "numpy_version": np.__version__,
                "trajectory": trajectory_reference(),
                "structure_1j7t": structure_reference(arguments.pdb),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
