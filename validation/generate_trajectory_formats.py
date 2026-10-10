"""Generate deterministic XTC/TRR/DCD fixtures for trajectory validation."""

import argparse
import json
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import calc_bonds


COORDINATES_ANGSTROM = np.array(
    [
        [[0.5, 0.5, 0.5], [9.7, 0.5, 0.5], [2.0, 2.0, 2.0], [2.8, 2.0, 2.0]],
        [[0.5, 0.5, 0.5], [9.4, 0.5, 0.5], [2.1, 2.0, 2.0], [2.9, 2.0, 2.0]],
        [[0.5, 0.5, 0.5], [9.8, 0.5, 0.5], [2.2, 2.0, 2.0], [3.0, 2.0, 2.0]],
        [[0.5, 0.5, 0.5], [9.2, 0.5, 0.5], [2.3, 2.0, 2.0], [3.1, 2.0, 2.0]],
        [[0.5, 0.5, 0.5], [9.9, 0.5, 0.5], [2.4, 2.0, 2.0], [3.2, 2.0, 2.0]],
    ],
    dtype=np.float32,
)

DIMENSIONS_ANGSTROM_DEGREES = np.tile(
    np.array([10.0, 10.0, 10.0, 90.0, 90.0, 90.0], dtype=np.float32),
    (len(COORDINATES_ANGSTROM), 1),
)


def write_topology(path: Path, atom_count: int = 4) -> None:
    lines = [
        "CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1",
    ]
    names = ["C1", "C2", "N1", "O1"]
    for index, (name, coordinate) in enumerate(
        zip(names[:atom_count], COORDINATES_ANGSTROM[0, :atom_count]), start=1
    ):
        x, y, z = coordinate
        element = name[0]
        lines.append(
            f"ATOM  {index:5d} {name:^4s} MOL A   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
        )
    lines.extend(["TER", "END"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    arguments.output.mkdir(parents=True, exist_ok=True)

    topology = arguments.output / "topology.pdb"
    write_topology(topology)
    write_topology(arguments.output / "wrong-topology.pdb", atom_count=3)

    universe = mda.Universe(topology)
    universe.load_new(
        COORDINATES_ANGSTROM.copy(),
        order="fac",
        dimensions=DIMENSIONS_ANGSTROM_DEGREES.copy(),
        dt=2.0,
    )
    for suffix in ("xtc", "trr", "dcd"):
        with mda.Writer(
            str(arguments.output / f"trajectory.{suffix}"),
            n_atoms=len(universe.atoms),
            dt=2.0,
        ) as writer:
            for _ in universe.trajectory:
                writer.write(universe.atoms)

    periodic_distances = []
    for frame, dimensions in zip(COORDINATES_ANGSTROM, DIMENSIONS_ANGSTROM_DEGREES):
        periodic_distances.append(
            float(calc_bonds(frame[0], frame[1], box=dimensions))
        )
    reference = {
        "mdanalysis_version": mda.__version__,
        "numpy_version": np.__version__,
        "coordinates_angstrom": COORDINATES_ANGSTROM.tolist(),
        "dimensions_angstrom_degrees": DIMENSIONS_ANGSTROM_DEGREES.tolist(),
        "times_ps": [2.0 * index for index in range(len(COORDINATES_ANGSTROM))],
        "periodic_distance_1_2_angstrom": periodic_distances,
    }
    (arguments.output / "reference.json").write_text(
        json.dumps(reference, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
