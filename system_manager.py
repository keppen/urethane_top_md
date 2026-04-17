from pathlib import Path
import parmed as pmd
from parmed import Structure
from polymer_manager import PolymerManager

# Supported file extensions
AMBER_EXT = [".rst7", ".parm7", ".nc"]
GROMACS_EXT = [".gro", ".top", ".xtc", ".trr"]


def load_structure(coordinate_file: Path, parameter_file: Path) -> Structure:
    """Load coordinates and parameters into a single Structure."""
    structure = _load_parameters(parameter_file)
    coords, box = _load_coordinates(coordinate_file)
    structure.positions, structure.box = coords, box
    return structure


def _load_coordinates(file: Path):
    if file.suffix in AMBER_EXT:
        coords = pmd.amber.Rst7(str(file))
    elif file.suffix in GROMACS_EXT:
        coords = pmd.gromacs.GromacsGroFile.parse(str(file))
    else:
        raise KeyError(f"Unknown extension: {file.suffix}")
    if coords.positions is None or coords.box is None:
        raise ValueError(f"Incomplete data in {file}")
    return coords.positions, coords.box


def _load_parameters(file: Path) -> Structure:
    if file.suffix in AMBER_EXT:
        parm = pmd.amber.AmberParm(str(file))
    elif file.suffix in GROMACS_EXT:
        parm = pmd.gromacs.GromacsTopologyFile(str(file))
    else:
        raise KeyError(f"Unknown extension: {file.suffix}")
    return parm.copy(Structure)


class ManageSystem:
    """Top-level system handler"""

    def __init__(self, coord: Path, parm: Path):
        self.system = load_structure(coord, parm)

    def process_first_monomer(self) -> None:
        """Extract first monomer and process its polymer chain."""
        monomer = self.system.split()[0][0]
        polymer = PolymerManager(monomer)


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <coordinate_file> <parameter_file>")
        sys.exit(1)
    coord_file = Path(sys.argv[1])
    parm_file = Path(sys.argv[2])
    ms = ManageSystem(coord_file, parm_file)
    ms.process_first_monomer()
