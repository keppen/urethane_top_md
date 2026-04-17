from copy import copy
from pathlib import Path
from numpy import int_, zeros
import numpy as np
from numpy.typing import NDArray
import parmed as pmd
from parmed import Atom, AtomList, Structure, modeller
from parmed.modeller.residue import ResidueTemplate
from graph import Graph
from collections import deque


class PolymerManager:
    """
    1. Find all urethane groups
    2.
    """

    PROJECT_ROOT = Path(__file__).resolve().parent
    RESIDUE_PATH: Path = Path(f"{PROJECT_ROOT}/residues/residues")

    def __init__(self, structure: pmd.Structure, resname: str | None = None) -> None:
        """
        Manipulate a single strucutre.
        Load starting structure and convet to parmed.Structure object.
        """
        self.residue_templates: dict[str, ResidueTemplate] = self._init_residues()

        self.structure: Structure = structure
        self.urethane_nitrogen_atoms: list[Atom] = [
            atom for atom in self.structure.atoms if self._is_urethane_nitrogen(atom)
        ]
        self.carbonyl_carbon_atoms: list[Atom] = [
            atom for atom in self.structure.atoms if self._is_carbonyl_carbon(atom)
        ]

        self.polymer_mer_list: list[ResidueTemplate] = self._init_mers()
        self.backbone: list[Atom] = self._init_backbone()

        new_structure = self.generate_new_structure()
        # new_structure.save("new_str.pdb", overwrite=True)
        # print(new_structure.atoms)
        self.atom_order_idx = self.get_atom_order(new_structure)

        self.print_backbone()

    def __repr__(self) -> str:
        return str(self.structure)

    def print_backbone(self):
        print(
            "Backbone indexes: VMD visualization, no HN and O atom types",
            " ".join([str(a.idx) for a in self.backbone]),
        )
        print(
            "Backbone indexes: GMX clustering, no HN and O atom types",
            " ".join([str(a.idx + 1) for a in self.backbone]),
        )
        backbone = copy(self.backbone)
        for a in self.backbone:
            bond_elements = [b.element_name for b in a.bond_partners]
            if self._is_carbonyl_carbon(a):
                backbone.extend([b for b in a.bond_partners if b.element_name == "O"])
            if self._is_terminal_oxygen(a):
                backbone.extend([b for b in a.bond_partners if b.element_name == "H"])
            if self._is_urethane_nitrogen(a):
                backbone.extend([b for b in a.bond_partners if b.element_name == "H"])
            if (
                bond_elements.count("C") == 2
                and bond_elements.count("N") == 1
                and bond_elements.count("H") == 1
            ):
                backbone.extend([b for b in a.bond_partners if b.element_name == "C"])
            if (
                bond_elements.count("C") == 2
                and bond_elements.count("O") == 1
                and bond_elements.count("H") == 1
            ):
                backbone.extend([b for b in a.bond_partners if b.element_name == "C"])

        print(
            "Backbone indexes: VMD visualization, with HN and O atom types",
            " ".join([str(a.idx) for a in set(backbone)]),
        )
        print(
            "Backbone indexes: GMX clustering, no HN and O atom types",
            " ".join([str(a.idx + 1) for a in set(backbone)]),
        )
        return " ".join([str(a.idx + 1) for a in set(backbone)])

    def get_structure_params(self):
        new_structure = self.generate_named_structure()
        return self.copy_parameters(new_structure)

    def get_structure_noparams(self):
        # n0 = self.generate_new_structure()
        # n1 = self.generate_new_structure()
        new_structure = self.generate_named_structure()
        # n0.atoms.claim()
        # n1.atoms.claim()
        new_structure.atoms.claim()

        new_order = {v: k for k, v in self.atom_order_idx.items()}

        new_structure.atoms.sort(key=lambda x: new_order[x.idx])
        new_structure.atoms.claim()

        # n0.atoms.sort(key=lambda x: new_order[x.idx])
        # n0.atoms.claim()
        #
        # print("src idx\tdest idx\tsrc_at\tto_sort\tresult")
        # for i in range(len(self.structure.atoms)):
        #     print(
        #         i,
        #         self.atom_order_idx[i],
        #         self.structure.atoms[i].name,
        #         n1.atoms[i].name,
        #         n0.atoms[i].name,
        #         # new_structure.atoms[i],
        #         sep="\t",
        #     )

        for i, atom in enumerate(new_structure.atoms, start=1):
            atom._idx = i - 1
            atom.number = i

        new_structure.coordinates = self.structure.coordinates

        return new_structure

    def extract_mers(self, resname) -> None:
        for i, res in enumerate(self.polymer_mer_list):
            if i == 0:
                res.save(f"Boc_res{resname}.pdb", overwrite=True)
            if i == 1:
                res.save(f"res{resname}.pdb", overwrite=True)
            if i == len(self.polymer_mer_list) - 1:
                res.save(f"Cterm_res{resname}.pdb", overwrite=True)
        # res.save(f"mer{i}.pdb", overwrite=True)

    # --- Init methods ---

    def _init_residues(self) -> dict[str, ResidueTemplate]:
        residue_templates = {}
        for pdb_file in self.RESIDUE_PATH.iterdir():
            if pdb_file.suffix != ".pdb":
                continue
            print("DEBUG: Residue loaded: ", pdb_file)
            structure = pmd.read_PDB(str(pdb_file))
            residue = structure.residues[0]

            template = ResidueTemplate.from_residue(residue)
            template = self._reorder_by_bond_priority(template)
            residue_templates[residue.name] = template

        return residue_templates

    def _init_mers(self) -> list[ResidueTemplate]:
        tail_mers = self._find_longest_path(self.urethane_nitrogen_atoms)
        graph = Graph(self.structure)
        borders = self.urethane_nitrogen_atoms + self.carbonyl_carbon_atoms

        # Determine mer order
        residue0 = graph.breath_first_search(tail_mers[0], borders=borders)
        first_mer, last_mer = (
            (tail_mers[0], tail_mers[1])
            if self._is_last_mer(residue0)
            else tail_mers[::-1]
        )

        path = self._find_path(first_mer, last_mer)
        graph.reset()

        # Collect mers along the path
        return [
            self._reorder_by_bond_priority(
                graph.breath_first_search(
                    nitrogen, borders=self.urethane_nitrogen_atoms
                )
            )
            for nitrogen in path
            if nitrogen.element_name == "N"
        ]

    def _init_backbone(self) -> list[Atom]:
        backbone_atoms = []
        urethane_nitrogen_names = {a.name for a in self.urethane_nitrogen_atoms}
        carbonyl_carbon_names = {a.name for a in self.carbonyl_carbon_atoms}

        for i, residue in enumerate(self.polymer_mer_list):
            residue_atom_names = {a.name for a in residue.atoms}

            # Determine start/end atoms
            start_name = next(iter(residue_atom_names & urethane_nitrogen_names))
            if i == len(self.polymer_mer_list) - 1:
                ends = [a.name for a in residue.atoms if self._is_terminal_oxygen(a)]
            else:
                ends = list(residue_atom_names & carbonyl_carbon_names)

            print("DEBUG: Starting and ending atom name", start_name, ends)

            # Find and extend path
            for end_name in ends:
                start_idx = [a.name for a in self.structure.atoms].index(start_name)
                end_idx = [a.name for a in self.structure.atoms].index(end_name)
                start = self.structure.atoms[start_idx]
                end = self.structure.atoms[end_idx]
                path = self._find_path(start, end)
                print("DEBUG: found path", path)
                if len(path) > 3:
                    break
            else:
                raise ValueError("Path not found")

            backbone_atoms.extend(path)
            # print("VMD", " ".join([str(a.idx) for a in backbone_atoms]))

        return backbone_atoms

    def _find_longest_path(self, all_ns: list[Atom]) -> list[Atom]:
        # Try to find path to another urethane nitrogen
        n_mers = len(all_ns)
        path_lens: NDArray[int_] = zeros([n_mers, n_mers], dtype=int_)
        for i in range(n_mers):
            start_n = all_ns[i]
            for j in range(n_mers):
                end_n = all_ns[j]
                if end_n == start_n:
                    continue

                path = self._find_path(start_n, end_n)
                if path:
                    path_lens[i, j] = len(path)

        i = path_lens.argmax() % n_mers
        j = path_lens.argmax() // n_mers

        return [all_ns[i], all_ns[j]]

    def _find_path(
        self,
        start: Atom,
        end: Atom,
    ) -> list[Atom]:
        """Find path between atoms matching element pattern using BFS."""
        queue = deque([(start, [start])])

        while queue:
            current, path = queue.popleft()

            if current == end:
                return path

            # Explore neighbors
            for neighbor in current.bond_partners:
                if neighbor in path:
                    continue

                new_path = path + [neighbor]
                queue.append((neighbor, new_path))

        return []

    # --- Topology and atom definition methods ---

    def _is_last_mer(self, residue: ResidueTemplate) -> bool:
        for o_atom in [atom for atom in residue.atoms if atom.element_name == "O"]:
            if self._is_terminal_oxygen(o_atom):
                continue
            path = self._find_path(residue.atoms[0], o_atom)
            path_element = "".join([a.element_name for a in path])
            if path_element == "NCCO" or path_element == "NCCCO":
                return True

        return False

    def _is_terminal_oxygen(self, o_atom: Atom) -> bool:
        if o_atom.element_name != "O":
            return False
        return len(o_atom.bond_partners) == 2 or o_atom.bond_partners.count("H") == 1

    def _is_residue_nitrogen(self, n_atom: Atom) -> bool:
        if n_atom.element_name != "N":
            return False
        bonded_elements = [a.element_name for a in n_atom.bond_partners]
        return bonded_elements.count("C") == 1 and bonded_elements.count("H") == 1

    def _is_urethane_nitrogen(self, n_atom: Atom) -> bool:
        if n_atom.element_name != "N":
            return False
        bonded_elements = [a.element_name for a in n_atom.bond_partners]
        return bonded_elements.count("C") >= 2 and any(
            self._is_carbonyl_carbon(a) for a in n_atom.bond_partners
        )

    def _is_carbonyl_carbon(self, c_atom: Atom) -> bool:
        if c_atom.element_name != "C":
            return False
        bonded_elements = [a.element_name for a in c_atom.bond_partners]
        return bonded_elements.count("O") == 2 and len(c_atom.bond_partners) == 3

    # --- Structure methods ---

    def _reorder_by_bond_priority(
        self, structure: ResidueTemplate | Structure
    ) -> ResidueTemplate:
        """Reorder atoms based on bond priority starting from a nitrogen atom."""
        # Convert to Structure if needed
        if isinstance(structure, ResidueTemplate):
            structure = structure.to_structure()

        # Find nitrogen atom
        nitrogen = next(
            (
                atom
                for atom in structure.atoms
                if self._is_residue_nitrogen(atom) or self._is_urethane_nitrogen(atom)
            ),
            None,
        )

        if nitrogen is None or nitrogen.element_name != "N":
            raise ValueError("Nitrogen not found in residue!")

        graph = Graph(
            structure=structure,
            testing=True,
        )

        iter1 = graph.breath_first_search(nitrogen)
        iter1 = iter1.to_structure()

        # Find nitrogen atom
        nitrogen = next(
            (
                atom
                for atom in iter1.atoms
                if self._is_residue_nitrogen(atom) or self._is_urethane_nitrogen(atom)
            ),
            None,
        )

        if nitrogen is None or nitrogen.element_name != "N":
            raise ValueError("Nitrogen not found in residue!")
        graph = Graph(
            structure=iter1,
            testing=True,
        )
        iter2 = graph.breath_first_search(nitrogen)
        iter2 = iter2.to_structure()

        # Find nitrogen atom
        nitrogen = next(
            (
                atom
                for atom in iter2.atoms
                if self._is_residue_nitrogen(atom) or self._is_urethane_nitrogen(atom)
            ),
            None,
        )

        if nitrogen is None or nitrogen.element_name != "N":
            raise ValueError("Nitrogen not found in residue!")
        graph = Graph(
            structure=iter2,
            testing=True,
        )
        iter3 = graph.breath_first_search(nitrogen)
        return iter3

    def get_atom_order(self, new_structure: Structure) -> dict[int, int]:
        """Get atom index mapping between new structu.re and template."""
        return {
            atom.idx: [a.name for a in new_structure.atoms].index(atom.name)
            for atom in self.structure.atoms
        }

    def generate_named_structure(self) -> Structure:
        """Standardize atom names across all residues."""
        return self._rejoin_residues(
            [self.detect_residue(residue) for residue in self.polymer_mer_list]
        )

    def generate_new_structure(self) -> Structure:
        """Generate complete polymer structure from residue templates."""
        return self._rejoin_residues(
            [residue.to_structure() for residue in self.polymer_mer_list]
        )

    def detect_residue(self, residue: ResidueTemplate) -> Structure:
        """Identify residue type and return standardized structure."""
        for res_name, template in self.residue_templates.items():
            if self._is_same_residue(residue, template):
                structure = template.to_structure()
                structure.residues[0].name = res_name
                return structure
        raise ValueError(f"Residue {residue} not found in templates")

    def _is_same_residue(
        self,
        res1: ResidueTemplate,
        res2: ResidueTemplate,
    ) -> bool:
        """Check if two residues have identical chemical structure."""
        # Quick check using chemical formula
        if res1.empirical_chemical_formula != res2.empirical_chemical_formula:
            return False

        # Convert to structures for detailed comparison
        struct1 = res1.to_structure()
        struct2 = res2.to_structure()

        # print(struct1, struct2)

        # Compare adjacency matrices
        return self._compare_residue_graphs(struct1, struct2)

    def _compare_residue_graphs(self, struct1: Structure, struct2: Structure) -> bool:
        """Compare residues using graph isomorphism. Struc2 is the template!"""
        # Build graphs for both structures
        adj_matrix1 = self._build_adjacency_matrix(struct1.atoms)
        adj_matrix2 = self._build_adjacency_matrix(struct2.atoms)

        # print(adj_matrix1)
        # print(adj_matrix2)

        return np.array_equal(adj_matrix1, adj_matrix2)

    def _build_adjacency_matrix(self, atoms: list[Atom]) -> NDArray[np.int8]:
        """Create adjacency matrix for atom list."""
        n = len(atoms)
        adj_matrix = np.zeros((n, n), dtype=np.int8)
        idx_map = {atom.idx: i for i, atom in enumerate(atoms)}

        for i, atom in enumerate(atoms):
            for neighbor in atom.bond_partners:
                if neighbor.idx in idx_map:  # Only consider neighbors in the atom list
                    j = idx_map[neighbor.idx]
                    adj_matrix[i, j] = 1
                    adj_matrix[j, i] = 1  # Undirected graph
        return adj_matrix

    def _rejoin_residues(self, residue_structures: list[Structure]) -> Structure:
        """Combine multiple residue structures into a single polymer structure."""
        if not residue_structures:
            raise ValueError("Cannot rejoin empty residue list")

        # Start with first residue and accumulate others
        polymer = residue_structures[0].copy(Structure)
        for residue in residue_structures[1:]:
            polymer += residue

        return polymer

    def copy_parameters(self, dest: Structure) -> Structure:
        src = self.structure
        order = self.atom_order_idx
        src_to_dest = {src.atoms[k]: dest.atoms[v] for k, v in order.items()}

        dest.bonds.clear()
        for bond in src.bonds:
            a1 = src_to_dest[bond.atom1]
            a2 = src_to_dest[bond.atom2]
            dest.bonds.append(pmd.Bond(a1, a2, type=bond.type))

        dest.angles.clear()
        for angle in src.angles:
            a1 = src_to_dest[angle.atom1]
            a2 = src_to_dest[angle.atom2]
            a3 = src_to_dest[angle.atom3]
            dest.angles.append(pmd.Angle(a1, a2, a3, type=angle.type))

        dest.dihedrals.clear()
        for dih in src.dihedrals:
            a1 = src_to_dest[dih.atom1]
            a2 = src_to_dest[dih.atom2]
            a3 = src_to_dest[dih.atom3]
            a4 = src_to_dest[dih.atom4]
            dest.dihedrals.append(pmd.Dihedral(a1, a2, a3, a4, type=dih.type))

        dest.impropers.clear()
        for imp in src.impropers:
            a1 = src_to_dest[imp.atom1]
            a2 = src_to_dest[imp.atom2]
            a3 = src_to_dest[imp.atom3]
            a4 = src_to_dest[imp.atom4]
            dest.impropers.append(pmd.Improper(a1, a2, a3, a4, type=imp.type))

        if hasattr(src, "urey_bradleys"):
            dest.urey_bradleys.clear()
            for ub in src.urey_bradleys:
                a1 = src_to_dest[ub.atom1]
                a2 = src_to_dest[ub.atom2]
                dest.urey_bradleys.append(pmd.UreyBradley(a1, a2, type=ub.type))

        if hasattr(src, "cmaps"):
            dest.cmaps.clear()
            dest.cmaps.extend(src.cmaps)

        dest.bond_types.clear()
        dest.bond_types.extend(src.bond_types)
        dest.angle_types.clear()
        dest.angle_types.extend(src.angle_types)
        dest.dihedral_types.clear()
        dest.dihedral_types.extend(src.dihedral_types)
        dest.improper_types.clear()
        dest.improper_types.extend(src.improper_types)
        if hasattr(src, "urey_bradley_types"):
            dest.urey_bradley_types.clear()
            dest.urey_bradley_types.extend(src.urey_bradley_types)
        if hasattr(src, "cmap_types"):
            dest.cmap_types.clear()
            dest.cmap_types.extend(src.cmap_types)

        for a_src, a_dest in src_to_dest.items():
            a_dest.charge = a_src.charge
            a_dest.type = a_src.type
            a_dest.atomic_number = a_src.atomic_number

            a_dest.xx, a_dest.xy, a_dest.xz = a_src.xx, a_src.xy, a_src.xz

            a_dest.nb_idx = a_src.nb_idx
            a_dest.rmin = a_src.rmin
            a_dest.epsilon = a_src.epsilon
            a_dest.rmin14 = getattr(a_src, "rmin14", None)
            a_dest.epsilon14 = getattr(a_src, "epsilon14", None)

        dest.box = src.box

        return dest


if __name__ == "__main__":
    import sys
    import os

    # Used just to update atom names, residues names and ids.
    # Coordinates are retained, but I do not gaurantee if other properties stay the same

    # To create new atom group in ndx file
    update_ndx_file = False
    # PDB file HAS to have each atom named UNIQULY, ParmEd limitation
    pdb = sys.argv[1]

    pdb_path = Path(pdb)
    basename_path = pdb_path.stem
    structure = pmd.read_PDB(pdb)
    polymer = PolymerManager(structure)
    os.makedirs(f"{basename_path}", exist_ok=True)

    ######
    # Topology and configuration generator using ParmEd
    # Limited up to no use - parmed limitations
    # gro - per residue sorted, seems ok, not tested
    # top - failure
    # pdb - per residue sorted, seems ok, not tested
    poly1 = polymer.get_structure_params()
    poly1.save(f"{basename_path + '/params.top'}", format="gromacs", overwrite=True)
    poly1.save(f"{basename_path + '/params.gro'}", format="gro", overwrite=True)
    poly1.save(
        f"{basename_path + '/params.pdb'}", format="pdb", overwrite=True, renumber=True
    )
    #####
    # Polymer object with atom sorted by hand. No topolgy generating
    # pdb - per residue sorted, wired atom indexing
    # gro - per residue sorted, indexing seems ok, yet different from poly1 gro
    poly2 = polymer.get_structure_noparams()
    poly2.save(
        f"{basename_path + '/noparams.pdb'}",
        format="pdb",
        overwrite=True,
        renumber=False,
    )
    poly2.save(f"{basename_path + '/noparams.gro'}", format="gro", overwrite=True)
    #####

    #####
    # Simple atom name, resid and resname substitution
    # Actually used solution
    poly3 = polymer.get_structure_noparams()
    new_name = Path(pdb.replace("_NEW", "")).stem
    with open(f"{basename_path}/NAMED-{new_name}.pdb", "w") as f:
        # sort your atoms globally however you like:
        sorted_atoms = sorted(poly3.atoms, key=lambda a: a._idx)
        for i, atom in enumerate(sorted_atoms, start=1):
            x, y, z = atom.xx, atom.xy, atom.xz
            f.write(
                f"HETATM{i:5d} {atom.name:^4s}{atom.residue.name:>3s} "
                f"{atom.residue.number:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                "  1.00  0.00           "
                f"{atom.element_name:>2s}\n"
            )
    #####

    #####
    # Update GROMACS index file for clustering
    if update_ndx_file:
        backbone_ndx = polymer.print_backbone()
        index_file = sys.argv[2]
        with open(index_file, "a") as file:
            title_line = "[ BACKBONE+CZ ]\n"
            file.write(title_line)
            file.write(f"{backbone_ndx}\n")
