import numpy as np
import hashlib
from numpy.typing import NDArray
import parmed as pmd
from parmed import Atom, AtomList, Structure
from parmed.modeller import ResidueTemplate
from queue import Queue


class Graph:
    def __init__(
        self,
        structure: Structure,
        visited_atoms: list[int] = [],
        testing=False,
    ) -> None:
        self.structure: Structure = structure
        self.structure_copy: Structure = structure.copy(Structure)

        # Use WL to get labels
        wl_labels = self.weisfeiler_lehman_labels(10)
        self.bond_priority = {}
        for atom in self.structure.atoms:
            partners = atom.bond_partners
            sorted_partners = sorted(partners, key=lambda a: wl_labels[a.idx])

            self.bond_priority[atom] = tuple(sorted_partners)

        for bond in self.structure_copy.bonds:
            bond.delete()

        self.visited_atoms: NDArray[np.bool_] = np.zeros(
            [len(self.structure.atoms)], dtype=bool
        )
        for atom_ndx in visited_atoms:
            self.visited_atoms[atom_ndx] = True

        self.__visited_atoms_bakup = np.copy(self.visited_atoms)

    def reset(self):
        self.structure_copy: Structure = self.structure.copy(Structure)
        for bond in self.structure_copy.bonds:
            bond.delete()
        self.visited_atoms = np.copy(self.__visited_atoms_bakup)

    def breath_first_search(
        self,
        start_atom: Atom,
        borders: list[Atom] | None = None,
        depth: int = -1,
    ) -> ResidueTemplate:
        neighbor: Atom
        neighbor_copy: Atom
        current_atom: Atom
        current_atom_copy: Atom

        residue = ResidueTemplate("UNK")
        residue.add_atom(self.structure_copy.atoms[start_atom.idx])

        path = [start_atom]
        atom_path: tuple[Atom, list[Atom]] = (start_atom, path)
        queue: Queue[tuple[Atom, list[Atom]]] = Queue()
        queue.put(atom_path)

        while not queue.empty():
            current_atom, path = queue.get()
            current_atom_copy = self.structure_copy.atoms[current_atom.idx]
            self.visited_atoms[current_atom.idx] = True

            neighbors = self.bond_priority[current_atom]

            for neighbor in neighbors:
                neighbor_copy = self.structure_copy.atoms[neighbor.idx]
                atom_path = (neighbor, path + [neighbor])

                if self.visited_atoms[neighbor.idx]:
                    continue
                if depth >= 0 and len(atom_path[1]) > depth:
                    self.visited_atoms[neighbor.idx] = True
                    continue
                if borders and neighbor in borders:
                    self.visited_atoms[neighbor.idx] = True
                    continue

                self.visited_atoms[neighbor.idx] = True
                self.update_residuetemplate(residue, current_atom_copy, neighbor_copy)

                queue.put(atom_path)
        return residue

    def update_residuetemplate(
        self, res: ResidueTemplate, atom: Atom, bond_to: Atom
    ) -> None:
        res.add_atom(bond_to)
        res.add_bond(atom, bond_to)

    def weisfeiler_lehman_labels(self, iterations: int = 3) -> dict[int, str]:
        """Computes Weisfeiler-Lehman labels for a ParmEd structure."""
        # Step 1: Initial labels using atomic number or type
        labels = {atom.idx: str(atom.atomic_number) for atom in self.structure.atoms}

        for _ in range(iterations):
            new_labels = {}
            for atom in self.structure.atoms:
                neighbors = atom.bond_partners
                neighbor_labels = sorted(
                    [labels[neighbor.idx] for neighbor in neighbors]
                )
                combined = labels[atom.idx] + "|" + "|".join(neighbor_labels)
                # Hash the label to keep it compact and consistent
                hash_label = hashlib.sha1(combined.encode()).hexdigest()
                new_labels[atom.idx] = hash_label
                # new_labels[atom.idx] = combined
            labels = new_labels

        return labels  # maps atom.idx → label
