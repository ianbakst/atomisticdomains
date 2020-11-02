import numpy as np
import time
import os
from atomistic_domains import Atom
from typing import List


class Lattice:
    """
    The class "lattice" defines methods which manipulate lattice structures.
    Lattice structures have the following attributes:
    :param type: Lattice type: e.g. fcc, hcp, bcc, etc.
    :param box: Unit cell box dimensions (including tilt)
    :param basis_atoms: List of atoms (of object class: atom) which make up the
                                                primitive unit cell (in relative coordinates)
    :param chem: List of chemical elements in alloy
    :param num_el: Number of different elements present
    :param num_atoms: number of atoms per each element (dictionary type)
    :param tot_atoms: Total number of atoms in unit cell
    """

    def __init__(
        self,
        x_vector: np.ndarray = np.zeros(shape=(1, 3)),
        y_vector: np.ndarray = np.zeros(shape=(1, 3)),
        z_vector: np.ndarray = np.zeros(shape=(1, 3)),
        basis_atoms: List[Atom] = None,
    ):
        self.x_vector = x_vector
        self.y_vector = y_vector
        self.z_vector = z_vector
        if not basis_atoms:
            basis_atoms = []
        else:
            self.basis_atoms = basis_atoms
        return

    def add_atom(self, new_atom):
        """
        Adds an atom to the list of basis atoms of the lattice. The new atom
        must be of the type: "atom",
        """
        self.basis_atoms.append(new_atom)
        self.tot_atoms += 1
        new = False
        for i in range(self.num_el):
            if new_atom.element == self.chem[i]:
                new = True
        if new == False:
            self.chem.append(new_atom.element)
            self.num_atoms[new_atom.element] = 0
        self.num_atoms[new_atom.element] += 1

        return

    def genfromposcar(self, filename):
        scale = np.genfromtxt(filename, skip_header=1, max_rows=1)
        box = np.genfromtxt(filename, delimiter="", skip_header=2, max_rows=3)
        self.box = box
        self.chem = np.genfromtxt(filename, delimiter="", skip_header=5, max_rows=1, dtype=str)
        self.num_el = len(self.chem)
        atom_numbers = np.genfromtxt(filename, delimiter="", skip_header=6, max_rows=1, dtype=int)
        for i in range(self.num_el):
            self.num_atoms[self.chem[i]] = atom_numbers[i]
        self.tot_atoms = sum(atom_numbers)
        ptype = np.genfromtxt(filename, delimiter=">", skip_header=7, max_rows=1, dtype=str)
        if ptype == "Selective Dynamics":
            eln = 0
            j = 0
            for i in range(self.tot_atoms):
                if j < self.num_atoms[self.chem[eln]]:
                    atomel = self.chem[eln]
                    j += 1
                else:
                    eln += 1
                    atomel = self.chem[eln]
                    j = 1
                acoord = np.genfromtxt(filename, delimiter="", skip_header=9 + i, max_rows=1)[
                    :, 0:3
                ]
                adyn = np.genfromtxt(
                    filename, delimiter="", skip_header=9 + i, max_rows=1, dtype=str
                )[:, 3:]
                newatom = atom(id=i + 1, element=atomel, pos=acoord, dyn=adyn)
                self.basis_atoms.append(newatom)
        else:
            eln = 0
            j = 0
            for i in range(self.tot_atoms):
                if j < self.num_atoms[self.chem[eln]]:
                    atomel = self.chem[eln]
                    j += 1
                else:
                    eln += 1
                    atomel = self.chem[eln]
                    j = 1
                acoord = np.genfromtxt(filename, delimiter="", skip_header=8 + i, max_rows=1)
                newatom = atom(id=i + 1, element=atomel, pos=acoord)
                self.basis_atoms.append(newatom)

        mag_header = int(self.tot_atoms + 9)
        self.mag_mom = np.genfromtxt(filename, delimiter="", skip_header=mag_header)
        if self.mag_mom.shape[0] != self.tot_atoms:
            print(
                "WARNING: Insufficient number of magnetic moments. \
					Magnetism ignored."
            )
            self.mag_mom = np.zeros(shape=(self.tot_atoms, 3))

        for i in range(self.tot_atoms):
            self.basis_atoms[i].mag_mom = self.mag_mom[i]

        self.tot_atoms = len(self.basis_atoms)
        print(self.tot_atoms)
        return

    def repeated_length(self, direction=np.zeros(shape=(1, 3))):
        """
        Calculates the repeatable euclidian distance length of the unit cell in
        the Miller Incex direction specified.
        :param direction: Crystallographic direction of axis of interest
        :return: repeatable length
        """
        lengths = np.zeros(shape=(3, 1))

        for i in range(3):
            if direction[i] == 0:
                lengths[i] = 1e6
            else:
                t = 0.5 * self.box[i, i] / float(direction[i])
                lengths[i] = 2 * sqrt(
                    (t * direction[0]) ** 2 + (t * direction[1]) ** 2 + (t * direction[2]) ** 2
                )
        rl = min(lengths)
        return rl

    def hexCij(c11, c12, c13, c33, c44):
        """
        Populates 6x6 elastic tensor for hexagonal crystals
        :param c11: c11 elastic constant
        :param c12: c12 elastic constant
        :param c13: c13 elastic constant
        :param c44: c44 elastic constant
        :return C: 6x6 elastic tensor for Voigt notation calculations
        :return S: 6x6 compliance tensor for Voigt notation calculations
        """
        C = np.zeros(shape=(6, 6))
        C[0, 0] = C[1, 1] = c11
        C[0, 1] = C[1, 0] = c12
        C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = c13
        C[2, 2] = c33
        C[3, 3] = C[4, 4] = c44
        C[5, 5] = 0.5 * (c11 - c12)
        S = np.linalg.inv(C)
        return {"C": C, "S": S}
