from copy import deepcopy
from math import acos, sqrt

import numpy as np
from typing import List, Optional, Union
from atomistic_domains import Atom


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
        vector_1: np.ndarray = np.zeros(shape=(1, 3)),
        vector_2: np.ndarray = np.zeros(shape=(1, 3)),
        vector_3: np.ndarray = np.zeros(shape=(1, 3)),
        basis_atoms: List[Atom] = None,
    ):
        self._x1_vector = vector_1
        self._x2_vector = vector_2
        self._x3_vector = vector_3
        if not basis_atoms:
            basis_atoms = []
        self.basis_atoms = basis_atoms

    @property
    def x1_vector(self):
        return self._x1_vector

    @x1_vector.setter
    def x1_vector(self, new_vector: np.ndarray = np.zeros(shape=(1, 3))):
        self._x1_vector = new_vector

    @property
    def x2_vector(self):
        return self._x2_vector

    @x2_vector.setter
    def x2_vector(self, new_vector: np.ndarray = np.zeros(shape=(1, 3))):
        self._x2_vector = new_vector

    @property
    def x3_vector(self):
        return self._x3_vector

    @x3_vector.setter
    def x3_vector(self, new_vector: np.ndarray = np.zeros(shape=(1, 3))):
        self._x3_vector = new_vector

    @property
    def x1_length(self):
        return np.linalg.norm(self._x1_vector)

    @property
    def x2_length(self):
        return np.linalg.norm(self._x2_vector)

    @property
    def x3_length(self):
        return np.linalg.norm(self._x3_vector)

    @property
    def alpha(self):
        return acos(np.dot(self.x1_vector, self.x2_vector) / self.x1_length / self.x2_length)

    @property
    def beta(self):
        return acos(np.dot(self.x1_vector, self.x3_vector) / self.x1_length / self.x3_length)

    @property
    def gamma(self):
        return acos(np.dot(self.x3_vector, self.x2_vector) / self.x3_length / self.x2_length)

    @property
    def volume(self):
        abc = self.x1_length * self.x2_length * self.x3_length
        cosa = np.dot(self.x1_vector, self.x2_vector) / self.x1_length / self.x2_length
        cosb = np.dot(self.x1_vector, self.x3_vector) / self.x1_length / self.x3_length
        cosg = np.dot(self.x3_vector, self.x2_vector) / self.x3_length / self.x2_length
        return abc * sqrt(1 - cosa**2 - cosb**2 - cosg**2 + 2 * cosa * cosb * cosg)


def create(
    *,
    x1_vector: Union[List[Union[float, int]], np.ndarray] = np.zeros(shape=(1, 3)),
    x2_vector: Union[List[Union[float, int]], np.ndarray] = np.zeros(shape=(1, 3)),
    x3_vector: Union[List[Union[float, int]], np.ndarray] = np.zeros(shape=(1, 3)),
    basis_atoms: List[Atom] = None,
) -> Lattice:

    return Lattice(
        x1_vector if isinstance(x1_vector, np.ndarray) else np.array(x1_vector),
        x2_vector if isinstance(x2_vector, np.ndarray) else np.array(x2_vector),
        x3_vector if isinstance(x3_vector, np.ndarray) else np.array(x3_vector),
        basis_atoms,
    )


class BCC(Lattice):
    def __init__(
        self, side_length: float = 1, atom1: Optional[Atom] = None, atom2: Optional[Atom] = None
    ):
        x1 = side_length * np.array([1, 0, 0])
        x2 = side_length * np.array([0, 1, 0])
        x3 = side_length * np.array([0, 0, 1])
        if atom1 is None:
            atom1 = Atom()
        atom1.position = np.array([0, 0, 0])
        if atom2 is None:
            atom2 = deepcopy(atom1)
        atom2.position = np.array([0.5 * side_length] * 3)
        super().__init__(x1, x2, x3, [atom1, atom2])
