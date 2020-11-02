import numpy as np
from math import *
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


def create(
        x_vector: np.ndarray = np.zeros(shape=(1, 3)),
        y_vector: np.ndarray = np.zeros(shape=(1, 3)),
        z_vector: np.ndarray = np.zeros(shape=(1, 3)),
        basis_atoms: List[Atom] = None,
) -> Lattice:
    return Lattice(
        x_vector,
        y_vector,
        z_vector,
        basis_atoms,
    )

