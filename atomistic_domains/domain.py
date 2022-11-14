import numpy as np
from typing import List, Optional
from atomistic_domains.atom import Atom, displace_and_add
from atomistic_domains.boundary import BoundaryConditions
from atomistic_domains.lattice import Lattice


class Domain:
    """
    The class "domain" defines methods which manipulate atomistic domain structures.
    Domain structures have the attributes:
    :param box: High and low coordinates of simulation box boundaries. (3x2 array)
    :param chem: list of chemical elements
    :param num_el: number of element
    :param num_atoms: number of atoms of each element (dictionary type)
    :param tot_atoms: Total number of atoms in domain.
    :param atoms: list of atoms in the simulation domain (see class "atom")
    :param x: crystallographic direction of the positive x-axis
    :param y: crystallographic direction of the positive y-axis
    :param z: crystallographic direction of the positive z-axis
    :param p: periodicity of simulation domain
    """

    def __init__(
        self,
        *,
        atoms: Optional[List[Atom]] = None,
        x1_vector: np.ndarray = np.zeros(shape=(1, 3)),
        x2_vector: np.ndarray = np.zeros(shape=(1, 3)),
        x3_vector: np.ndarray = np.zeros(shape=(1, 3)),
        boundary_conditions: Optional[BoundaryConditions] = None,
    ):
        if atoms is None:
            self.atoms = []
        else:
            self.atoms = atoms
        self._x1_vector = x1_vector
        self._x2_vector = x2_vector
        self._x3_vector = x3_vector
        if boundary_conditions is None:
            self.boundary_conditions = BoundaryConditions(1, 1, 1)
        else:
            self.boundary_conditions = boundary_conditions

    @property
    def x1_vector(self):
        return self._x1_vector

    @x1_vector.setter
    def x1_vector(self, new_vector):
        self._x1_vector = new_vector

    @property
    def x2_vector(self):
        return self._x1_vector

    @x2_vector.setter
    def x2_vector(self, new_vector):
        self._x2_vector = new_vector

    @property
    def x3_vector(self):
        return self._x3_vector

    @x3_vector.setter
    def x3_vector(self, new_vector):
        self._x3_vector = new_vector


def create_from_lattice(
    lattice: Lattice,
    nx: int = 1,
    ny: int = 1,
    nz: int = 1,
    boundary_conditions: Optional[BoundaryConditions] = None,
):
    atoms = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                displacement = (
                    i * lattice.x1_vector + j * lattice.x2_vector + k * lattice.x3_vector
                )
                atoms.extend(
                    [displace_and_add(atom, displacement) for atom in lattice.basis_atoms]
                )
    return create(
        atoms=atoms,
        x1_vector=nx * lattice.x1_vector,
        x2_vector=ny * lattice.x2_vector,
        x3_vector=nz * lattice.x3_vector,
        boundary_conditions=boundary_conditions,
    )


def create(
    atoms: List[Atom] = None,
    x1_vector: np.ndarray = np.zeros(shape=(1, 3)),
    x2_vector: np.ndarray = np.zeros(shape=(1, 3)),
    x3_vector: np.ndarray = np.zeros(shape=(1, 3)),
    boundary_conditions: Optional[BoundaryConditions] = None,
):
    return Domain(
        atoms=atoms,
        x1_vector=x1_vector,
        x2_vector=x2_vector,
        x3_vector=x3_vector,
        boundary_conditions=boundary_conditions,
    )
