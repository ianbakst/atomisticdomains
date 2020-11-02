import numpy as np
from typing import List
from atomistic_domains.atom import Atom, displace_and_add
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
        atoms: List[Atom] = None,
        x_vector: np.ndarray = np.zeros(shape=(1, 3)),
        y_vector: np.ndarray = np.zeros(shape=(1, 3)),
        z_vector: np.ndarray = np.zeros(shape=(1, 3)),
        boundary_conditions: dict = None,
    ):
        if not atoms:
            self.atoms = []
        else:
            self.atoms = atoms

        self.x_vector = x_vector
        self.y_vector = y_vector
        self.z_vector = z_vector
        if not boundary_conditions:
            self.boundary_conditions = {
                "x_boundary": "p",
                "y_boundary": "p",
                "z_boundary": "p",
            }
        else:
            self.boundary_conditions = boundary_conditions
        return


def create_from_lattice(
    lattice: Lattice,
    nx: int = 1,
    ny: int = 1,
    nz: int = 1,
    boundary_conditions: dict = None,
):
    atoms = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                displacement = i * lattice.x_vector + j * lattice.y_vector + k * lattice.z_vector
                atoms.extend(
                    [displace_and_add(atom, displacement) for atom in lattice.basis_atoms]
                )
    return create(
        atoms=atoms,
        x_vector=nx * lattice.x_vector,
        y_vector=ny * lattice.y_vector,
        z_vector=nz * lattice.z_vector,
        boundary_conditions=boundary_conditions,
    )


def create(
    atoms: List[Atom] = None,
    x_vector: np.ndarray = np.zeros(shape=(1, 3)),
    y_vector: np.ndarray = np.zeros(shape=(1, 3)),
    z_vector: np.ndarray = np.zeros(shape=(1, 3)),
    boundary_conditions: dict = None,
):

    return Domain(
        atoms=atoms,
        x_vector=x_vector,
        y_vector=y_vector,
        z_vector=z_vector,
        boundary_conditions=boundary_conditions,
    )
