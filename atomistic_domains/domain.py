import numpy as np
from typing import List
from atomistic_domains.atom import Atom, displace_and_add
from atomistic_domains.lattice import Lattice
import pandas as pd


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
        x_vector: np.ndarray = np.zeros(shape=(1, 3), dtype=float),
        y_vector: np.ndarray = np.zeros(shape=(1, 3), dtype=float),
        z_vector: np.ndarray = np.zeros(shape=(1, 3), dtype=float),
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

    def write_poscar(
            self,
            file_name: str,
            comment: str = 'Structure',
            scale: float = 1.0,
            selective_dynamics: bool = False,
            direct: bool = True,
    ):
        r = [{'element': atom.element, 'position': atom.position, 'velocity': atom.velocity} for atom in self.atoms]
        df = pd.DataFrame(r)
        with open(file_name, 'w') as f:
            f.write(f"{comment}\n")
            f.write(f"{scale}\n")
            for vec in (self.x_vector, self.y_vector, self.z_vector):
                f.write(f"{' '.join([str(float(i)) for i in vec])}\n")

            element_list = df.element.unique()

            for el in element_list:
                f.write(f" {el}")
            f.write("\n")

            for el in element_list:
                f.write(f" {len(df[df.element == el])}")
            f.write("\n")

            if selective_dynamics:
                f.write('Selective Dynamics\n')

            f.write('Direct\n') if direct else f.write('Cartesian\n')

            for el in element_list:
                for _, row in df[df.element == el].iterrows():
                    pos = fractionalize(row.position) if direct else row.position
                    f.write(f"{' '.join([str(float(i)) for i in pos])}\n")
                    if selective_dynamics:
                        f.write("")
                    f.write("\n")
        return


    def fractionalize(self, position: np.ndarray) -> np.ndarray:
        


def create_from_lattice(
    lattice: Lattice,
    nx: int = 1,
    ny: int = 1,
    nz: int = 1,
    boundary_conditions: dict = None,
) -> Domain:
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
) -> Domain:

    return Domain(
        atoms=atoms,
        x_vector=x_vector,
        y_vector=y_vector,
        z_vector=z_vector,
        boundary_conditions=boundary_conditions,
    )
