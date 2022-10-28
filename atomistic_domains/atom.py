import numpy as np


class Atom:
    """
    The class "Atom" defines atom objects which have the attributes:
    :param id: atom id number
    :param element: chemical element
    :param pos: x,y,z-coordinates of atom position
    :param mag_mom: magnetic moment of atom (only used in class: vasp)
    :param dyn: selective dynamics to be applied to this atom (only used in class:vasp)
    """

    def __init__(
        self,
        atom_id: int = 0,
        element: str = "Nu",
        position: np.ndarray = np.zeros(shape=(1, 3)),
        velocity: np.ndarray = np.zeros(shape=(1, 3)),
        magnetic_moment: np.ndarray = np.zeros(shape=(1, 3)),
    ):
        self._id = atom_id
        self._element = element
        self._position = position
        self._velocity = velocity
        self._magnetic_moment = magnetic_moment

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, new_id: int = 0):
        """
        Changes or assigns a new atom ID to an atom.
        :param new_id: new atom id
        :return: none
        """
        self._id = new_id

    @property
    def element(self):
        return self._element

    @element.setter
    def element(self, new_element: str = "Nu"):
        """
        Changes chemical element of atom
        :param new_element: new element
        :return: none
        """
        self._element = new_element

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, new_position: np.ndarray = np.zeros(shape=(1, 3))):
        """
                Changes or updates atomic position to new specified atomic position.
        :param new_pos: new atomic position
        :return: none
        """
        self._position = new_position

    @property
    def velocity(self):
        return self._velocity

    @velocity.setter
    def velocity(self, new_velocity: np.ndarray = np.zeros(shape=(1, 3))):
        self._velocity = new_velocity

    @property
    def magnetic_moment(self):
        return self._magnetic_moment

    @magnetic_moment.setter
    def magnetic_moment(self, new_magnetic_moment: np.ndarray = np.zeros(shape=(1, 3))):
        self._magnetic_moment = new_magnetic_moment


def create(
    atom_id: int = 0,
    element: str = "Nu",
    position: np.ndarray = np.zeros(shape=(1, 3)),
    velocity: np.ndarray = np.zeros(shape=(1, 3)),
    magnetic_moment: np.ndarray = np.zeros(shape=(1, 3)),
) -> Atom:
    return Atom(
        atom_id=atom_id,
        element=element,
        position=position,
        velocity=velocity,
        magnetic_moment=magnetic_moment,
    )


def displace_and_add(
    atom: Atom,
    displacement: np.ndarray = np.zeros(shape=(1, 3)),
) -> Atom:
    new_atom = create(
        atom.id,
        atom.element,
        atom.position,
        atom.velocity,
        atom.magnetic_moment,
    )
    new_atom.position += displacement
    return new_atom
