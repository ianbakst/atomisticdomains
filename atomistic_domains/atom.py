import numpy as np
from math import *
import time
import os

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
			id: int = 0,
			element: str = 'H',
			position: np.ndarray = np.zeros(shape=(1,3)),
			velocity: np.ndarray = np.zeros(shape=(1,3)),
			magnetic_moment: np.ndarray = np.zeros(shape=(1,3)),
	):
		self.id = id
		self.element = element
		self.position = position
		self.velocity = velocity
		self.magnetic_moment = magnetic_moment
		return


	def set_id(self, new_id: int = 0):
		"""
		Changes or assigns a new atom ID to an atom.
        :param new_id: new atom id
        :return: none
        """
		self.id = new_id
		return self


	def set_element(self, new_element: str = 'H'):
		"""
		Changes chemical element of atom
		:param new_element: new element
		:return: none
		"""
		self.element = new_element
		return

	def set_position(self, new_position: np.ndarray = np.zeros(shape=(1,3))):
		"""
		Changes or updates atomic position to new specified atomic position.
        :param new_pos: new atomic position
        :return: none
        """
		self.position = new_position
		return


	def displace(self, displacement: np.ndarray = np.zeros(shape=(1,3))):
		"""
        Displace an atom by certain values. Default is (0,0,0)
        :param displacement: Displacement vector to be applied to atom
        :return: none
        """
		new_position = self.position + displacement
		return self.set_position(new_position)


	def set_velocity(self, new_velocity: np.ndarray = np.zeros(shape=(1,3))):
		self.velocity = new_velocity
		return


	def accelerate(self, jolt: np.ndarray = np.zeros(shape=(1,3))):
		new_velocity = self.velocity + jolt
		return self.set_velocity(new_velocity)


	def set_magmom(self, new_magnetic_moment: np.ndarray = np.zeros(shape=(1,3))):
		self.magnetic_moment = new_magnetic_moment
		return
