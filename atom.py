import numpy as np
from math import *
import time
import os
from base import *

class atom:
	"""
    The class "atom" defines atom objects which have the attributes:
    :param id: atom id number
    :param element: chemical element
    :param pos: x,y,z-coordinates of atom position
    :param mag_mom: magnetic moment of atom (only used in class: vasp)
    :param dyn: selective dynamics to be applied to this atom (only used in class:vasp)
    """
	def __init__(self,id=0,element='H',pos=np.zeros(shape=(1,3)),mag_mom=np.zeros(shape=(1,3)),dyn=['T','T','T']):
		self.id=id
		self.element=element
		self.pos=pos
		self.rpos=rpos
		self.vel=vel
		self.mag=mag
		self.dyn=dyn
		return


	def set_dyn(self,new_dyn=['T','T','T']):
		"""
		Changes dynamics associated with atom (for vasp structures)
        :param new_dyn: New dynamics specifications
        :return: none
        """
		self.dyn=new_dyn
		return

	def set_element(self,new_element):
		"""
		Changes chemical element of atom
        :param new_element: new element
        :return: none
        """
		self.element=new_element
		return

	def set_id(self,new_id=0):
		"""
		Changes or assigns a new atom ID to an atom.
        :param new_id: new atom id
        :return: none
        """
		self.id=new_id
		return

	def set_pos(self,new_pos=np.zeros(shape=(1,3))):
		"""
		Changes or updates atomic position to new specified atomic position.
        :param new_pos: new atomic position
        :return: none
        """
		self.pos=new_pos
		return
	def set_rpos(self,new_pos=np.zeros(shape=(1,3))):
		"""
		Changes or updates atomic position to new specified atomic position.
        :param new_pos: new atomic position
        :return: none
        """
		self.rpos=new_pos
		return

	def displace(self,displacement=np.zeros(shape=(1,3))):
		"""
        Displace an atom by certain values. Default is (0,0,0)
        :param displacement: Displacement vector to be applied to atom
        :return: none
        """
		self.pos += displacement
		return

	def scale(self,scale_factor=np.array([1.0])):
		"""
        Scales coordinates of an atom by either constant value or mapping box
        :param scale_factor: either scalar or 3x3 mapping matrix
        :return: none
        """
		self.pos=np.dot(self.pos,scale_factor)
		return
