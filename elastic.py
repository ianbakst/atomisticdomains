import numpy as np
from math import *
import time
import os
from .base import *

class elastic:
    """
    Elastic tensor object which will be used for calculations/methods requiring
    the elastic constants of a material. Methods of this object will manipulate
    all elements of the elastic tensor.
    """
    def __init__(self, cijkl=np.zeros(shape=(3,3,3,3)),Cij=np.zeros(shape=(6,6))):

        return

    def Cijkl(C):
    	"""
    	Populates 4D elastic tensor from 6x6 elastic tensor
    	:param C: 6x6 elastic constants tensor
    	:return: 3x3x3x3 elastic Cijkl tensor
    	"""
    	c=np.zeros(shape=(3,3,3,3))
    	CC=np.zeros(shape=(9,9))
    	CC[0:6,0:6] = C[0:6,0:6]
    	CC[6:9,6:9] = C[3:6,3:6]
    	CC[0:6,6:9] = C[0:6,3:6]
    	CC[6:9,0:6] = C[3:6,0:6]

    	c[0,0,0,0] = CC[0,0]
    	c[0,0,1,1] = CC[0,1]
    	c[0,0,2,2] = CC[0,2]
    	c[0,0,1,2] = CC[0,3]
    	c[0,0,2,0] = CC[0,4]
    	c[0,0,0,1] = CC[0,5]
    	c[0,0,2,1] = CC[0,6]
    	c[0,0,0,2] = CC[0,7]
    	c[0,0,1,0] = CC[0,8]

    	c[1,1,0,0] = CC[1,0]
    	c[1,1,1,1] = CC[1,1]
    	c[1,1,2,2] = CC[1,2]
    	c[1,1,1,2] = CC[1,3]
    	c[1,1,2,0] = CC[1,4]
    	c[1,1,0,1] = CC[1,5]
    	c[1,1,2,1] = CC[1,6]
    	c[1,1,0,2] = CC[1,7]
    	c[1,1,1,0] = CC[1,8]

    	c[2,2,0,0] = CC[2,0]
    	c[2,2,1,1] = CC[2,1]
    	c[2,2,2,2] = CC[2,2]
    	c[2,2,1,2] = CC[2,3]
    	c[2,2,2,0] = CC[2,4]
    	c[2,2,0,1] = CC[2,5]
    	c[2,2,2,1] = CC[2,6]
    	c[2,2,0,2] = CC[2,7]
    	c[2,2,1,0] = CC[2,8]

    	c[1,2,0,0] = CC[3,0]
    	c[1,2,1,1] = CC[3,1]
    	c[1,2,2,2] = CC[3,2]
    	c[1,2,1,2] = CC[3,3]
    	c[1,2,2,0] = CC[3,4]
    	c[1,2,0,1] = CC[3,5]
    	c[1,2,2,1] = CC[3,6]
    	c[1,2,0,2] = CC[3,7]
    	c[1,2,1,0] = CC[3,8]

    	c[2,0,0,0] = CC[4,0]
    	c[2,0,1,1] = CC[4,1]
    	c[2,0,2,2] = CC[4,2]
    	c[2,0,1,2] = CC[4,3]
    	c[2,0,2,0] = CC[4,4]
    	c[2,0,0,1] = CC[4,5]
    	c[2,0,2,1] = CC[4,6]
    	c[2,0,0,2] = CC[4,7]
    	c[2,0,1,0] = CC[4,8]

    	c[0,1,0,0] = CC[5,0]
    	c[0,1,1,1] = CC[5,1]
    	c[0,1,2,2] = CC[5,2]
    	c[0,1,1,2] = CC[5,3]
    	c[0,1,2,0] = CC[5,4]
    	c[0,1,0,1] = CC[5,5]
    	c[0,1,2,1] = CC[5,6]
    	c[0,1,0,2] = CC[5,7]
    	c[0,1,1,0] = CC[5,8]

    	c[2,1,0,0] = CC[6,0]
    	c[2,1,1,1] = CC[6,1]
    	c[2,1,2,2] = CC[6,2]
    	c[2,1,1,2] = CC[6,3]
    	c[2,1,2,0] = CC[6,4]
    	c[2,1,0,1] = CC[6,5]
    	c[2,1,2,1] = CC[6,6]
    	c[2,1,0,2] = CC[6,7]
    	c[2,1,1,0] = CC[6,8]

    	c[0,2,0,0] = CC[7,0]
    	c[0,2,1,1] = CC[7,1]
    	c[0,2,2,2] = CC[7,2]
    	c[0,2,1,2] = CC[7,3]
    	c[0,2,2,0] = CC[7,4]
    	c[0,2,0,1] = CC[7,5]
    	c[0,2,2,1] = CC[7,6]
    	c[0,2,0,2] = CC[7,7]
    	c[0,2,1,0] = CC[7,8]

    	c[1,0,0,0] = CC[8,0]
    	c[1,0,1,1] = CC[8,1]
    	c[1,0,2,2] = CC[8,2]
    	c[1,0,1,2] = CC[8,3]
    	c[1,0,2,0] = CC[8,4]
    	c[1,0,0,1] = CC[8,5]
    	c[1,0,2,1] = CC[8,6]
    	c[1,0,0,2] = CC[8,7]
    	c[1,0,1,0] = CC[8,8]
    	return (c)

    def hexCij(c11,c12,c13,c33,c44):
		"""
		Populates 6x6 elastic tensor for hexagonal crystals
		:param c11: c11 elastic constant
		:param c12: c12 elastic constant
		:param c13: c13 elastic constant
		:param c44: c44 elastic constant
		:return C: 6x6 elastic tensor for Voigt notation calculations
		:return S: 6x6 compliance tensor for Voigt notation calculations
		"""
		C=np.zeros(shape=(6,6))
		C[0,0]=C[1,1]=c11
		C[0,1]=C[1,0]=c12
		C[0,2]=C[2,0]=C[1,2]=C[2,1]=c13
		C[2,2]=c33
		C[3,3]=C[4,4]=c44
		C[5,5]=0.5*(c11-c12)
		S=np.linalg.inv(C)
		return {'C':C,'S':S}

    def cubeCij(c11,c13,c44):
    	"""
    	Populates 6x6 elastic tensor for cubic crystals
    	:param c11:
    	:param c13:
    	:param c44:
    	:return: 6x6 elastic stiffness tensor
    	"""
    	C=np.zeros(shape=(6,6))
    	C[0,0]=C[1,1]=C[2,2]=c11
    	C[0,1]=C[0,2]=C[1,0]=C[1,2]=C[2,0]=C[2,1]=c13
    	C[3,3]=C[4,4]=C[5,5]=c44
    	S=np.linalg.inv(C)
    	return {'C':C,'S':S}
