import numpy as np
from math import *
import time
import os
def mbplane(A,c=sqrt(8/3)):
	"""
	Converts Miller-Bravais plane indices to cartesian normal vectors
	:param A:
	:param c:
	:return:
	"""
	#for Sun Mg Potential: c=1.6281689374348
	if A[0]==0 and A[1]==0 and A[2]==0:
		B=np.array([0,0,1])

	else:
		if A[0]==0:
			p1=np.array([-0.5/A[1],0.5*sqrt(3)/A[1],0])
			p2=np.array([-0.5/A[2],-0.5*sqrt(3)/A[2],0])
		elif A[1]==0:
			p1=np.array([-0.5/A[2],-0.5*sqrt(3)/A[2],0])
			p2=np.array([1/A[0],0,0])
		else:
			p1=np.array([1/A[0],0,0])
			p2=np.array([-0.5/A[1],0.5*sqrt(3)/A[1],0])
		if A[3]==0:
			z=p1+np.array([0,0,1])
		else:
			z=np.array([0,0,c/A[3]])

		P1=p1-z
		P2=p2-z
		B=np.cross(P1,P2)
	return (B)

def mbvector(A,c=sqrt(8/3)):
	"""
	Converts non-cubic crystal indices to cartesian-coordinate vectors
	whilst preserving vector length.
	:param A: list of row vectors in crystal-index form
	:param c: constant(s) or basis vectors which describe the crystal unit cell.
	Default is HCP
	:return: list of row vectors in cartesian-corrdinates
	"""
	la=len(A)
	sa=A.size
	if la==sa:
		B=np.array([0.0,0.0,0.0])
		a1=A[0]*np.array([1.0,0.0])
		a2=A[1]*np.array([-0.5,0.5*sqrt(3)])
		a3=A[2]*np.array([-0.5,-0.5*sqrt(3)])
		B[0]=a1[0]+a2[0]+a3[0]
		B[1]=a1[1]+a2[1]+a3[1]
		B[2]=c*A[3]
	else:
		sa=A.shape
		B=np.zeros(shape=(sa[0],3))
		for i in range(sa[0]):
			B[i,0]=a1[0]+a2[0]+a3[0]
			B[i,1]=a1[1]+a2[1]+a3[1]
			B[i,2]=c*A[i,3]
	return (B)

def mvector(B,c):
	"""Converts cartesian normal vectors into Miller-Bravais plane indices"""
	#for Sun Mg Potential: c=1.6281689374348
	A=np.zeros(shape=4)
	A[0]=(2/3)*B[0]
	A[1]=0.5*((2/sqrt(3))*B[1]-A[0])
	A[2]=-A[0]-A[1]
	A[3]=B[2]/c
	return (A)

def mplane(B,c):
	"""Converts cartesian vectors into Miller-Bravais vector indices"""
	#for Sun Mg Potential: c=1.6281689374348
	A=np.zeros(shape=4)
	s=A.shape
	phi=0
	l=0
	b=0
	for i in range(s[0]):
		A[0]=(2/3)*B[0]
		A[1]=0.5*((2/sqrt(3))*B[1]-A[0])
		A[2]=-A[0]-A[1]
		phi=atan(B[1]/B[0])
		b=0.5*sqrt(3)/sin(2*pi/3-phi)
		l=sqrt(B[0]**2+B[1]**2)
		A[3]=c*B[2]/(l*b)
	return (A)

def rotation(x1,z1,x2,z2):
	"""
	Rotates coordinate system from x1,z1 to x2,z2
	:param x1: direction of initial x-axis
	:param z1: direction of initial z-axis
	:param x2: direction of final x-axis
	:param z2: direction of final z-axis
	:return: rotation matrix between the initial and final coordinate systems
	"""
	e1=np.zeros(shape=(3,3))
	e2=np.zeros(shape=(3,3))
	e1[0,:]=x1/np.linalg.norm(x1)
	e1[2,:]=z1/np.linalg.norm(z1)
	e1[1,:]=np.cross(e1[2,:],e1[0,:])
	e2[0,:]=x2/np.linalg.norm(x2)
	e2[2,:]=z2/np.linalg.norm(z2)
	e2[1,:]=np.cross(e2[2,:],e2[0,:])
	R=np.zeros(shape=(3,3))
	for i in range(3):
		for j in range(3):
			R[i,j]=np.dot(e1[i,:],e2[j,:])
	R=np.transpose(R)
	return (R)

def strain2disp(e,A):
	"""
	Converts strain state to displacements
	With the convention that ux=f(x,y), uy=f(y,z), uz=f(z,x)
	:param e: 3x3 symmetric strain tensor
	:param A: list of atomic positions (rows)
	:return: displacements of atoms given uniform strain state.
	"""
	l=A.shape
	u=np.zeros(shape=l)
	for i in range(l[0]):
		u[i,0]=e[0,0]*A[i,0]+e[0,1]*A[i,1]
		u[i,1]=e[1,1]*A[i,1]+e[1,2]*A[i,2]
		u[i,2]=e[2,2]*A[i,2]+e[0,2]*A[i,0]
	return (u)

def stress2strain(sigma,S):
	"""
	Converts stress state to strain state
	:param sigma: 6x1 Voigt vector of stresses
	:param S: 6x6 compliance tensor (inverse of C)
	:return: 6x1 Voigt vector of strains
	"""
	s6=voigt(sigma)
	e6=np.dot(S,s6)
	e6[3]=0.5*e6[3]
	e6[4]=0.5*e6[4]
	e6[5]=0.5*e6[5]
	e=unvoigt(e6)
	return (e)

def unvoigt(A):
	"""
	Converts from 6x1 to 3x3
	:param A: 6x1 Voigt vector (strain or stress)
	:return: 3x3 symmetric tensor (strain or stress)
	"""
	a=np.zeros(shape=(3,3))
	a[0,0]=A[0]
	a[0,1]=A[5]
	a[0,2]=A[4]
	a[1,0]=A[5]
	a[1,1]=A[1]
	a[1,2]=A[3]
	a[2,0]=A[4]
	a[2,1]=A[3]
	a[2,2]=A[2]
	return (a)

def vect_contract(m,c,n):
	"""Performs vector contraction in the form of m C n"""
	a=np.tensordot(m,c,(0,0))
	mn=np.tensordot(a,n,(2,0))
	return (mn)

def voigt(a):
	"""
	Converts from 3x3 to 6x1
	:param a: 3x3 symmetric tensor (strain or stress)
	:return: 6x1 Voigt vector (strain or stress)
	"""
	A=np.zeros(shape=(6,1))
	A[0]=a[0,0]
	A[1]=a[1,1]
	A[2]=a[2,2]
	A[3]=a[1,2]
	A[4]=a[0,2]
	A[5]=a[0,1]
	return (A)
