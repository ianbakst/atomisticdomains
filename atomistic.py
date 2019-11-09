#!/usr/bin/python
"""
Notes on the CW Group python library:
	When specifying dimensions (to make periodic or center points around, etc.)
		the convention is: 1 = x, 2 = y, 4 = z. Combinations of those dimensions are the
		sums of their component dimensions (e.g. x and y = 3).
	Unless otherwise specified, elastic constants should be in Pa. Atom positions should
		be in Angstroms
"""
import numpy as np
from math import *
import time
import os

class vasp:
	"""
	The class "vasp" defines methods which manipulate VASP structures.
	VASP structures have the attributes:
	:param scale: Scale factor of atom positions
	:param box: Box dimensions and tilts
	:param a: x-side box length
	:param b: y-side box length
	:param c: z-side box length
	:param chem: List of elements in this simulation box
	:param num_atoms: Lists of number of atoms for each element
	:param tot_atoms: Total number of atoms in simulation box
	:param coord: List of coordinates
	:param dyn: List of selective dynamics attributes.
	:param mag_mom: List of magnetic moments
	"""

	def __init__(self,scale=1.0,box=np.zeros(shape=(3,3)),chem=np.zeros(shape=(1,1),dtype=str),num_atoms=np.zeros(shape=(1,1)),coord=np.zeros(shape=(1,1)),mag_mom=np.zeros(shape=(1,1)),dyn=np.zeros(shape=(1,1))):

		self.scale=scale
		self.dyn=dyn

		self.box=box
		self.a=box[0,0]
		self.b=box[1,1]
		self.c=box[2,2]
		self.num_atoms = num_atoms
		self.tot_atoms = sum(num_atoms)

		if chem.shape[1]!=num_atoms.shape[1]:
			print 'WARNING: Inconsistent number of atoms and number of element types.'
			self.chem = np.zeros(shape=(1, num_atoms.shape[1]), dtype=str)
			for i in range(num_atoms.shape[1]):
				aname = 'atom' + str(i)
				self.chem[0, i] = aname
		else:
			self.chem=chem

		if coord.shape[1]==3 and coord.shape[0]==self.tot_atoms:
			self.coord=coord
		elif coord.shape[1]==1 and coord.shape[0]==1:
			print 'WARNING: Initializing blank VASP object.'
			self.coord=coord
		else:
			print 'FATAL ERROR: Insufficient number of coordinates.'

		if mag_mom.shape[0]==self.tot_atoms:
			self.mag_mom=mag_mom
		else:
			print 'WARNING: Insufficient number of magnetic moments. Magnetism ignored.'
			self.mag_mom=np.zeros(shape=(self.tot_atoms,3))
		return

	def readposcar(self,filename):
		"""
		Populates a VASP object directly from a POSCAR
		:param filename: name of POSCAR
		:return: none
		"""
		self.scale=np.genfromtxt(filename,skip_header=1,max_rows=1)
		self.box=np.genfromtxt(filename,delimiter='',skip_header=2,max_rows=3)
		self.a=self.box[0,0]
		self.b=self.box[1,1]
		self.c=self.box[2,2]
		self.chem=np.genfromtxt(filename,delimiter='',skip_header=5,max_rows=1,dtype=str)
		self.num_atoms=np.genfromtxt(filename,delimiter='',skip_header=6,max_rows=1,dtype=int)
		self.tot_atoms=sum(self.num_atoms)
		ptype = np.genfromtxt(filename, delimiter='>', skip_header=7, max_rows=1, dtype=str)
		if ptype == 'Selective Dynamics':
			self.coord=np.genfromtxt(filename,delimiter='',skip_header=9,max_rows=self.tot_atoms)[:,0:3]
			self.dyn=np.genfromtxt(filename,delimiter='',skip_header=9,max_rows=self.tot_atoms,dtype=str)[:,3:]
		else:
			self.coord=np.genfromtxt(filename,delimiter='',skip_header=8,max_rows=self.tot_atoms)
		mag_header=int(self.tot_atoms+9)
		self.mag_mom=np.genfromtxt(filename,delimiter='',skip_header=mag_header)
		if self.mag_mom.shape[0]!=self.tot_atoms:
			print 'WARNING: Insufficient number of magnetic moments. Magnetism ignored.'
			self.mag_mom=np.zeros(shape=(self.tot_atoms,3))
		return

	def shift(self,cut=0.0,dx=0.0,dy=0.0,dz=0.0,plane='z'):
		"""
		Shifts atoms above a specified plane by specified displacements. Useful for generating GSF curves.
		:param cut: Height of plane above which atoms will be displaced
		:param dx: Distance atoms will be moved in x-direction
		:param dy: Distance atoms will be moved in y-direction
		:param dz: Distance atoms will be moved in z-direction
		:param plane: Plane above which atoms will be displaced
		:return: New vasp object with updated coordinates.
		"""
		new=self
		if plane=='x' or plane=='X':
			for i in range(self.tot_atoms):
				if self.coord[i,0]>=cut:
					new.coord[i,1]+=dy
					new.coord[i,2]+=dz
				else:
					pass
		elif plane=='y' or plane=='Y':
			for i in range(self.tot_atoms):
				if self.coord[i, 1] >= cut:
					new.coord[i, 0] += dx
					new.coord[i, 2] += dz
				else:
					pass
		elif plane=='z' or plane=='Z':
			for i in range(self.tot_atoms):
				if self.coord[i, 2] >= cut:
					new.coord[i, 0] += dx
					new.coord[i, 1] += dy
				else:
					pass
		return (new)

	def writeposcar(self,filename,comment='Structure',ptype='Direct'):
		"""
		Writes vasp object to POSCAR
		:param filename: name of POSCAR to be written
		:param comment: Opening line comment
		:param ptype: POSCAR type, based on dynamics
		:return: NONE
		"""
		f=open(filename, 'w')
		f.write(comment+'\n')
		f.write('%8.12f \n' %(self.scale))
		f.write('%8.12f    %8.12f    %8.12f\n' %(self.box[0,0],self.box[0,1],self.box[0,2]))
		f.write('%8.12f    %8.12f    %8.12f\n' % (self.box[1,0],self.box[1,1],self.box[1,2]))
		f.write('%8.12f    %8.12f    %8.12f\n' % (self.box[2,0],self.box[2,1],self.box[2,2]))
		num_el=self.num_atoms.shape[0]
		for i in range(num_el):
			f.write('%s   ' %(self.chem[i])),
		f.write('\n')
		for i in range(num_el):
			f.write('%d   ' % (self.num_atoms[i])),
		f.write('\n')
		f.write('%s\n' % (ptype))
		if ptype=='Direct':
			for i in range(self.tot_atoms):
				f.write('%8.12f %8.12f %8.12f\n' % (self.coord[i,0], self.coord[i,1], self.coord[i,2]))
		elif ptype=='Selective Dynamics':
			for i in range(self.tot_atoms):
				f.write('Cartesian\n')
				f.write('%8.12f %8.12f %8.12f %s %s %s\n' % (self.coord[i, 0],self.coord[i, 1],self.coord[i, 2],self.dyn[i,0],self.dyn[i,1],self.dyn[i,2]))
		else:
			pass
		f.write('\n')
		for i in range(self.tot_atoms):
			f.write('%8.12f %8.12f %8.12f\n' %(self.mag_mom[i,0],self.mag_mom[i,1],self.mag_mom[i,2]))
		return ('POSCAR written to file: '+filename)


"""
class lmp:
	#""
	The class "lmp" defines methods which manipulate LAMMPS structures.
	LAMMPS structures have the attributes:
	:param num_atoms: number of atoms in the simulation box
	:param box_bounds: 3x2 list, (xlo xhi, ylo yhi, zlo zhi)
	:param box_lengths:
	:param id: list of atom IDs
	#""
	def __init__(self,num_atoms=np.zeros(shape=(1,1)),box=np.zeros(shape=(3,2)),id=np.zeros(shape=(1,1)),):

	def read_dump(self,filename):

	def read_data(self,filename):

	def write_dump(self,filename):

	def write_data(self,filename):

	def write_neb_end(self,filename):
"""
def add_atoms(basis_atoms,unit_cell,sim_box,o=np.array([0.0,0.0,0.0])):
	"""
	Constructs list of atom positions based on basis atoms and box dimensions
	:param basis_atoms: List of basis atoms (row vectors)
	:param unit_cell: List of unit cell box-side lengths
	:param sim_box: List of side lengths for simulation box in which atoms will be added
	:param o: location of desired origin, considering the box is initially centered at (0,0,0).
	"""
	lb=basis_atoms.shape

	nx=ceil(sim_box[0]/unit_cell[0])
	ny=ceil(sim_box[1]/unit_cell[1])
	nz=ceil(sim_box[2]/unit_cell[2])

	print nx
	print ny
	print nz

	xstart=-int(nx/2)
	ystart=-int(ny/2)
	zstart=-int(nz/2)
	if (nx%2)==1:
		xend=int(nx/2+1)
	else:
		xend=int(nx/2)
	if (ny%2)==1:
		yend=int(ny/2+1)
	else:
		yend=int(ny/2)
	if (nz%2)==1:
		zend=int(nz/2+1)
	else:
		zend=int(nz/2)

	n=nx*ny*nz*lb[0]
	A=np.zeros(shape=(n,3))
	cnt=0
	for i in range(xstart,xend):
		for j in range(ystart,yend):
			for k in range(zstart,zend):
				for l in range(lb[0]):
					A[cnt,0]=basis_atoms[l,0]+i*unit_cell[0]
					A[cnt,1]=basis_atoms[l,1]+j*unit_cell[1]
					A[cnt,2]=basis_atoms[l,2]+k*unit_cell[2]
					cnt+=1
	for i in range(int(n)):
		for j in range(3):
			A[i,j]+=o[j]
	return (A)

def bjs(l,c):
	"""
	Computes Bjs matrix for a defect in a given anisotropic medium
	:param l: line direction of defect
	:param c: 4th order elastic tensor
	:return: Bjs matrix
	"""
	if len(l)==4:
		l=mbvector(l)
	elif len(l)==3:
		pass
	else:
		return (0)
	v=np.array([1,pi,e])
	r=l/np.linalg.norm(l)
	m=np.cross(r,v)
	n=np.cross(r,m)
	m=m/np.linalg.norm(m)
	n=n/np.linalg.norm(n)
	w=np.arange(0,2*pi,0.001)
	s=len(w)

	mm=vect_contract(m,c,m)
	mn=vect_contract(m,c,n)
	nm=vect_contract(n,c,m)
	nn0=vect_contract(n,c,n)
	nn=np.linalg.inv(nn0)
	
	val1=mm-np.dot(np.dot(mn,nn),nm)
	R=BB=np.zeros(shape=(3,3))
	for i in range(1,s):
		t=1-cos(w[i])
		CO=cos(w[i])
		SI=sin(w[i])
		R[0,0]=t*r[0]**2+CO
		R[0,1]=t*r[0]*r[1]-SI*r[2]
		R[0,2]=t*r[0]*r[2]+SI*r[1]
		R[1,0]=t*r[0]*r[1]+SI*r[2]
		R[1,1]=t*r[1]**2+CO
		R[1,2]=t*r[1]*r[2]-SI*r[0]
		R[2,0]=t*r[0]*r[2]-SI*r[1]
		R[2,1]=t*r[1]*r[2]+SI*r[0]
		R[2,2]=t*r[2]^2+CO

		mr=np.dot(R,np.transpose(m))
		nr=np.dot(R,np.transpose(n))

		mm=vect_contract(mr,c,mr)
		mn=vect_contract(mr,c,nr)
		nm=vect_contract(nr,c,mr)
		nn0=vect_contract(nr,c,nr)
		nn=np.linalg.inv(nn0)
		val2=mm-np.dot(np.dot(mn,nn),nm)
		BB=BB+0.5*(val2+val1)*(w[i]-w[i-1])

		val1=val2
	B=BB/(8*pi**2)
	return (B)

def center(A,p):
	"""
	Centers atomic positions around average point
	:param A: list of atomic points
	:param p: desired center point
	:return: atomic positions centered at p
	"""
	x1=0
	x2=0
	s=A.shape
	B=np.transpose(np.zeros(shape=s))
	if p==1:
		x1=sum(np.transpose(A)[1])/s[0]
		x2=sum(np.transpose(A)[2])/s[0]
		B[0]=A[0]
		B[1]=A[1]-x1
		B[2]=A[2]-x2
	elif p==2:
		x1=sum(np.transpose(A)[0])/s[0]
		x2=sum(np.transpose(A)[2])/s[0]
		B[1]=A[1]
		B[0]=A[0]-x1
		B[2]=A[2]-x2
	elif p==4:
		x1=sum(np.transpose(A)[0])/s[0]
		x2=sum(np.transpose(A)[1])/s[0]
		B[2]=A[2]
		B[0]=A[0]-x1
		B[1]=A[1]-x2
	B=np.transpose(B)
	return (B)

def Cijkl(C):
	"""
	Populates 4D elastic tensor from 6x6 elastic tensor
	:param C: 6x6 elastic constants tensor
	:return: 3x3x3x3 elastic Cijkl tensor
	"""
	c=np.zeros(shape=(3,3,3,3))
	CC=np.zeros(shape=(9,9))
	CC[0:6,0:6]=C[0:6,0:6]
	CC[6:9,6:9]=C[3:6,3:6]
	CC[0:6,6:9]=C[0:6,3:6]
	CC[6:9,0:6]=C[3:6,0:6]

	c[0,0,0,0]=CC[0,0]
	c[0,0,1,1]=CC[0,1]
	c[0,0,2,2]=CC[0,2]
	c[0,0,1,2]=CC[0,3]
	c[0,0,2,0]=CC[0,4]
	c[0,0,0,1]=CC[0,5]
	c[0,0,2,1]=CC[0,6]
	c[0,0,0,2]=CC[0,7]
	c[0,0,1,0]=CC[0,8]

	c[1,1,0,0]=CC[1,0]
	c[1,1,1,1]=CC[1,1]
	c[1,1,2,2]=CC[1,2]
	c[1,1,1,2]=CC[1,3]
	c[1,1,2,0]=CC[1,4]
	c[1,1,0,1]=CC[1,5]
	c[1,1,2,1]=CC[1,6]
	c[1,1,0,2]=CC[1,7]
	c[1,1,1,0]=CC[1,8]

	c[2,2,0,0]=CC[2,0]
	c[2,2,1,1]=CC[2,1]
	c[2,2,2,2]=CC[2,2]
	c[2,2,1,2]=CC[2,3]
	c[2,2,2,0]=CC[2,4]
	c[2,2,0,1]=CC[2,5]
	c[2,2,2,1]=CC[2,6]
	c[2,2,0,2]=CC[2,7]
	c[2,2,1,0]=CC[2,8]

	c[1,2,0,0]=CC[3,0]
	c[1,2,1,1]=CC[3,1]
	c[1,2,2,2]=CC[3,2]
	c[1,2,1,2]=CC[3,3]
	c[1,2,2,0]=CC[3,4]
	c[1,2,0,1]=CC[3,5]
	c[1,2,2,1]=CC[3,6]
	c[1,2,0,2]=CC[3,7]
	c[1,2,1,0]=CC[3,8]

	c[2,0,0,0]=CC[4,0]
	c[2,0,1,1]=CC[4,1]
	c[2,0,2,2]=CC[4,2]
	c[2,0,1,2]=CC[4,3]
	c[2,0,2,0]=CC[4,4]
	c[2,0,0,1]=CC[4,5]
	c[2,0,2,1]=CC[4,6]
	c[2,0,0,2]=CC[4,7]
	c[2,0,1,0]=CC[4,8]

	c[0,1,0,0]=CC[5,0]
	c[0,1,1,1]=CC[5,1]
	c[0,1,2,2]=CC[5,2]
	c[0,1,1,2]=CC[5,3]
	c[0,1,2,0]=CC[5,4]
	c[0,1,0,1]=CC[5,5]
	c[0,1,2,1]=CC[5,6]
	c[0,1,0,2]=CC[5,7]
	c[0,1,1,0]=CC[5,8]

	c[2,1,0,0]=CC[6,0]
	c[2,1,1,1]=CC[6,1]
	c[2,1,2,2]=CC[6,2]
	c[2,1,1,2]=CC[6,3]
	c[2,1,2,0]=CC[6,4]
	c[2,1,0,1]=CC[6,5]
	c[2,1,2,1]=CC[6,6]
	c[2,1,0,2]=CC[6,7]
	c[2,1,1,0]=CC[6,8]

	c[0,2,0,0]=CC[7,0]
	c[0,2,1,1]=CC[7,1]
	c[0,2,2,2]=CC[7,2]
	c[0,2,1,2]=CC[7,3]
	c[0,2,2,0]=CC[7,4]
	c[0,2,0,1]=CC[7,5]
	c[0,2,2,1]=CC[7,6]
	c[0,2,0,2]=CC[7,7]
	c[0,2,1,0]=CC[7,8]

	c[1,0,0,0]=CC[8,0]
	c[1,0,1,1]=CC[8,1]
	c[1,0,2,2]=CC[8,2]
	c[1,0,1,2]=CC[8,3]
	c[1,0,2,0]=CC[8,4]
	c[1,0,0,1]=CC[8,5]
	c[1,0,2,1]=CC[8,6]
	c[1,0,0,2]=CC[8,7]
	c[1,0,1,0]=CC[8,8]
	return(c)

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

def cubeperm(m):
	"""
	Permutes all equivalent vectors or planes in cubic crystals
	:param m: vector or plane family to permute
	:return: complete list of equivalent vectors or planes
	"""
	mm=np.zeros(shape=(24,3))
	mmm=np.zeros(shape=4)
	for i in range(4):
		mmm[0]=0
		mmm[1:3]=m
		mmm[i]=-mmm[i]
		mm[i*6,:]=(mmm[2],mmm[3],mmm[4])
		mm[i*6+1,:]=(mmm[2],mmm[4],mmm[3])
		mm[i*6+2,:]=(mmm[3],mmm[2],mmm[4])
		mm[i*6+3,:]=(mmm[3],mmm[4],mmm[2])
		mm[i*6+4,:]=(mmm[4],mmm[2],mmm[3])
		mm[i*6+5,:]=(mmm[4],mmm[3],mmm[2])
	return(mm)

def disl_maker_a(A,c,b,o=np.array([0.0,0.0,0.0]),x=np.array([1.0,0.0,0.0]),z=np.array([0.0,0.0,1.0]),l=np.array([0.0]),p=np.array([0.0])):
	"""
	Inserts a dislocation based on the anisotropic displacement solution
	:param A: atomic positions in Angstroms
	:param c: elasticity constants in Pascals
	:param b: burgers vector in Angstroms
	:param l: line direction (unit invariant)
	:param p: dislocation plane normal (unit invariant)
	:param o: dislocation core position in Angstroms
	:param x: Direction of global x axis (unit invariant)
	:param z: Direction of global z axis (unit invariant)
	:return: atomic positions including dislocation
	"""
	c_a=1.6281689374348

	sx=len(x)
	if sx==4:
		x=mbvector(x,c=c_a)
	elif sx==3:
		pass
	else:
		print 'x-direction vector is nonsensical'
		return (0)

	sz=len(z)
	if sz==4:
		z=mbvector(z,c=c_a)
	elif sz==3:
		pass
	else:
		print 'z-direction vector is nonsensical'
		return (0)

	sl=len(l)
	if sl==1:
		l=x
	elif sl==4:
		l=mbvector(l,c=c_a)
	elif sl==3:
		pass
	else:
		print 'line-direction vector is nonsensical'
		return (0)

	sp=len(p)
	if sp==1:
		p=np.array([0,0,1])
	elif sp==4:
		p=mbplane(p,c=c_a)
	elif sp==3:
		pass
	else:
		print 'dislocation plane is nonsensical'
		return (0)

	#b=b*10**-10
	#oo=o*10**-10
	#A=A*10**-10

	ex=x/np.linalg.norm(x)
	ez=z/np.linalg.norm(z)
	x0=np.array([1.0, 0.0, 0.0])
	z0=np.array([0.0, 0.0, 1.0])

	R=rotation(ex,ez,x0,z0)
	Rp=np.transpose(R)

	ep=p/np.linalg.norm(p)
	el=l/np.linalg.norm(l)

	n=ep
	m=np.cross(ep,el)
	nn=vect_contract(n,c,n)
	inn=np.linalg.inv(nn)
	nm=vect_contract(n,c,m)
	mn=np.transpose(nm)

	mm=vect_contract(m,c,m)

	N=np.zeros(shape=(6,6))
	N[0:3,0:3]=-np.dot(inn,nm)
	N[0:3,3:6]=-inn
	N[3:6,0:3]=-np.dot(np.dot(mn,inn),nm)+mm
	N[3:6,3:6]=-np.dot(mn,inn)

	V=np.zeros(shape=(6,6),dtype=np.cfloat)

	EI=np.linalg.eig(N)

	S=np.array([EI[0][0],EI[0][2],EI[0][4],EI[0][1],EI[0][3],EI[0][5]])

	V[:,0]=EI[1][:,0]
	V[:,1]=EI[1][:,2]
	V[:,2]=EI[1][:,4]
	V[:,3]=EI[1][:,1]
	V[:,4]=EI[1][:,3]
	V[:,5]=EI[1][:,5]

	#aa=np.zeros(shape=(3,6),dtype=np.cfloat)
	#ll=np.zeros(shape=(3,6),dtype=np.cfloat)
	aa=V[0:3,:]
	ll=V[3:6,:]

	"""
	for i in range(6):
		AA[:,i]=aa[:,i]/np.sqrt(2.0*np.dot(aa[:,i],ll[:,i]))
		LL[:,i]=ll[:,i]/np.sqrt(2.0*np.dot(aa[:,i],ll[:,i]))
	"""
	s=A.shape
	D=np.zeros(shape=s,dtype=np.cfloat)
	B=np.zeros(shape=s,dtype=np.cfloat)
	U=np.zeros(shape=s,dtype=np.cfloat)
	
	for i in range(s[0]):
		B[i,0]=A[i,0]-o[0]
		B[i,1]=A[i,1]-o[1]
		B[i,2]=A[i,2]-o[2]
		B[i,:]=np.dot(R,B[i,:])

	writelammps('rotated.apf',B)

	d=np.zeros(shape=(1,3))
	for i in range(s[0]):
		d=B[i,:]
		X=np.dot(m,d)
		Y=np.dot(n,d)
		ln=np.array([0.0],dtype=np.cfloat)
		for j in range(3):
			ln=np.log(X+S[j]*Y)#*np.dot(b,LL[:,j])
			U[i,:]=U[i,:]+ln*0.5*np.dot(b,ll[:,j])/np.dot(ll[:,j],aa[:,j])*aa[:,j]
		for j in range(3,6):
			ln=np.log(X+S[j]*Y)#*np.dot(b,LL[:,j])
			U[i,:]=U[i,:]-ln*0.5*np.dot(b,ll[:,j])/np.dot(ll[:,j],aa[:,j])*aa[:,j]

	U=0.5/(1j*pi)*U
	f=open('displacements.txt','w')
	for i in range(s[0]):
		f.write('%8.6f %8.6f %8.6f\n' % (U[i,0],U[i,1],U[i,2]))
	f.close()

	D=B+np.real(U)
	writelammps('with_disl.apf',D)

	for i in range(s[0]):
		D[i,:]=np.dot(Rp,D[i,:])
		D[i,0]=D[i,0]+o[0]
		D[i,1]=D[i,1]+o[1]
		D[i,2]=D[i,2]+o[2]

	return (D)

def disl_maker_i(A,b,l,oo,x,z,nu):
	"""
	Inserts a dislocation based on the isotropic displacement solution
	:param A: atomic positions before dislocation
	:param b: burgers vector of dislocaiton
	:param l: dislocation line direction
	:param oo: dislocation core location
	:param x: direction of the simulation box x-axis
	:param z: direction of the simulation box y-axis
	:param nu: Poisson's ratio
	:return: atomic positions including dislocation
	"""
	tol=10**-15
	x=x/np.linalg.norm(x)
	z=z/np.linalg.norm(z)
	ez=l/np.linalg.norm(l)
	ex=b/np.linalg.norm(b)

	s=A.shape
	D=B=np.zeros(shape=s)
	R=rotation(x,z,ex,ez)
	bs=np.dot(b,ez)
	be=abs(np.linalg.norm(np.cross(b,ez)))

	B[:,0]=A[:,0]-oo[0]
	B[:,1]=A[:,1]-oo[1]
	B[:,2]=A[:,2]-oo[2]

	B=np.dot(B,R)

	for i in range(s[0]):
		theta=atan2(B[i,1],-B[i,0])
		rad=sqrt(B[i,0]**2+B[i,1]**2)
		D[i,0]=B[i,0]+be/(2*pi)*(theta+0.5*B[i,0]*B[i,1]/(1-nu)/rad**2)
		D[i,1]=B[i,1]-be/(8*pi)/(1-nu)*((1-2*nu)*log(rad)+(B[i,0]**2-B[i,1]**2)/rad**2)
		D[i,2]=B[i,2]+bs/(2*pi)*theta
	D=np.dot(D,np.transpose(R))
	D[:,0]=D[:,0]+oo[0]
	D[:,1]=D[:,1]+oo[1]
	D[:,2]=D[:,2]+oo[2]
	return (D)

def domain_construct(basis_atoms,lx,ly,lz,rx,ry,x,zp):
	"""
	Makes an hcp lattice with an arbitrary orientation
	:param a: hcp lattice constant
	:param lx: length of the x-dimesnion of the simulation box
	:param ly: length of the y-dimesnion of the simulation box
	:param lz: length of the z-dimesnion of the simulation box
	:param x: cartesian direction of the desired a1-axis oriention
	:param y: direction of the desired
	:param zp: cartesian direction of the desired c-axis oriention
	:param c_a: c over a ratio (ideal ratio is default)
	:param p: dimensions which much be periodic
	:return: atom list of the desired simulation box
	"""
	x0=np.array([1,0,0])
	z0=np.array([0,0,1])
	b=a*sqrt(3)
	c=a*c_a

	if p==0:
		unit_cell_box=np.zeros(shape=(3,3))
	else:
		unit_cell_box=np.zeros(shape=(8,3))
		for i in range(8):
			n1=i/4
			n2=(i-n1*4)/2
			n3=i-n1*4-n2*2
			unit_cell_box[i,0]=n1*a
			unit_cell_box[i,1]=n2*b
			unit_cell_box[i,2]=n3*c

	if len(zp)==4:
		z=a*mbplane(zp,c=c_a)
	else:
		z=zp
	if len(x)==4:
		x=a*mbvector(x,c=c_a)
	else:
		pass
	if len(y)==4:
		y=a*mbvector(y,c=c_a)
	else:
		pass

	R=rotation(x0,z0,x,z)
	Rp=rotation(x,z,x0,z0)

	simulation_box=np.zeros(shape=(8,3))
	for i in range(8):
		n1=i/4
		n2=(i-n1*4)/2
		n3=i-n1*4-n2*2
		simulation_box[i,0]=(n1-0.5)*lx
		simulation_box[i,1]=(n2-0.5)*ly
		simulation_box[i,2]=(n3-0.5)*lz
	rotated_simulation_box=np.dot(simulation_box,Rp)
	sim_x=0
	sim_y=0
	sim_z=0
	for i in range(8):
		for j in range(i,8):
			dx=abs(rotated_simulation_box[i,0]-rotated_simulation_box[j,0])
			if dx>sim_x:
				sim_x=dx
			dy=abs(rotated_simulation_box[i,1]-rotated_simulation_box[j,1])
			if dy>sim_y:
				sim_y=dy
			dz=abs(rotated_simulation_box[i,2]-rotated_simulation_box[j,2])
			if dz>sim_z:
				sim_z=dz
	nx=ceil(sim_x/a)+2
	ny=ceil(sim_y/b)+2
	nz=ceil(sim_z/c)+2
	#print sim_y

	# calculate repeatable lattice dimensions in final orientation
	rotated_unit_cell=np.dot(unit_cell_box,R)
	rd=np.array([0.0,0.0,0.0])
	"""
	for i in range(8):
		for j in range(8):
			for k in range(3):
				dxx=abs(rotated_unit_cell[i,k]-rotated_unit_cell[j,k])
				if dxx>rd[k]:
					rd[k]=dxx
	"""
	rx=np.dot(x,R)
	ry=np.dot(y,R)
	rz=np.dot(z,R)

	xs=np.array([rx[0],ry[0],rz[0]])
	ys=np.array([rx[1],ry[1],rz[1]])
	zs=np.array([rx[2],ry[2],rz[2]])
	dx=np.amax(xs)
	dy=np.amax(ys)
	dz=np.amax(zs)
	if p==0:
		pass
	elif p==1:
		dx=np.dot(x,Rp)
		max_x=dx[0]
		cut_x=0.5*max_x*ceil(lx/max_x)
		cut_y=0.5*ly
		cut_z=0.5*lz
		rd[0]=max_x
		rd[1]=1e6
		rd[2]=1e6
	elif p==2:
		dy=np.dot(y,Rp)
		max_y=dy[1]
		cut_x=0.5*lx
		cut_y=0.5*max_y*ceil(ly/max_y)
		cut_z=0.5*lz
		rd[1]=max_y
	elif p==3:
		#dx=np.dot(x,Rp)
		max_x=dx
		#dy=np.dot(y,Rp)
		max_y=dy
		cut_x=0.5*max_x*ceil(lx/max_x)
		cut_y=0.5*max_y*ceil(ly/max_y)
		cut_z=0.5*lz
		rd[0]=max_x
		rd[1]=max_y
		rd[2]=1e6
		print rd
	elif p==4:
		dz=np.dot(z,Rp)
		max_z=dz[2]
		cut_x=0.5*lx
		cut_y=0.5*ly
		cut_z=0.5*max_z*ceil(lz/max_z)
		rd[2]=max_z
	elif p==5:
		dx=np.dot(x,Rp)
		max_x=dx[0]
		dz=np.dot(z,Rp)
		max_z=dz[2]
		cut_x=0.5*max_x*ceil(lx/max_x)
		cut_y=0.5*ly
		cut_z=0.5*max_z*ceil(lz/max_z)
		rd[0]=max_x
		rd[2]=max_z
	elif p==6:
		dy=np.dot(y,Rp)
		max_y=dy
		dz=np.dot(z,Rp)
		max_z=dz
		cut_x=0.5*lx
		cut_y=0.5*max_y*ceil(ly/max_y)
		cut_z=0.5*max_z*ceil(lz/max_z)
		rd[1]=max_y
		rd[2]=max_z
	elif p==7:
		dx=np.dot(x,Rp)
		max_x=dx
		dy=np.dot(y,Rp)
		max_y=dy
		dz=np.dot(z,Rp)
		max_z=dz
		cut_x=0.5*max_x*ceil(lx/max_x)
		cut_y=0.5*max_y*ceil(ly/max_y)
		cut_z=0.5*max_z*ceil(lz/max_z)
		rd[0]=max_x
		rd[1]=max_y
		rd[2]=max_z
	else:
		print 'FATAL ERROR: Periodic dimensions were defined with an unknown identifier.'
		A=0
		return (A)
	rd[1]=a*sqrt(c_a**2+3)
	N=4*nx*ny*nz
	xstart=-int(nx/2)
	ystart=-int(ny/2)
	zstart=-int(nz/2)
	if (nx%2)==1:
		xend=int(nx/2+1)
	else:
		xend=int(nx/2)
	if (ny%2)==1:
		yend=int(ny/2+1)
	else:
		yend=int(ny/2)
	if (nz%2)==1:
		zend=int(nz/2+1)
	else:
		zend=int(nz/2)

	B=np.zeros(shape=(N,3))

	aa=np.array([[0,0,0],[0.5*a,0.5*b,0],[0.5*a,b/6,0.5*c],[0,2*b/3,0.5*c]])
	n=0
	ts=time.time()
	for i in range(xstart,xend):
		for j in range(ystart,yend):
			for k in range(zstart,zend):
				B[n,0]=aa[0,0]+i*a
				B[n,1]=aa[0,1]+j*b
				B[n,2]=aa[0,2]+k*c
				B[n+1,0]=aa[1,0]+i*a
				B[n+1,1]=aa[1,1]+j*b
				B[n+1,2]=aa[1,2]+k*c
				B[n+2,0]=aa[2,0]+i*a
				B[n+2,1]=aa[2,1]+j*b
				B[n+2,2]=aa[2,2]+k*c
				B[n+3,0]=aa[3,0]+i*a
				B[n+3,1]=aa[3,1]+j*b
				B[n+3,2]=aa[3,2]+k*c
				n=n+4
				per=100.0*float(n)/float(N)
				tc=time.time()-ts
				print 'Initial lattice %8.2f%% complete. Elapsed time: %s seconds.\r' % (per, tc),
	print "Initial lattice created in %s seconds... Removing excess atoms\n" %(tc)
	print "%d atoms created before cut\n" % (N)

	Br=np.zeros(shape=(N,4))
	Br[:,0:3]=np.dot(B,R)
	num_cut=0.0
	ts=time.time()
	for i in range(int(N)):
		if Br[i,0]>(cut_x+0.05):
			Br[i,:]=np.array([0,0,0,0])
			num_cut+=1
		elif Br[i,0]<-(cut_x-0.05):
			Br[i,:]=np.array([0,0,0,0])
			num_cut+=1
		elif Br[i,1]>(cut_y+0.05):
			Br[i,:]=np.array([0,0,0,0])
			num_cut+=1
		elif Br[i,1]<-(cut_y-0.05):
			Br[i, :]=np.array([0,0,0,0])
			num_cut+=1
		elif Br[i,2]>(cut_z+0.05):
			Br[i,:]=np.array([0,0,0,0])
			num_cut+=1
		elif Br[i, 2]<-(cut_z-0.05):
			Br[i, :]=np.array([0,0,0,0])
			num_cut+=1
		Br[i,3]=sqrt(Br[i,0]**2+Br[i,1]**2+Br[i,2]**2)
		per=100.0*i/float(N)
		tf=time.time()-ts
		print 'Flagging atoms to be cut %8.2f%% complete. Elapsed time: %s seconds.\r' % (per,tf),
	print "\n%d atoms marked for removal in %s seconds.\n" % (num_cut,tf)
	I=np.argsort(Br[:,3])
	Br=Br[I,:]
	print "Removing atoms..."
	A=Br[int(num_cut):,0:3]
	if o[0]==0 and o[1]==0 and o[2]==0:
		pass
	else:
		print "Re-centering atoms."
		num_atoms=A.shape
		for i in range(num_atoms[0]):
			for j in range(3):
				A[i,j]+=o[j]

	print "HCP box construction is complete."
	box=np.zeros(shape=(3,2))
	box[0,0]=-cut_x+o[0]
	box[0,1]=cut_x+o[0]
	box[1,0]=-cut_y+o[1]
	box[1,1]=cut_y+o[1]
	box[2,0]=-cut_z+o[2]
	box[2,1]=cut_z+o[2]
	return {'A':A,'box':box,'RD':rd}

def dump_to_apf(in_name,out_name,p=7,dl=10.0):
	"""
	Converts LAMMPS dump file to loadable LAMMPS data file
	:param in_name: name of dump file to be read
	:param out_name: name of loadable data file to be created
	:param p: periodic dimensions of the box
	:return: 0
	"""
	# Column of x-coordinate in dump file
	l=np.genfromtxt(in_name,delimiter=' ',skip_header=8,max_rows=1,dtype='str')
	n=len(l)
	for i in range(n):
		if l[i]=='id':
			idc=i-2
		elif l[i]=='x':
			xc=i-2
			zc=xc+3
	
	# Conversion
	b=np.genfromtxt(in_name,delimiter='',skip_header=5,max_rows=3)
	box=np.zeros(shape=(3,2))
	box=b[:,0:2]
	if p==0:
		box[0,0]-=dl
		box[0,1]+=dl
		box[1,0]-=dl
		box[1,1]+=dl
		box[2,0]-=dl
		box[2,1]+=dl
	elif p==1:
		box[1,0]-=dl
		box[1,1]+=dl
		box[2,0]-=dl
		box[2,1]+=dl
	elif p==2:
		box[0,0]-=dl
		box[0,1]+=dl
		box[2,0]-=dl
		box[2,1]+=dl
	elif p==3:
		box[2,0]-=dl
		box[2,1]+=dl
	elif p==4:
		box[0,0]-=dl
		box[0,1]+=dl
		box[1,0]-=dl
		box[1,1]+=dl
	elif p==5:
		box[1,0]-=dl
		box[1,1]+=dl
	elif p==6:
		box[0,0]-=dl
		box[0,1]+=dl
	elif p==7:
		pass
	else:
		print 'FATAL ERROR: Periodic specification not recognized.'
		return ()
	A=np.genfromtxt(in_name,delimiter='',skip_header=9)
	writelammps(out_name,A[:,xc:zc],id=A[:,idc],p=p,box=box)
	return ()

def duplicate(A,m,box,id=np.array([0.0])):
	"""
	Duplicates atoms in periodic directions based on specified multiples along repeated dimensions
	:param A: Initial atom list
	:param m: list of multiplications in each direction 
	:param box: Repeated simulation box dimensions
	:param o: Relative box origin [-1,1] 
	:return: duplicated atomic position list
	"""
	la=A.shape
	
	lid=len(id)
	
	n=int(m[0]*m[1]*m[2])
	B=np.zeros(shape=(n*la[0],3))
	new_id=np.zeros(shape=(n*la[0],1))
	if lid==la[0] or lid==0:
		pass
	else:
		print 'WARNING: Inappropriate atom ID specification. New atom IDs will be assigned.'
	cnt=0
	dup=0
	for i in range(int(m[0])):
		print 'i is %d \n' %(i)
		dx=i*box[0]
		for j in range(int(m[1])):
			print 'j is %d \n' %(j)
			dy=j*box[1]
			for k in range(int(m[2])):
				print 'k is %d \n' %(k)
				print dup
				dz=k*box[2]
				for l in range(la[0]):
					B[cnt,0]=A[l,0]+dx
					B[cnt,1]=A[l,1]+dy
					B[cnt,2]=A[l,2]+dz
					if lid==la[0]:
						new_id[cnt,0]=id[l]+dup*la[0]
					else:
						new_id[cnt,0]=cnt+1
					cnt+=1
				dup+=1
	return {'A':B,'id':new_id}
	
def elastic_disl_energy(A,a,sfe):
	"""
	Calculates energy of a dislocation as function of distance from core
	:param A:
	:param a:
	:param sfe:
	:return:
	"""
	s=A.shape
	B=np.zeros(shape=(s[0],2))
	E=np.zeros(shape=s)

	B[:,0]=sqrt(A[:,1]**2+A[:,2]**2)
	B[:,1]=A[:,3]
	C=B[np.argsort(B[:,0])]
	Ei=0
	xmin=0
	xmax=0

	for i in range(s[0]):
		Ei=Ei+C[i,1]+a
		E[i,0]=C[i,1]
		E[i,1]=Ei
		if A[i,0]>xmax:
			xmax=A[i,0]
		elif A[i,0]<xmin:
			xmin=A[i,0]
	t=xmax-xmin
	E[:,1]=E[:,1]+sfe*t*10**-20*E[:,0]
	return (E)

def gf_l(A,f,c,potential,ld=1,B=0):
	"""
	Generates correction displacements for Green's functions solutions given atomic positions and line forces
	Calculates the Lattice GF for distances smaller than "tol" and Elastic GF for distances greater than "tol"
	Ignores effects with distances greater than "cutoff"
	:param A: atomic positions in meters
	:param f: line forces in newtons/meter
	:param c: 4th order elastic tensor
	:param ld: line direction of force (x=1, y=2, z=4)
	:param potential: file name of lgf table
	:return: displacements in meters
	"""
	if ld==1:
		m=np.array([0,1,0])
		n=np.array([0,0,1])
	elif ld==2:
		m=np.array([0,0,1])
		n=np.array([1,0,0])
	elif ld==4:
		m=np.array([1,0,0])
		n=np.array([0,1,0])
	else:
		u=0
		return (u)

	s=A.shape
	sf=f.shape
	sb=B.shape
	bb=sb[0]*sb[0]
	if bb==1:
		B=A[0:sf[0],:]

	u=np.zeros(shape=(s[0],3))

	nn=vect_contract(n,c,n)
	inn=np.linalg.inv(nn)
	nm=vect_contract(n,c,m)
	mn=np.transpose(nm)
	mm=vect_contract(m,c,m)

	N=np.zeros(shape=(6,6))
	N[0:3,0:3]=-inn*nm
	N[0:3,3:6]=-inn
	N[3:6,0:3]=-(mn)*inn*nm+mm
	N[3:6,3:6]=-mn*inn

	V=np.zeros(shape=(6,6),dtype=np.cfloat)

	EI=np.linalg.eig(N)
	S=np.array([EI[0][0],EI[0][2],EI[0][4],EI[0][1],EI[0][3],EI[0][5]])
	V[:,0]=EI[1][:,0]
	V[:,1]=EI[1][:,2]
	V[:,2]=EI[1][:,4]
	V[:,3]=EI[1][:,1]
	V[:,4]=EI[1][:,3]
	V[:,5]=EI[1][:,5]

	AA=aa=ll=np.zeros(shape=(3,6),dtype=np.cfloat)
	aa=V[0:3,:]
	ll=V[3:6,:]

	for i in range(6):
		AA[:,i]=aa[:,i]/np.sqrt(2*np.dot(aa[:,i],ll[:,i]))

	ls=np.genfromtxt(potential,max_rows=1)
	ls=int(ls)
	r0=np.genfromtxt(potential,skip_header=1,max_rows=1)
	R=np.genfromtxt(potential,skip_header=2,max_rows=1)
	lgf=np.genfromtxt(potential,delimiter="",skip_header=3,max_rows=ls)

	m=m/R
	n=n/R
	itt=0
	ts=time.time()
	for j in range(sf[0]):				# effect of force j
		for i in range(s[0]):			# on atom i
			r=A[i,:]-B[j,:]
			#t1=time.time()
			#print 't1=%s' % (t1-ts)
			rr=sqrt(r[1]**2+r[2]**2)
			a2=np.array([r[1],r[2]])
			#t2=time.time()
			#print 't2=%s' % (t2-t1)
			if rr<r0:
				u[i,:]=u[i,:]+lgf_l(a2,f[j,:],lgf)
				#t3=time.time()
				#print 't3=%s\n' % (t3-t2)
			elif rr>R:
				pass
			else:
				#t5=time.time()
				U=np.zeros(shape=(3,1))
				X=(m[0]*r[0]+m[1]*r[1]+m[2]*r[2])
				Y=(n[0]*r[0]+n[1]*r[1]+n[2]*r[2])
				ln=0
				ln=np.dtype(np.cfloat)
				for k in range(3):
					ln=np.log(X+S[k]*Y)
					U=U+AA[:,k]*ln*np.dot(f[j,:],AA[:,k])
				for k in range(4,6):
					ln=np.log(X+S[k]*Y)
					U=U-AA[:,k]*ln*np.dot(f[j,:],AA[:,k])
				U=0.5/pi*1j*U
				u[i,:]=u[i,:]+U[:,0]
				#u[i,1]=u[i,1]+U[1,0]
				#u[i,2]=u[i,2]+U[2,0]
				#t4=time.time()
				#print 'choice=%s\n' % (t5-t2)
				#print 't4=%s\n' % (t4-t2)
				#print 't5=%s\n' % (t4-t5)
			itt+=1
			per=100*itt/(sf[0]*s[0])
			tc=time.time()-ts
			if per>0:
				tr=tc/per*(100-per)
				print '%8.2f%% complete. Elapsed time: %s seconds. Estimated time remaining: %s seconds.\r' % (per,tc,tr),
			else:
				print '%8.2f%% complete. Elapsed time: %s seconds. Estimated time remaining: CALCULATING\r' % (per,tc),
	print '\n'
	return (u)

def gf_p(A,f,c):
	"""
	Calculates the displacements of atoms based on the Green's Function solutions to an array of point forces.
	:param A: list of atomic positions in meters
	:param f: forces on atoms in Newtons
	:param c: 4th order elastic tensor
	:return: displacement of atoms
	"""

def hcpmaker(a,lx,ly,lz,x,y,zp,c_a=sqrt(8/3),p=7,o=np.array([0.0,0.0,0.0])):
	"""
	Makes an hcp lattice with an arbitrary orientation
	:param a: hcp lattice constant
	:param lx: length of the x-dimesnion of the simulation box
	:param ly: length of the y-dimesnion of the simulation box
	:param lz: length of the z-dimesnion of the simulation box
	:param x: cartesian direction of the desired a1-axis oriention
	:param y: direction of the desired
	:param zp: cartesian direction of the desired c-axis oriention
	:param c_a: c over a ratio (ideal ratio is default)
	:param p: dimensions which much be periodic
	:return: atom list of the desired simulation box
	"""
	x0=np.array([1,0,0])
	z0=np.array([0,0,1])
	b=a*sqrt(3)
	c=a*c_a

	if p==0:
		unit_cell_box=np.zeros(shape=(3,3))
	else:
		unit_cell_box=np.zeros(shape=(8,3))
		for i in range(8):
			n1=i/4
			n2=(i-n1*4)/2
			n3=i-n1*4-n2*2
			unit_cell_box[i,0]=n1*a
			unit_cell_box[i,1]=n2*b
			unit_cell_box[i,2]=n3*c

	if len(zp)==4:
		z=a*mbplane(zp,c=c_a)
	else:
		z=zp
	if len(x)==4:
		x=a*mbvector(x,c=c_a)
	else:
		pass
	if len(y)==4:
		y=a*mbvector(y,c=c_a)
	else:
		pass

	R=rotation(x0,z0,x,z)
	Rp=rotation(x,z,x0,z0)

	simulation_box=np.zeros(shape=(8,3))
	for i in range(8):
		n1=i/4
		n2=(i-n1*4)/2
		n3=i-n1*4-n2*2
		simulation_box[i,0]=(n1-0.5)*lx
		simulation_box[i,1]=(n2-0.5)*ly
		simulation_box[i,2]=(n3-0.5)*lz
	rotated_simulation_box=np.dot(simulation_box,Rp)
	sim_x=0
	sim_y=0
	sim_z=0
	for i in range(8):
		for j in range(i,8):
			dx=abs(rotated_simulation_box[i,0]-rotated_simulation_box[j,0])
			if dx>sim_x:
				sim_x=dx
			dy=abs(rotated_simulation_box[i,1]-rotated_simulation_box[j,1])
			if dy>sim_y:
				sim_y=dy
			dz=abs(rotated_simulation_box[i,2]-rotated_simulation_box[j,2])
			if dz>sim_z:
				sim_z=dz
	nx=ceil(sim_x/a)+2
	ny=ceil(sim_y/b)+2
	nz=ceil(sim_z/c)+2
	#print sim_y

	# calculate repeatable lattice dimensions in final orientation
	rotated_unit_cell=np.dot(unit_cell_box,R)
	rd=np.array([0.0,0.0,0.0])
	"""
	for i in range(8):
		for j in range(8):
			for k in range(3):
				dxx=abs(rotated_unit_cell[i,k]-rotated_unit_cell[j,k])
				if dxx>rd[k]:
					rd[k]=dxx
	"""
	rx=np.dot(x,R)
	ry=np.dot(y,R)
	rz=np.dot(z,R)

	xs=np.array([rx[0],ry[0],rz[0]])
	ys=np.array([rx[1],ry[1],rz[1]])
	zs=np.array([rx[2],ry[2],rz[2]])
	dx=np.amax(xs)
	dy=np.amax(ys)
	dz=np.amax(zs)
	if p==0:
		pass
	elif p==1:
		dx=np.dot(x,Rp)
		max_x=dx[0]
		cut_x=0.5*max_x*ceil(lx/max_x)
		cut_y=0.5*ly
		cut_z=0.5*lz
		rd[0]=max_x
		rd[1]=1e6
		rd[2]=1e6
	elif p==2:
		dy=np.dot(y,Rp)
		max_y=dy[1]
		cut_x=0.5*lx
		cut_y=0.5*max_y*ceil(ly/max_y)
		cut_z=0.5*lz
		rd[1]=max_y
	elif p==3:
		#dx=np.dot(x,Rp)
		max_x=dx
		#dy=np.dot(y,Rp)
		max_y=dy
		cut_x=0.5*max_x*ceil(lx/max_x)
		cut_y=0.5*max_y*ceil(ly/max_y)
		cut_z=0.5*lz
		rd[0]=max_x
		rd[1]=max_y
		rd[2]=1e6
		print rd
	elif p==4:
		dz=np.dot(z,Rp)
		max_z=dz[2]
		cut_x=0.5*lx
		cut_y=0.5*ly
		cut_z=0.5*max_z*ceil(lz/max_z)
		rd[2]=max_z
	elif p==5:
		dx=np.dot(x,Rp)
		max_x=dx[0]
		dz=np.dot(z,Rp)
		max_z=dz[2]
		cut_x=0.5*max_x*ceil(lx/max_x)
		cut_y=0.5*ly
		cut_z=0.5*max_z*ceil(lz/max_z)
		rd[0]=max_x
		rd[2]=max_z
	elif p==6:
		dy=np.dot(y,Rp)
		max_y=dy
		dz=np.dot(z,Rp)
		max_z=dz
		cut_x=0.5*lx
		cut_y=0.5*max_y*ceil(ly/max_y)
		cut_z=0.5*max_z*ceil(lz/max_z)
		rd[1]=max_y
		rd[2]=max_z
	elif p==7:
		dx=np.dot(x,Rp)
		max_x=dx
		dy=np.dot(y,Rp)
		max_y=dy
		dz=np.dot(z,Rp)
		max_z=dz
		cut_x=0.5*max_x*ceil(lx/max_x)
		cut_y=0.5*max_y*ceil(ly/max_y)
		cut_z=0.5*max_z*ceil(lz/max_z)
		rd[0]=max_x
		rd[1]=max_y
		rd[2]=max_z
	else:
		print 'FATAL ERROR: Periodic dimensions were defined with an unknown identifier.'
		A=0
		return (A)
	rd[1]=a*sqrt(c_a**2+3)
	N=4*nx*ny*nz
	xstart=-int(nx/2)
	ystart=-int(ny/2)
	zstart=-int(nz/2)
	if (nx%2)==1:
		xend=int(nx/2+1)
	else:
		xend=int(nx/2)
	if (ny%2)==1:
		yend=int(ny/2+1)
	else:
		yend=int(ny/2)
	if (nz%2)==1:
		zend=int(nz/2+1)
	else:
		zend=int(nz/2)

	B=np.zeros(shape=(N,3))

	aa=np.array([[0,0,0],[0.5*a,0.5*b,0],[0.5*a,b/6,0.5*c],[0,2*b/3,0.5*c]])
	n=0
	ts=time.time()
	for i in range(xstart,xend):
		for j in range(ystart,yend):
			for k in range(zstart,zend):
				B[n,0]=aa[0,0]+i*a
				B[n,1]=aa[0,1]+j*b
				B[n,2]=aa[0,2]+k*c
				B[n+1,0]=aa[1,0]+i*a
				B[n+1,1]=aa[1,1]+j*b
				B[n+1,2]=aa[1,2]+k*c
				B[n+2,0]=aa[2,0]+i*a
				B[n+2,1]=aa[2,1]+j*b
				B[n+2,2]=aa[2,2]+k*c
				B[n+3,0]=aa[3,0]+i*a
				B[n+3,1]=aa[3,1]+j*b
				B[n+3,2]=aa[3,2]+k*c
				n=n+4
				per=100.0*float(n)/float(N)
				tc=time.time()-ts
				print 'Initial lattice %8.2f%% complete. Elapsed time: %s seconds.\r' % (per, tc),
	print "Initial lattice created in %s seconds... Removing excess atoms\n" %(tc)
	print "%d atoms created before cut\n" % (N)

	Br=np.zeros(shape=(N,4))
	Br[:,0:3]=np.dot(B,R)
	num_cut=0.0
	ts=time.time()
	for i in range(int(N)):
		if Br[i,0]>(cut_x+0.05):
			Br[i,:]=np.array([0,0,0,0])
			num_cut+=1
		elif Br[i,0]<-(cut_x-0.05):
			Br[i,:]=np.array([0,0,0,0])
			num_cut+=1
		elif Br[i,1]>(cut_y+0.05):
			Br[i,:]=np.array([0,0,0,0])
			num_cut+=1
		elif Br[i,1]<-(cut_y-0.05):
			Br[i, :]=np.array([0,0,0,0])
			num_cut+=1
		elif Br[i,2]>(cut_z+0.05):
			Br[i,:]=np.array([0,0,0,0])
			num_cut+=1
		elif Br[i, 2]<-(cut_z-0.05):
			Br[i, :]=np.array([0,0,0,0])
			num_cut+=1
		Br[i,3]=sqrt(Br[i,0]**2+Br[i,1]**2+Br[i,2]**2)
		per=100.0*i/float(N)
		tf=time.time()-ts
		print 'Flagging atoms to be cut %8.2f%% complete. Elapsed time: %s seconds.\r' % (per,tf),
	print "\n%d atoms marked for removal in %s seconds.\n" % (num_cut,tf)
	I=np.argsort(Br[:,3])
	Br=Br[I,:]
	print "Removing atoms..."
	A=Br[int(num_cut):,0:3]
	if o[0]==0 and o[1]==0 and o[2]==0:
		pass
	else:
		print "Re-centering atoms."
		num_atoms=A.shape
		for i in range(num_atoms[0]):
			for j in range(3):
				A[i,j]+=o[j]

	print "HCP box construction is complete."
	box=np.zeros(shape=(3,2))
	box[0,0]=-cut_x+o[0]
	box[0,1]=cut_x+o[0]
	box[1,0]=-cut_y+o[1]
	box[1,1]=cut_y+o[1]
	box[2,0]=-cut_z+o[2]
	box[2,1]=cut_z+o[2]
	return {'A':A,'box':box,'RD':rd}

def hcpperm(m):
	"""
	Permutes all equivalent vectors or planes in hexagonal crystals
	:param m: hcp vector or plane family
	:return: list of all equivalent vectors or planes
	"""
	mm=np.zeros(shape=(12,4))
	mm[0,:]=(m[0],m[1],m[2],m[3])
	mm[1,:]=(m[0],m[2],m[1],m[3])
	mm[2,:]=(m[1],m[0],m[2],m[3])
	mm[3,:]=(m[1],m[2],m[0],m[3])
	mm[4,:]=(m[2],m[1],m[0],m[3])
	mm[5,:]=(m[2],m[0],m[1],m[3])
	mm[6,:]=(m[0],m[1],m[2],-m[3])
	mm[7,:]=(m[0],m[2],m[1],-m[3])
	mm[8,:]=(m[1],m[0],m[2],-m[3])
	mm[9,:]=(m[1],m[2],m[0],-m[3])
	mm[10,:]=(m[2],m[1],m[0],-m[3])
	mm[11,:]=(m[2],m[0],m[1],-m[3])
	return (mm)

def hcptwin(a,lx,ly,lz,c=sqrt(8/3),ex=0,ey=0,pz=0,p=3,type='extension'):
	"""
	Build an hcp twin structure of given dimensions and crystallographic paramteters.
	Periodic dimensions will be rounded up to the nearest repeatable length.
	XY plane is the default twin plane.
	:param a: hcp lattice constant
	:param lx: x-length of the produced simulation box 
	:param ly: y-length of the produced simulation box
	:param lz: z-length of the produced simulation box
	:param c: c over a ratio (ideal c over a ratio is default)
	:param ex:
	:param ey:
	:param pz:
	:param p: periodic dimensions (x=1, y=2, z=4, combinations=sums)
	:return:
	"""
	if np.all(ex)==0 and np.all(ey)==0 and np.all(pz)==0:
		if type=='extension' or type=='Extension' or type=='tension' or type=='Tension' or type=='{10-12}':
			ex=np.array([2,-1,-1,0])/3.0
			ey=np.array([0,-1,1,1])
			pz=np.array([0,-1,1,2])
			gap=0.5*a
			basis=np.array([[0.5,0.1,0.0],[0.5,0.65,0.0]])
		elif type == 'contraction' or type == 'Contraction' or type == 'compression' or type == 'Compression' or type == '{10-11}':
			ex=np.array([2,-1,-1,0])/3.0
			ey=np.array([0,-1,1,2])
			pz=np.array([0,-1,1,1])
			gap=0.5*a
			basis=np.array([[0.0,0.0,0.0],[0.0,1.0/4.0,0.0],[0.0,1.0/2.0,0.0],[0.0,3.0/4.0,0.0]])
			nbsx=1.0
			nbsy=4.0
			nbsz=1.0
	elif np.all(ex)==0:
		print 'FATAL ERROR: Not enough inputs to describe the twin. 2D periodic direction missing'
		return (0)
	elif np.all(ey)==0:
		print 'FATAL ERROR: Not enough inputs to describe the twin. Twin slip direction missing'
		return (0)
	elif np.all(pz)==0:
		print 'FATAL ERROR: Not enough inputs to describe the twin. Twin plane missing'
		return(0)

	if len(ex)==4:
		Ex=mbvector(ex,c=c)
	elif len(ex)==3:
		pass
	else:
		print 'x-direction given is nonsensical'
		return (0)

	if len(ey) == 4:
		Ey=mbvector(ey, c=c)
	elif len(ey) == 3:
		pass
	else:
		print 'y-direction given is nonsensical'
		return (0)

	if len(pz)==4:
		Ez=mbplane(pz,c=c)
	elif len(pz)==3:
		pass
	else:
		print 'z-direction (plane normal) given is nonsensical'
		return (0)

	lz*=0.5

	zmir = np.array([[1.0, 0, 0], [0, 1.0, 0], [0, 0, -1.0]])
	HCP=hcpmaker(a,lx,ly,lz,ex,ey,pz,c_a=c,p=p,o=np.array([0.0,0.0,-0.5*lz]))
	bottom_half=HCP['A']
	sim_box=HCP['box']
	unit_cell=HCP['RD']
	print unit_cell
	print sim_box
	n_bottom=bottom_half.shape

	top_atom=np.amax(bottom_half[:,2])
	dz=-(gap+top_atom)

	for i in range(n_bottom[0]):
		bottom_half[i,2]+=dz

	#mirror the lattice across twin plane
	top_half=np.dot(bottom_half,zmir)

	# add missing twin-plane atoms
	twin_plane_box=np.array([(sim_box[0,1]-sim_box[0,0]),(sim_box[1,1]-sim_box[1,0]),1.0])
	print twin_plane_box
	UC=np.array([[unit_cell[0],0.0,0.0],[0.0,unit_cell[1],0.0],[0.0,0.0,unit_cell[2]]])
	basis_atoms=np.dot(basis,UC)
	#unit_cell[0]/=nbsx
	#unit_cell[1]/=nbsy
	#unit_cell[2]/=nbsz
	twin_plane_atoms=add_atoms(basis_atoms,unit_cell,twin_plane_box)
	n_twinplane=twin_plane_atoms.shape
	print n_twinplane[0]

	if p>=4:
		top_twin_plane_atoms=add_atoms(basis_atoms,unit_cell,twin_plane_box,o=np.array([0.0,0.0,lz]))
		n_tot=2*n_bottom[0]+2*n_twinplane[0]
		A=np.zeros(shape=(n_tot,3))
		A[:n_bottom[0],:]=bottom_half
		A[n_bottom[0]:2*n_bottom[0],:]=top_half
		A[2*n_bottom[0]:2*n_bottom[0]+n_twinplane[0],:]=twin_plane_atoms
		A[2*n_bottom[0]+n_twinplane[0]:,:]=top_twin_plane_atoms
	else:
		n_tot=2*n_bottom[0]+n_twinplane[0]
		A=np.zeros(shape=(n_tot,3))
		A[:n_bottom[0],:]=bottom_half
		A[n_bottom[0]:2*n_bottom[0],:]=top_half
		A[2*n_bottom[0]:,:]=twin_plane_atoms
	return {'A':A,'box':sim_box}

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

def insertCplane(A,a,lx,o=np.zeros(shape=(1,3))):
	"""
	inserts a "C"-type plane in an hcp structure.
	:param A:
	:param a:
	:param lx:
	:param o:
	:return:
	"""
	sa=A.shape
	b=a*sqrt(3)
	n=ceil(lx/b)
	print n
	B=np.zeros(shape=((sa[0]+2*n),3))
	for i in range(sa[0]):
		B[i,:]=A[i,:]-o
	for i in range(int(n)):
		B[(sa[0]+2*i+1),:]=np.array([0.0,2.0*b/3.0+float(i)*b,0.0])
		B[(sa[0]+2*i),:]=np.array([0.5*a,1.0*b/6.0+float(i)*b,0.0])
	for i in range((sa[0]+2*int(n))):
		B[i,:]=B[i,:]+o
	return (B)

def lgf_l(a,f,lgf):
	"""
	Calculates displacements due to an applied line force based on the lattice Green's Function of a
	:param a: relative atom position
	:param f: vector of line force
	:param LGF: lgf table
	:return: displacement of atom 'a' due to line force 'f' at origin
	"""
	s=lgf.shape
	lgfl=lgf*1e-10
	ds=np.zeros(shape=(s[0],3))
	ds[:,0]=lgfl[:,0]-a[0]
	ds[:,1]=lgfl[:,1]-a[1]
	for i in range(s[0]):
		ds[i,2]=sqrt(ds[i,0]**2+ds[i,1]**2)
	u=np.zeros(shape=(1,3))
	j=ds[:,2].argmin(axis=0)
	u[0,0]=lgfl[j,2]*f[0]+lgfl[j,3]*f[1]+lgfl[j,4]*f[2]
	u[0,1]=lgfl[j,5]*f[0]+lgfl[j,6]*f[1]+lgfl[j,7]*f[2]
	u[0,2]=lgfl[j,8]*f[0]+lgfl[j,9]*f[1]+lgfl[j,10]*f[2]
	return (u)

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

def norepeats(A):
	"""Removes all redundant points from a given list of points"""
	s=A.shape
	for i in range(s[0]-1):
		j=i+1
		while j<=s[0]:
			if A[i,:]==A[j,:]:
				for k in range(j,(s[0]-1)):
					A[k,:]=A[k+1,:]
				s[0]=s[0]-1
				j=j-1
			j=j+1
	B=np.zeros(shape=s)
	B=A[1:s[0],:]
	return (B)

def ogf_l(A,f,c,ld=1,R=1.0e-9):
	"""
	Generates correction displacements for Green's functions solutions given atomic positions and line forces
	Calculates the Lattice GF for distances smaller than "tol" and Elastic GF for distances greater than "tol"
	Ignores effects with distances greater than "cutoff"
	:param A: atomic positions in meters
	:param f: line forces in newtons/meter
	:param c: 4th order elastic tensor
	:param ld: line direction of force (x=1, y=2, z=4)
	:return: displacements in meters
	"""
	cutoff=R
	if ld==1:
		m=np.array([0,1,0])
		n=np.array([0,0,1])
	elif ld==2:
		m=np.array([0,0,1])
		n=np.array([1,0,0])
	elif ld==4:
		m=np.array([1,0,0])
		n=np.array([0,1,0])
	else:
		u=0
		return (u)

	s=A.shape
	u=np.zeros(shape=(s[0],3))

	nn=vect_contract(n,c,n)
	inn=np.linalg.inv(nn)
	nm=vect_contract(n,c,m)
	mn=np.transpose(nm)
	mm=vect_contract(m,c,m)

	N=np.zeros(shape=(6,6))
	N[0:3,0:3]=-inn*nm
	N[0:3,3:6]=-inn
	N[3:6,0:3]=-(mn)*inn*nm+mm
	N[3:6,3:6]=-mn*inn

	V=np.zeros(shape=(6,6),dtype=np.cfloat)

	EI=np.linalg.eig(N)
	S=np.array([EI[0][0],EI[0][2],EI[0][4],EI[0][1],EI[0][3],EI[0][5]])
	V[:,0]=EI[1][:,0]
	V[:,1]=EI[1][:,2]
	V[:,2]=EI[1][:,4]
	V[:,3]=EI[1][:,1]
	V[:,4]=EI[1][:,3]
	V[:,5]=EI[1][:,5]

	AA=aa=ll=np.zeros(shape=(3,6),dtype=np.cfloat)
	aa=V[0:3,:]
	ll=V[3:6,:]

	for i in range(6):
		AA[:,i]=aa[:,i]/np.sqrt(2*np.dot(aa[:,i],ll[:,i]))

	m=m/cutoff
	n=n/cutoff
	itt=0
	for i in range(s[0]):
		r=A[i,:]
		if ld==1:
			rr=sqrt(r[1]**2+r[2]**2)
		elif ld==2:
			rr=np.linalg.norm([r[0],r[2]])
		elif ld==4:
			rr=np.linalg.norm(r[0:2])
		if rr>cutoff:
			u[i,:]=np.zeros(shape=(1,3))
		else:
			U=np.zeros(shape=(3,1))
			X=np.dot(m,r)
			Y=np.dot(n,r)
			ln=0
			ln=np.dtype(np.complex)
			for k in range(3):
				ln=np.log(X+S[k]*Y)
				U=U+ln*AA[:,k]*ln*np.dot(f,AA[:,k])
			for k in range(4,6):
				ln=np.log(X+S[k]*Y)
				U=U-ln*AA[:,k]*ln*np.dot(f,AA[:,k])
			u[i,:]=0.5/pi*1j*U[:,0]
	return (u)

def perm(m,n):
	"""Permutes all equivalent vectors of planes calls hcpperm or cubeperm"""
	p=0
	q=0
	lm=len(m)
	ln=len(n)
	if lm!=ln:
		return
	elif lm==4:
		cm=m[0]+m[1]+m[2]
		cn=n[0]+n[1]+n[3]
		if cm!=0 or cn!=0:
			return
		mm=hcpperm(m)
		nn=hcpperm(n)
	elif lm==3:
		mm=cubeperm(m)
		nn=cubeperm(n)
	MM=norepeats(mm)
	NN=norepeats(nn)
	p=NN.shape
	q=MM.shape
	M=np.zeros(shape=(p*q,lm))
	N=np.zeros(shape=(p*q,lm))
	for i in range (p[0]):
		for j in range(q[0]):
			N[i*q+j,:]=NN[i,:]
			M[i*q+j,:]=MM[j,:]

	for ii in range(p*q):
		d=np.dot(M[ii,:],N[ii,:])
		if d!=0:
			M[ii,:]=0
			N[ii,:]=0
	C=np.zeros(shape=(p[0]*q[0],2*lm))
	C[:,0:lm]=M
	C[:,lm:2*lm]=N
	CC=norepeats(C)
	M1=CC[:,0:lm]
	N1=CC[:,lm:2*lm]
	return {'M':M1,'N':N1}

def remap(A,B,box,p=7,idA=0,idB=0):
	"""
	Remaps atom positions of 'B' which may have been remapped by LAMMPS across periodic boundary conditions.
	Useful for expanding simulation boxes of NEB simulations.
	:param A:
	:param B:
	:param box: Box bounds: xlo xhi; ylo yhi; zlo zhi
	:param p:
	:param idA:
	:param idB:
	:return:
	"""
	sA=A.shape
	sB=B.shape
	if sA[0]==sB[0]:
		if sA[1]==4:
			print 'ID info for initial configuration taken from atom list. Additional inputs ignored.'
			del idA
			idA=A[:,0]
			Ia=np.argsort(idA)
			Apos=A[Ia,1:4]
			if sB[1]==4:
				print 'ID info for final configuration taken from atom list. Additional inputs ignored.'
				del idB
				idB=B[:,0]
				Ib=np.argsort(idB)
				Bpos=B[Ib,1:4]
			elif sB[1]==3:
				if len(idB)==sB[1]:
					Ib=np.argsort(idB)
					Bpos=B[Ib,:]
				else:
					print 'WARNING: No ID info for final configuration found. Order consistency assumed.'
					Bpos=B[Ia,:]
			else:
				print 'FATAL ERROR: Insufficent number of coordinates for final configuration.'
				return (0)
		elif sA[1]==3:
			if len(idA)==sA[0]:
				Ia=np.argsort(idA)
				Apos=A[Ia,:]
				if sB[1]==4:
					del idB
					idB=B[:,0]
					Ib=np.argsort(idB)
					Bpos=B[Ib,1:4]
				elif sB[1]==3:
					if len(idB)==sB[0]:
						Ib=np.argsort(idB)
						Bpos=B[Ib,:]
					else:
						Bpos=B[Ia,:]
				else:
					print 'FATAL ERROR: Insufficient number of coordinates for final configuration.'	
				
			else:
				print 'WARNING: No initial atom IDs found. Order consistency assumed.'
				del idA
				idA=np.linspace(1,sA[0],sa[0])
		else:
			print 'FATAL ERROR: Insufficient number of coordinates for initial configuration.'
			return (0)	
					
	else:
		print 'FATAL ERROR: Configurations have different number of atoms. Aborting.'
		return (0)
		
	# Calculate box-lenths
	lx=box[0,1]-box[0,0]
	ly=box[1,1]-box[1,0]
	lz=box[2,1]-box[2,0]
		
	# Analysis
	for i in range(sA[0]):
		D=Bpos[i,:]-Apos[i,:]
		if p==1 or p==3 or p==5 or p==7:
			if D[0]>0.5*lx:
				Bpos[i,0]-=lx
			elif D[0]<-0.5*lx:
				Bpos[i,0]+=lx
			else:
				pass
		else:
			pass
		if p==2 or p==3 or p==6 or p==7:
			if D[1]>0.5*ly:
				Bpos[i,1]-=ly
			elif D[1]<-0.5*ly:
				Bpos[i,1]+=ly
			else:
				pass
		else:
			pass
		if p==4 or p==5 or p==6 or p==7:
			if D[2]>0.5*lz:
				Bpos[i,2]-=lz
			elif D[2]<-0.5*lz:
				Bpos[i,2]+=lz
			else:
				pass
		else:
			pass
		
		 
	
	return{'A':Apos,'B':Bpos}
	

#def removeplane(A,p=np.array([0,0,1]),t=0.5,d=np.array([0,1,0]),o=np.zeros(shape=(1,3))):

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
	"""
	if abs(np.dot(z1,z2))>(1-tol):
		b=0
		g=0
		a=np.arccos(np.dot(x2,x1))
		a=a.real
	else:
		n=np.cross(z1,z2)
		n=n/np.linalg.norm(n)
		a=np.arccos(np.dot(n,x1))
		b=np.arccos(np.dot(z1,z2))
		g=np.arccos(np.dot(n,x2))
		g=g.real

	R=np.zeros(shape=(3,3))
	R[0,0]=np.cos(a)*np.cos(g)-np.cos(b)*np.sin(a)*np.sin(g)
	R[0,1]=-np.cos(a)*np.sin(g)-np.cos(b)*np.cos(g)*np.sin(a)
	R[0,2]=np.sin(a)*np.sin(b)
	R[1,0]=np.cos(g)*np.sin(a)+np.cos(a)*np.cos(b)*np.sin(g)
	R[1,1]=np.cos(a)*np.cos(b)*np.cos(g)-np.sin(a)*np.sin(g)
	R[1,2]=-np.cos(a)*np.sin(b)
	R[2,0]=np.sin(b)*np.sin(g)
	R[2,1]=np.cos(g)*np.sin(b)
	R[2,2]=np.cos(b)
	"""
	R=np.zeros(shape=(3,3))
	for i in range(3):
		for j in range(3):
			R[i,j]=np.dot(e1[i,:],e2[j,:])
	R=np.transpose(R)
	return (R)

def sort_atoms(infile,outfile):
	"""
	Sorts atoms by ID number. LAMMPS outputs are usually not in order.
	:param infile: name of input file
	:param outfile: name of output file
	:return: none
	"""

def split(A,d,o=np.zeros(shape=(1,3))):
	sa=A.shape
	B=np.zeros(shape=sa)
	for i in range(sa[0]):
		B[i,:]=A[i,:]-o
		B[i,2]=B[i,2]+0.5*d*np.arctan2(B[i,2],-B[i,1])/pi
		B[i,:]=B[i,:]+o
	return (B)

def squeeze(A,d,o=np.zeros(shape=(1,3))):
	sa=A.shape
	B=np.zeros(shape=sa)
	for i in range(sa[0]):
		B[i,:]=A[i,:]-o
		B[i,2]=B[i,2]-0.5*d*np.arctan2(B[i,2],-B[i,1])/pi
		B[i,:]=B[i,:]+o
	return (B)

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

def stroh(m,n,c):
	"""
	Computes the Stroh eigenvalues and eigenvectors of the sextic formalism. Ordered such that 1-3 have positive
	imaginary components and 4-6 have negative imaginary components. 1 and 4, 2 and 5, 3 and 6 are complex conjugate
	pairs.
	:param m: unit vector orthogonal to the line direction
	:param n: unit vector mutually orthogonal to the line direction and m
	:param c: 4th order elastic tensor
	:return: 'S':Eigenvectors, 'A': Polarization vectors, 'L': Auxilliary vectors
	"""
	nn=vect_contract(n,c,n)
	inn=np.linalg.inv(nn)
	nm=vect_contract(n,c,m)
	mn=np.transpose(nm)
	mm=vect_contract(m,c,m)

	N=np.zeros(shape=(6,6))
	N[0:3,0:3]=-inn*nm
	N[0:3,3:6]=-inn
	N[3:6,0:3]=-(mn)*inn*nm+mm
	N[3:6,3:6]=-mn*inn

	V=np.zeros(shape=(6,6),dtype=np.cfloat)
	S=np.zeros(shape=(6,1),dtype=np.cfloat)
	EI=np.linalg.eig(N)
	S=np.array([EI[0][0],EI[0][2],EI[0][4],EI[0][1],EI[0][3],EI[0][5]])
	V[:,0]=EI[1][:,0]
	V[:,1]=EI[1][:,2]
	V[:,2]=EI[1][:,4]
	V[:,3]=EI[1][:,1]
	V[:,4]=EI[1][:,3]
	V[:,5]=EI[1][:,5]

	a=l=np.zeros(shape=(3,6),dtype=np.cfloat)
	A=L=a
	a=V[0:3,:]
	l=V[3:6,:]
	for i in range(6):
		A[:,i]=a[:,i]/np.sqrt(2*np.dot(a[:,i],l[:,i]))
		L[:,i]=l[:,i]/np.sqrt(2*np.dot(a[:,i],l[:,i]))
	return {'S':S,'A':A,'L':L}

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

def writelammps(fname,A,p=7,id=np.array([0]),box=np.zeros(shape=(3,2)),tri=0,type=1):
	"""
	Write atomic positions to a file which may be inputted into LAMMPS
	:param fname: filename to write positions to
	:param A: list of (rows) atomic positions
	:param p: array of periodic directions
	:param box: desired dimensions of the simulation box, if not given, dimensions will be calculated by atom positions
	:param tri: is the simulation domain orthogonal (0, default) or triclinic (1)
	:return: none
	"""
	if p==1:
		xp=1
		yp=0
		zp=0
	elif p==2:
		xp=0
		yp=1
		zp=0
	elif p==3:
		xp=1
		yp=1
		zp=0
	elif p==4:
		xp=0
		yp=0
		zp=1
	elif p==5:
		xp=1
		yp=0
		zp=1
	elif p==6:
		xp=0
		yp=1
		zp=1
	elif p==7:
		xp=1
		yp=1
		zp=1
	else:
		xp=0
		yp=0
		zp=0
	N=A.shape
	n=id.shape

	if n[0]!=N[0]:
		id=np.array(range(1,N[0]+1))
	xmin=0
	xmax=0
	ymin=0
	ymax=0
	zmin=0
	zmax=0
	for i in range(N[0]):
		if xp==1:
			if abs(A[i,0])>xmax:
				xmax=abs(A[i,0])
				xmin=-xmax
		else:
			if A[i,0]<xmin:
				xmin=A[i,0]
			elif A[i,0]>xmax:
				xmax=A[i,0]
		if yp==1:
			if abs(A[i,1])>ymax:
				ymax=abs(A[i,1])
				ymin=-ymax
		else:
			if A[i,1]<ymin:
				ymin=A[i,1]
			elif A[i,1]>ymax:
				ymax=A[i,1]
		if zp==1:
			if abs(A[i,2])>zmax:
				zmax=abs(A[i,2])
				zmin=-zmax
		else:
			if A[i,2]<zmin:
				zmin=A[i,2]
			elif A[i,2]>zmax:
				zmax=A[i,2]
	if (box[0,1]-box[0,0])!=0:
		xmax=box[0,1]
		xmin=box[0,0]
	if (box[1,1]-box[1,0])!=0:
		ymax=box[1,1]
		ymin=box[1,0]
	if (box[2,1]-box[2,0])!=0:
		zmax=box[2,1]
		zmin=box[2,0]

	f=open(fname,'w')
	f.write('LAMMPS\n')
	f.write('%i atoms\n' % (N[0]))
	f.write('1 atom types\n')
	xy = 0.0
	xz = 0.0
	yz = 0.0
	f.write('%f %f xlo xhi \n' % (xmin,xmax))
	f.write('%f %f ylo yhi \n' % (ymin,ymax))
	f.write('%f %f zlo zhi \n' % (zmin,zmax))
	if tri==1:
		f.write('%f %f %f xy xz yz \n' % (xy,xz,yz))
	else:
		pass
	f.write(' \n')
	f.write('Masses \n')

	f.write(' \n')
	f.write('1 1.0 \n')

	f.write(' \n')
	f.write('Atoms \n')
	f.write(' \n')
	
	if type==1:
		for i in range(N[0]):
			f.write('%d 1 %8.12f %8.12f %8.12f \n' % (id[i],A[i,0],A[i,1],A[i,2]))
	else:
		for i in range(N[0]):
			f.write('%d %8.12f %8.12f %8.12f \n' % (id[i],A[i,0],A[i,1],A[i,2]))
	f.close()
	return ()







