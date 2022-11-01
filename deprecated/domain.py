import numpy as np
from math import *
import time
from typing import List
import os
from atomistic_domains.atom import Atom
from atomistic_domains.lattice import Lattice
from .elastic import elastic


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
			self,
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
					atoms.extend([atom.displace(displacement) for atom in lattice.basis_atoms])

		x_vector = nx * lattice.x_vector
		y_vector = ny * lattice.y_vector
		z_vector = nz * lattice.z_vector
		return self.create(
			atoms=atoms,
			x_vector=x_vector,
			y_vector=y_vector,
			z_vector=z_vector,
			boundary_conditions=boundary_conditions,
		)

	def create(
			self,
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


		"""
		Creates an orthorhombic, single-crystal domain based on a specified lattice, dimensions, and orientation.
		:param lattice: Repeatable lattice on which the domain will be based. (see class "lattice")
		:param lx: x-direction length of simulation domain
		:param ly: y-direction length of simulation domain
		:param lz: z-direction length of simulation domain
		:param x: crystallographic direction of the positive x-axis
		:param y: crystallographic direction of the positive y-axis
		:param z: cyrstallographic direction of the positive z-axis
		:param p: perodicity of the simulation domain
		:param origin: global coordinates of the center of the box
		:return: none
		
		self.p=p

		# Define global reference axes
		x0 = np.array([1,0,0])
		z0 = np.array([0,0,1])

		# build/define simulation box's crystallographic axes
		if np.linalg.norm(x) == 0:
			if np.linalg.norm(y)*np.linalg.norm(z) == 0:
				print ('FATAL ERROR: Not enough input directions')
				return
			else:
				x = np.cross(y,z)
				self.x = x / np.linalg.norm(x)
				self.y = y / np.linalg.norm(y)
				self.z = z / np.linalg.norm(z)
		elif np.linalg.norm(y) == 0:
			if np.linalg.norm(x)*np.linalg.norm(z) == 0:
				print ('FATAL ERROR: Not enough input directions')
				return
			else:
				y = np.cross(z,x)
				self.x = x / np.linalg.norm(x)
				self.y = y / np.linalg.norm(y)
				self.z = z / np.linalg.norm(z)
		elif np.linalg.norm(z) == 0:
			if np.linalg.norm(x) * np.linalg.norm(y) == 0:
				print ('FATAL ERROR: Not enough input directions')
				return
			else:
				z = np.cross(x, y)
				self.x = x / np.linalg.norm(x)
				self.y = y / np.linalg.norm(y)
				self.z = z / np.linalg.norm(z)

		# Rotation Matrices
		R = rotation(x0, z0, self.x, self.z)
		Rp = rotation(self.x, self.z, x0, z0)

		# Determine minimum lattice repetition for rotated box
		simulation_box = np.zeros(shape = (8,3))
		for i in range(8):
			n1 = i / 4
			n2 = (i - n1 * 4) / 2
			n3 = i - n1 * 4 - n2 * 2
			simulation_box[i,0] = (n1 - 0.5) * lx
			simulation_box[i,1] = (n2 - 0.5) * ly
			simulation_box[i,2] = (n3 - 0.5) * lz

		rotated_simulation_box=np.dot(simulation_box,Rp)

		sim_x = 0
		sim_y = 0
		sim_z = 0
		for i in range(8):
			for j in range(i,8):
				dx = abs(rotated_simulation_box[i, 0] - rotated_simulation_box[j, 0])
				if dx > sim_x:
					sim_x = dx
				dy = abs(rotated_simulation_box[i, 1] - rotated_simulation_box[j, 1])
				if dy > sim_y:
					sim_y = dy
				dz = abs(rotated_simulation_box[i, 2] - rotated_simulation_box[j, 2])
				if dz > sim_z:
					sim_z = dz

		# Repeated lattice multiples to be created before rotation
		nx = ceil(sim_x / lattice.box[0,0]) + 2
		ny = ceil(sim_y / lattice.box[1,1]) + 2
		nz = ceil(sim_z / lattice.box[2,2]) + 2

		# Calculate rotated box-cut dimensions for desired periodic conditions
		if p == 1:
			rx=lattice.repeated_length(self.x)
			cut_x = 0.5 * rx * ceil(lx / rx)
			cut_y = 0.5 * ly
			cut_z = 0.5 * lz
		elif p == 2:
			ry=lattice.repeated_length(self.y)
			cut_x = 0.5 * lx
			cut_y = 0.5 * ry * ceil(ly / ry)
			cut_z = 0.5 * lz
		elif p == 3:
			rx=lattice.repeated_length(self.x)
			ry=lattice.repeated_length(self.y)
			cut_x = 0.5 * rx * ceil(lx / rx)
			cut_y = 0.5 * ry * ceil(ly / ry)
			cut_z = 0.5 * lz
		elif p == 4:
			rz=lattice.repeated_length(self.z)
			cut_x = 0.5 * lx
			cut_y = 0.5 * ly
			cut_z = 0.5 * rz * ceil(lz / rz)
		elif p == 5:
			rx=lattice.repeated_length(self.x)
			rz=lattice.repeated_length(self.z)
			cut_x = 0.5 * rx * ceil(lx / rx)
			cut_y = 0.5 * ly
			cut_z = 0.5 * rz * ceil(lz / rz)
		elif p == 6:
			ry = lattice.repeated_length(self.y)
			rz = lattice.repeated_length(self.z)
			cut_x = 0.5 * lx
			cut_y = 0.5 * ry * ceil(ly / ry)
			cut_z = 0.5 * rz * ceil(lz / rz)
		elif p == 7:
			rx = lattice.repeated_length(self.x)
			ry = lattice.repeated_length(self.y)
			rz = lattice.repeated_length(self.z)
			cut_x = 0.5 * rx * ceil(lx / rx)
			cut_y = 0.5 * ry * ceil(ly / ry)
			cut_z = 0.5 * rz * ceil(lz / rz)
		else:
			print ('WARNING: Periodic directions were either undefined or, defined with an unknown identifier. No periodicity will be enforced.')
			cut_x = 0.5 * lx
			cut_y = 0.5 * ly
			cut_z = 0.5 * lz

		# Lattice constants extracted from box
		a=np.linalg.norm(np.dot(lattice.box,np.array([[1],[0],[0]])))
		b=np.linalg.norm(np.dot(lattice.box,np.array([[0],[1],[0]])))
		c=np.linalg.norm(np.dot(lattice.box,np.array([[0],[0],[1]])))

		N = lattice.tot_atoms * nx * ny * nz

		print (nx)
		print (ny)
		print (nz)
		print (lattice.tot_atoms)

		xstart = -int(nx / 2)
		ystart = -int(ny / 2)
		zstart = -int(nz / 2)
		if (nx % 2) == 1:
			xend = int(nx / 2 + 1)
		else:
			xend = int(nx / 2)
		if (ny % 2) == 1:
			yend = int(ny / 2 + 1)
		else:
			yend = int(ny / 2)
		if (nz % 2) == 1:
			zend = int(nz / 2 + 1)
		else:
			zend = int(nz / 2)

		#aa = np.array([[0, 0, 0], [0.5 * a, 0.5 * b, 0], [0.5 * a, b / 6, 0.5 * c], [0, 2 * b / 3, 0.5 * c]])
		n = 0
		self.atoms=[]
		self.chem=lattice.chem
		self.num_el=len(self.chem)
		for i in range(self.num_el):
			self.num_atoms[self.chem[i]] = 0
		list_of_atoms=[]
		ts = time.time()
		tc = time.time() - ts
		for i in range(xstart, xend):
			for j in range(ystart, yend):
				for k in range(zstart, zend):
					for l in range(lattice.tot_atoms):
						n += 1
						temp_atom = copy.deepcopy(lattice.basis_atoms[l])
						temp_atom.scale(lattice.box)
						disp = np.array([i*a, j*b, k*c])
						temp_atom.displace(disp)
						temp_atom.scale(R)
						list_of_atoms.append(temp_atom)
						del temp_atom
						per = 100.0 * float(n) / float(N)
						tc = time.time() - ts
						print ('Initial lattice %8.2f%% complete. Elapsed time: %s seconds.\r' % (per, tc),)
		print ("Initial lattice created in %s seconds... Removing excess atoms\n" % (tc))
		print ("%d atoms created before cut\n" % (N))

		atomid = 0
		ts = time.time()
		tf = time.time() - ts
		for i in range(len(list_of_atoms)):
			if list_of_atoms[i].pos[0] > (cut_x + 0.05):
				pass
			elif list_of_atoms[i].pos[0] < -(cut_x - 0.05):
				pass
			elif list_of_atoms[i].pos[1] > (cut_y + 0.05):
				pass
			elif list_of_atoms[i].pos[1] < -(cut_y - 0.05):
				pass
			elif list_of_atoms[i].pos[2] > (cut_z + 0.05):
				pass
			elif list_of_atoms[i].pos[2] < -(cut_z - 0.05):
				pass
			else:
				atomid += 1
				list_of_atoms[i].change_id(atomid)
				list_of_atoms[i].displace(origin)
				self.atoms.append(list_of_atoms[i])
				self.num_atoms[list_of_atoms[i].element] += 1
			per = 100.0 * i / float(N)
			tf = time.time() - ts
			print ('Cutting atoms: %8.2f%% complete. Elapsed time: %s seconds.\r' % (per, tf),)
		print ("\n%d atoms remain in %s seconds.\n" % (atomid, tf))
		self.tot_atoms = atomid
		print ("HCP box construction is complete.")
		self.box[0, 0] = -cut_x + origin[0]
		self.box[0, 1] = cut_x + origin[0]
		self.box[1, 0] = -cut_y + origin[1]
		self.box[1, 1] = cut_y + origin[1]
		self.box[2, 0] = -cut_z + origin[2]
		self.box[2, 1] = cut_z + origin[2]
		return
	"""
	def mirror(self,plane):
		"""
		Mirrors a domain about a plane and adds it to the original domain
		:param plane: plane by which the domain will be mirrored
		:return: mirrored domain
		"""
		return

	def readlammps(self,filename,type='dump'):
		"""
		Populates domain object from LAMMPS atom position file
		:param filename: name of file to be read
		:param type: file type
		:return: none
		"""
		return

	def shift(self, cut=0.0, dx=0.0, dy=0.0, dz=0.0, plane='z'):
		"""
		Shifts atoms above a specified plane by specified displacements. Useful for generating GSF curves.
		:param cut: Height of plane above which atoms will be displaced
		:param dx: Distance atoms will be moved in x-direction
		:param dy: Distance atoms will be moved in y-direction
		:param dz: Distance atoms will be moved in z-direction
		:param plane: Plane above which atoms will be displaced
		:return: New vasp object with updated coordinates.
		"""
		new = copy.deepcopy(self)
		if plane == 'x' or plane == 'X':
			for i in range(self.tot_atoms):
				if self.atoms[i].pos[0] >= cut:
					new.atoms[i].displace(np.array([0.0,dy,dz]))
				else:
					pass
		elif plane == 'y' or plane == 'Y':
			for i in range(self.tot_atoms):
				if self.atoms[i].pos[1] >= cut:
					new.atoms[i].displace(np.array([dx,0.0,dz]))
				else:
					pass
		elif plane == 'z' or plane == 'Z':
			for i in range(self.tot_atoms):
				if self.atoms[i].pos[2] >= cut:
					new.atoms[i].displace(np.array([dx,dy,0.0]))
				else:
					pass
		return (new)

	def subtract(self,shape):
		"""
		Removes atoms from domain where the domain overlaps shape.
		:param shape: Shape within which atoms will be removed.
		:return: NONE
		"""
		return

	def writeposcar(self, filename, comment='Structure', ptype='Direct'):
		"""
		Writes domain object to POSCAR
		:param filename: name of POSCAR to be written
		:param comment: Opening line comment
		:param ptype: POSCAR type, based on dynamics
		:return: NONE
		"""
		a = self.box[0,1] - self.box[0,0]
		b = self.box[1,1] - self.box[1,0]
		c = self.box[2,1] - self.box[2,0]
		f = open(filename, 'w')
		f.write(comment + '\n')
		f.write('1.000000 \n')
		f.write('%8.12f    %8.12f    %8.12f\n' % (a, 0.0, 0.0))
		f.write('%8.12f    %8.12f    %8.12f\n' % (0.0, b, 0.0))
		f.write('%8.12f    %8.12f    %8.12f\n' % (0.0, 0.0, c))

		for i in range(self.num_el):
			f.write('%s   ' % (self.chem[i])),
		f.write('\n')
		for i in range(self.num_el):
			f.write('%d   ' % (self.num_atoms[self.chem[i]])),
		f.write('\n')
		f.write('%s\n' % (ptype))
		if ptype == 'Direct':
			for i in range(self.num_el):
				for j in range(self.tot_atoms):
					if self.atoms[j].element == self.chem[i]:
						x = self.atoms[j].pos[0] / float(a)
						y = self.atoms[j].pos[1] / float(b)
						z = self.atoms[j].pos[2] / float(c)
						f.write('%8.12f %8.12f %8.12f\n' % (x, y, z))
		elif ptype == 'Selective Dynamics':
			f.write('Cartesian\n')
			for i in range(self.num_el):
				for j in range(self.tot_atoms):
					if self.atoms[j].element == self.chem[i]:
						x = self.atoms[j].pos[0]
						y = self.atoms[j].pos[1]
						z = self.atoms[j].pos[2]
						f.write('%8.12f %8.12f %8.12f %s %s %s\n' % (x, y, z, self.atoms[j].dyn[0], self.atoms[j].dyn[1], self.atoms[j].dyn[2]))
		else:
			pass
		f.write('\n')
		for i in range(self.tot_atoms):
			f.write('%8.12f %8.12f %8.12f\n' % (self.atoms[i].mag_mom[0], self.atoms[i].mag_mom[1], self.atoms[i].mag_mom[2]))
		return ('POSCAR written to file: ' + filename)

	def writelammps(self, filename, type='apf',tri=0):
		"""
		Writes domain to a LAMMPS input file
		:param filename: name of file to be written
		:param type: type of LAMMPS input file to be written
		:return: none
		"""
		f = open(filename, 'w')
		f.write('LAMMPS\n')
		f.write('%i atoms\n' % (self.tot_atoms))
		f.write('1 atom types\n')
		f.write('%f %f xlo xhi \n' % (self.box[0,0], self.box[0,1]))
		f.write('%f %f ylo yhi \n' % (self.box[1,0], self.box[1,1]))
		f.write('%f %f zlo zhi \n' % (self.box[2,0], self.box[2,1]))
		if tri == 1:
			f.write('0.0 0.0 0.0 xy xz yz \n')
		else:
			pass
		f.write(' \n')
		f.write('Masses \n')

		f.write(' \n')
		f.write('1 1.0 \n')

		f.write(' \n')
		f.write('Atoms \n')
		f.write(' \n')

		if type == 1:
			for i in range(self.tot_atoms):
				f.write('%d 1 %8.12f %8.12f %8.12f \n' % (self.atoms[i].id, self.atoms[i].pos[0], self.atoms[i].pos[1], self.atoms[i].pos[2]))
		else:
			for i in range(self.tot_atoms):
				f.write('%d %8.12f %8.12f %8.12f \n' % (self.atoms[i].id, self.atoms[i].pos[0], self.atoms[i].pos[1], self.atoms[i].pos[2]))
		f.close()
		return

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
			print('x-direction vector is nonsensical')
			return (0)

		sz=len(z)
		if sz==4:
			z=mbvector(z,c=c_a)
		elif sz==3:
			pass
		else:
			print('z-direction vector is nonsensical')
			return (0)

		sl=len(l)
		if sl==1:
			l=x
		elif sl==4:
			l=mbvector(l,c=c_a)
		elif sl==3:
			pass
		else:
			print('line-direction vector is nonsensical')
			return (0)

		sp=len(p)
		if sp==1:
			p=np.array([0,0,1])
		elif sp==4:
			p=mbplane(p,c=c_a)
		elif sp==3:
			pass
		else:
			print('dislocation plane is nonsensical')
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
			print('WARNING: Inappropriate atom ID specification. New atom IDs will be assigned.')
		cnt=0
		dup=0
		for i in range(int(m[0])):
			print('i is %d \n' %(i))
			dx=i*box[0]
			for j in range(int(m[1])):
				print('j is %d \n' %(j))
				dy=j*box[1]
				for k in range(int(m[2])):
					print('k is %d \n' %(k))
					print(dup)
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

	def remap(A,B,box,p=7,idA=0,idB=0):
		
		Remaps atom positions of 'B' which may have been remapped by LAMMPS across periodic boundary conditions.
		Useful for expanding simulation boxes of NEB simulations.
		:param A:
		:param B:
		:param box: Box bounds: xlo xhi; ylo yhi; zlo zhi
		:param p:
		:param idA:
		:param idB:
		:return:
		
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
"""