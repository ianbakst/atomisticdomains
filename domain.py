import numpy as np
from math import *
import time
import os
class domain:
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
	def __init__(self,box=np.zeros(shape=(3,2)),chem=[],atoms=[],num_atoms={},x=np.zeros(shape=(1,3)),y=np.zeros(shape=(1,3)),z=np.zeros(shape=(1,3)),p=0):
		self.box = box
		self.atoms = atoms
		self.tot_atoms = len(atoms)
		self.chem = chem
		self.num_el = len(chem)
		self.num_atoms = num_atoms
		self.x = x
		self.y = y
		self.z = z
		self.p = p
		return

	def add(self,domain2):
		"""

		:param domain2:
		:return:
		"""
		return

	def create(self,lattice,lx,ly,lz,x=np.zeros(shape=(1,3)),y=np.zeros(shape=(1,3)),z=np.zeros(shape=(1,3)),p=0,origin=np.zeros(shape=(1,3))):
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
		"""
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

	def mirror(self,plane):
		"""
		Mirrors a domain about a plane and adds it to the original domain
		:param plane: plane by which the domain will be mirrored
		:return: mirrored domain
		"""

	def readlammps(self,filename,type='dump'):
		"""
		Populates domain object from LAMMPS atom position file
		:param filename: name of file to be read
		:param type: file type
		:return: none
		"""

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
