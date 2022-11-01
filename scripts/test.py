import numpy as np
import atomistic_domains as ad

x = np.array([1, 0, 0])
y = np.array([0, 1, 0])
z = np.array([0, 0, 1])
a1 = ad.Atom(id=0, position=np.array([0.5, 0.5, 0.5]))
a2 = ad.Atom(id=1, position=np.array([0.75, 0.75, 0.75]))
lat = ad.Lattice(x_vector=x, y_vector=y, z_vector=z, basis_atoms=[a1, a2])
D = ad.domain.create_from_lattice(lat, 2, 2, 2)
