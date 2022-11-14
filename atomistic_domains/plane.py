import numpy as np


class Plane:
    _normal_vector: np.ndarray
    _point: np.ndarray

    def __init__(self, normal_vector: np.ndarray, point: np.ndarray = np.zeros(shape=(1, 3))):
        self._normal_vector = normal_vector
        self._point = point

    @property
    def normal_vector(self):
        return self._normal_vector

    @normal_vector.setter
    def normal_vector(self, normal_vector: np.ndarray):
        self._normal_vector = normal_vector
