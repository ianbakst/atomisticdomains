from enum import Enum
from typing import Union


class Boundary(Enum):
    PERIODIC = 1
    FIXED = 2
    SHRINK_WRAP = 3
    P = 1
    F = 2
    S = 3


class BoundaryConditions:
    _x1: Boundary
    _x2: Boundary
    _x3: Boundary

    @staticmethod
    def _parse_bc(bc):
        if isinstance(bc, Boundary):
            return bc
        elif isinstance(bc, int):
            return Boundary(bc)
        elif isinstance(bc, str):
            return Boundary[bc]

    def __init__(
        self,
        x1: Union[Boundary, int, str],
        x2: Union[Boundary, int, str],
        x3: Union[Boundary, int, str],
    ):
        self._x1 = self._parse_bc(x1)
        self._x2 = self._parse_bc(x2)
        self._x3 = self._parse_bc(x3)

    def __str__(self):
        return f"<Boundary Conditions: x1 = {self._x1}, x2 = {self._x2}, x3 = {self._x3}>"

    @property
    def x1(self):
        return self._x1

    @x1.setter
    def x1(self, new_bc):
        self._x1 = self._parse_bc(new_bc)

    @property
    def x2(self):
        return self._x2

    @x2.setter
    def x2(self, new_bc):
        self._x2 = self._parse_bc(new_bc)

    @property
    def x3(self):
        return self._x3

    @x3.setter
    def x3(self, new_bc):
        self._x3 = self._parse_bc(new_bc)
