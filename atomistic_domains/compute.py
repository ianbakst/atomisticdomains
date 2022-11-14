from abc import ABC, abstractmethod


class Compute(ABC):
    _executable: str

    @abstractmethod
    def setup(self):
        pass

    @abstractmethod
    def run(self):
        pass

    @property
    def executable(self):
        return self._executable


class Vasp(Compute):
    def __init__(self):
        self._executable =

    def setup(self):
        pass

    def run(self):
