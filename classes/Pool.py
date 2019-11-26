from typing import List, Set
from copy import deepcopy
from classes.Species import Species


class Pool:
    def __init__(self, species: List[Species]):
        self.species: List[Species] = deepcopy(species)
        self.elements: Set[str] = self.__get_unique_elements(species)

    @staticmethod
    def __get_unique_elements(species: List[Species]) -> Set[str]:
        return set(elm.name for spc in species for elm in spc.elements)

    def __str__(self):
        return " | ".join(str(s) for s in self.species)

    def __repr__(self):
        return str(self)
