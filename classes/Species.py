from typing import List, Union

import numpy as np

from classes.Element import Element


class Species:
    def __init__(self, elms: List[Element], mol: int = 1):
        """
        Example
        -------
        Species([Element("C"), Element("H",4)])
        """
        self.elements = elms
        self.mol = mol
        self.__name_map = {e.name: index for index, e in enumerate(elms)}

    def __str__(self):
        return str('' if self.mol == 1 else self.mol) + "".join(str(e) for e in self.elements)

    def __repr__(self):
        return str(self)

    def get_string(self):
        """
        Get string without mol
        """
        return "".join(str(e) for e in self.elements)

    def get_element_count(self, name: str, use_mol: bool = False) -> Union[int, np.float]:
        """
        Find the element count in given element `name`
        if the element is not in the species, return 0

        :param use_mol: specify whether to count `mol` or not (mol * elm_count)
        :param str name: element name to get the `count`
        :rtype int|tuple
        :return element count
        """
        if name in self.__name_map:
            index = self.__name_map[name]
            return self.elements[index].count * (self.mol if use_mol else 1)
        return 0


def update_mol(sp: Species, mol: int) -> Species:
    sp.mol = mol
    return sp
