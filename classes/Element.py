from typing import Union
import numpy as np


class Element:
    def __init__(self, name: str, count: Union[int, np.float] = 1):
        self.count = count
        self.name = name

    def __str__(self):
        num = "" if self.count == 1 else f"{self.count}"
        return f"{self.name}{num}"

    def __repr__(self):
        return str(self)
