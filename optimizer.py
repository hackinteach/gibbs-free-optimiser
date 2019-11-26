from typing import List, Set, Dict, Any, Generator, Optional, Tuple, Union, Callable

import numpy as np
from scipy.optimize import minimize

from classes.Pool import Pool
from gibbs import gibbs_system


class Optimizer:
    def __init__(self, pool: Pool, temperature: np.float, pressure: np.float):
        self.N = len(pool.species)
        self.t = temperature
        self.p = pressure
        self.species = pool.species
        self.elements: Set[str] = pool.elements
        self.pool = pool
        self.bjs = self.__initial_bj()

    def minimize(self):
        init_cond = np.array(list(map(lambda spc: spc.mol, self.species)))
        return minimize(self.__objective_function, init_cond, constraints=self.__get_constraints())

    def __objective_function(self, xs: List[np.float]):
        for i, xi in enumerate(xs):
            self.species[i].mol = xi
        return gibbs_system(self.t, self.p, self.species, self.N)

    def __initial_bj(self) -> Dict[str, np.float]:
        bs: Dict[str, int] = dict()
        for elm in self.elements:
            bs[elm] = 0
            for spc in self.species:
                bs[elm] += spc.get_element_count(elm, use_mol=True)
        return bs

    def __get_constraints(self) -> List[dict]:
        gt_zeros = [{"type": "ineq", "fun": lambda x: x[i]} for i in range(self.N)]
        mass_balance = [self.__mass_balance_constraint(bj) for bj in self.bjs.items()]
        constraints = gt_zeros + mass_balance
        print("****** CONSTRAINTS ******")
        print("Number of constraints:", len(constraints))
        print("xi > 0:", len(gt_zeros))
        print("mass balance:", len(mass_balance))
        print("*"*20)
        return constraints

    def __mass_balance_constraint(self, bjx: Tuple[str, float]) -> Dict[str, Any]:
        elm, bj = bjx
        species = self.species

        def do_balance(xs: List[np.float]) -> np.float:
            """
            Calculate the mass balance constraints
            :param xs: number of mol that is optimizing
            :return: sum(aij*xi)-bj
            :rtype: np.float
            """
            sum_ax = - bj
            for i, spc in enumerate(species):
                sum_ax += xs[i] * spc.get_element_count(elm)
            return sum_ax

        return {"type": "eq", "fun": do_balance}
