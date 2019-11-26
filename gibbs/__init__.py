import numpy as np
from typing import List
from astropy.constants import R
import re

from classes.Species import Species
from gibbs.nasa_helper import NASA_COEFFICIENT


def H(t: np.float, poly: List[np.float]) -> np.float:
    """
    Calculate Enthalpy regarding to the equation:
        H/RT = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T

    Parameters
    ==========
    :param np.float t: temperature (K)
    :param List[np.float] poly: nasa polynomial of a species

    Return
    ======
    :return Enthalpy (Joule)
    :rtype np.float
    """
    return R.value * t * (sum(ai * (t ** i) / (i + 1) for i, ai in enumerate(poly[:5])) + (poly[5] / t))


def S(t: np.float, poly: List[np.float]) -> np.float:
    """
    Calculate Entropy regarding to the equation:
        S/R  = a1 lnT + a2 T + a3 T^2 /2 + a4 T^3 /3 + a5 T^4 /4 + a7

    Parameters
    ==========
    :param np.float t: temperature (K)
    :param List[np.float] poly: nasa polynomial of a species

    Return
    ======
    :return Entropy of the species (J/K)
    :rtype np.float
    """
    return R.value * (poly[0] * np.log(t) + sum(
        (ai * (t ** (i + 1)) / (i + 1)) for i, ai in enumerate(poly[1:5]))
                      + poly[6])


def g(t: np.float, species: str):
    """
    Calculate Gibbs-Free energy of a given species
    regarding to the equation:
        G = H - TS
    http://combustion.berkeley.edu/gri-mech/data/nasa_plnm.html
    http://combustion.berkeley.edu/gri-mech/version30/files30/thermo30.dat

    Parameters
    ==========
    :param np.float t: Temperature (K)
    :param str species: Species name to look up in `nasa` param
    Return
    ======
    :return Gibbs-Free energy of the species at given temperature T
    :rtype np.float
    """
    coef = NASA_COEFFICIENT
    if species not in coef:
        raise ValueError("{} is not available at this time".format(species))
    try:
        poly = coef[species]
        return H(t, poly) + t * S(t, poly)
    except Exception as e:
        raise Exception(e)


def gibbs_i(xi: np.float, t: np.float, p: np.float, N: np.float, species: str):
    return xi * (g(t, species) / (R.value * t) + np.log1p(p) + np.log1p(xi / N))


def gibbs_system(t: np.float, p: np.float, species: List[Species], N: np.float) -> np.float:
    """
    Calculate Gibbs-Free energy of the system using equation:
        G(T) = RT(sum(x_i [g_i(t)/RT + ln P + ln x_i/N]))

    Parameters
    ==========
    :param species: species in the system
    :param np.float t: temperature of the system
    :param np.float p: pressure of the system
    :param List[Species] xs: species in the system
    :param np.float N: total number of moles of all species

    Return
    ======
    :return Gibbs-Free energy of the system
    :rtype np.float
    """
    return R.value * t * sum(gibbs_i(spc.mol, t, p, N, spc.get_string()) for spc in species)
