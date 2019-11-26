import re
from typing import Dict, List


def find_num_string(line: str) -> List[float]:
    """
    Helper function for parsing NASA coefficient
    """
    return list(map(float, re.findall(r'[-| ][0-9]\.[0-9]*E[-|+][0-9]{2}', line)))


def parse_nasa_coef(fn: str) -> Dict[str, List[float]]:
    """
    Parse NASA coefficient file
    """
    nasa = dict()
    curr_name = None
    with open(fn, 'r') as f:
        # skip first 5 lines
        for i in range(5):
            f.readline()
        for line in f.readlines():
            if line.strip() == 'END': break
            if line.split()[-1].strip() == '1':
                curr_name = line.split()[0].strip()
                nasa[curr_name] = []
            else:
                nasa[curr_name] += find_num_string(line)
    #     print("Species list: ", " ".join(nasa.keys()))
    return nasa


NASA_COEFFICIENT = parse_nasa_coef("data/thermo30.dat")
