from classes.Element import Element
from classes.Pool import Pool
from classes.Species import Species
from optimizer import Optimizer

m0_H2O2 = Species([Element("H", 2), Element("O", 2)], mol=0)
m0_H2O = Species([Element("H", 2), Element("O", 1)], mol=0)
m0_CO2 = Species([Element("C"), Element("O", 2)], mol=0)
m0_OH = Species([Element("O"), Element("H")], mol=0)
m0_CO = Species([Element("C"), Element("O")], mol=0)
CH4 = Species([Element("C"), Element("H", 4)])
O2 = Species([Element("O", 2)])
H2 = Species([Element("H", 2)])

pool_1 = Pool([O2, H2, m0_OH, m0_H2O, m0_H2O2, CH4, m0_CO, m0_CO2])
pool_2 = Pool([O2, H2, m0_H2O])
pool_3 = Pool([O2, m0_H2O, CH4, m0_CO2])

pools = [pool_1, pool_2, pool_3]
pres = 100_000
t = 500  # Kelvin
if __name__ == '__main__':
    for p in pools:
        print("Processing:", p)
        proj = Optimizer(p, t, pres)
        result = proj.minimize()
        print(result)
        print("\n" + "=" * 50 + "\n")
