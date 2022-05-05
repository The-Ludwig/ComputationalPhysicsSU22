from mendeleev import element
from yaml import dump

elements = []


def get_l(orb):
    match orb:
        case "s":
            return 0
        case "p":
            return 1
        case "d":
            return 2
        case "f":
            return 3
        case _:
            raise ValueError("orbital has to be s, p, d, f")


for i in range(2, 100):
    e = element(i)
    el = {"name": e.name, "symbol": e.symbol, "Z": i}
    orbitals = [[value, key[0], get_l(key[1])] for key, value in e.ec.conf.items()]
    orbitals_ionized = [
        [value, key[0], get_l(key[1])] for key, value in e.ec.ionize(1).conf.items()
    ]
    el["orbitals"] = orbitals
    el["orbitals_ionized"] = orbitals_ionized

    elements.append(el)


with open("build/output/pse.yaml", "w") as f:
    dump(elements, f)
