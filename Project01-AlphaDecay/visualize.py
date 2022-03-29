import numpy as np
import matplotlib.pyplot as plt
from yaml import safe_load

with open("config.yaml") as f:
    config = safe_load(f)

for simulation in config["simulation_settings"]:
    for element in config["data"]:
        r, pot = np.genfromtxt(
            f"build/output/{simulation['Name']}-{element['Symbol']}-pot.npy",
            unpack=True,
        )

        plt.plot(r, pot)

        plt.xlabel("Radius r / fm")
        plt.ylabel("Potential V / MeV")
        plt.grid()

        plt.tight_layout()
        plt.savefig(f"build/{simulation['Name']}-{element['Symbol']}-pot.pdf")
