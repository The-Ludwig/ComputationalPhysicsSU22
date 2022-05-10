import numpy as np
import matplotlib.pyplot as plt
from yaml import safe_load

with open("config.yaml") as f:
    config = safe_load(f)


x, y = np.genfromtxt(f"build/output/{config['name']}.npy", unpack=True)

plt.plot(x, y)
plt.savefig(f"build/plots/{config['name']}.pdf")
