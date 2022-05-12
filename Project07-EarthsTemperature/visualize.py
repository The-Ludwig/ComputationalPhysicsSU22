import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()

hs, rho, e, t_in_v, ts = np.genfromtxt("build/output/test.npy", unpack=True)

plt.plot(ts[:-1], hs[:-1] / 1000)
plt.xlabel(r"Temperature $T / \mathrm{{}^oC}$")
plt.ylabel(r"Height $h / \mathrm{km}$")

# plt.gca().twiny().plot(rho, hs / 1000)
# plt.xlabel(r"Density $\rho / \mathrm{kg/m^3}$")

plt.tight_layout()
plt.savefig("build/plots/temperature.pdf")
