import numpy as np
import matplotlib.pyplot as plt
from yaml import safe_load

# with open("config.yaml") as f:
#     config = safe_load(f)

# x, y = np.genfromtxt(f"build/output/{config['name']}.npy", unpack=True)

# plt.plot(x, y)
# plt.savefig(f"build/plots/{config['name']}.pdf")

test_distributions = np.genfromtxt("build/output/test_normal_distribution.npy").T

print(test_distributions[0].shape)

for distr in test_distributions:
    plt.hist(distr, bins=10, density=True, histtype="step")

x = np.linspace(np.min(test_distributions), np.max(test_distributions), 1000)
plt.plot(x, 1 / np.sqrt(2 * np.pi) * np.exp(-(x**2 / 2)))
plt.savefig("build/plots/test_1D_distributions.pdf")
