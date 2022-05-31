import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme()

lens, means, stds = np.genfromtxt("build/output/times_lin.npy", unpack=True)
plt.plot(lens, means, "x")


plt.tight_layout()
plt.savefig("build/plots/times_lin.pdf")
plt.cla()

lens, means, stds = np.genfromtxt("build/output/times_p2.npy", unpack=True)
plt.plot(lens, means, "x")

plt.tight_layout()
plt.savefig("build/plots/times_p2.pdf")
