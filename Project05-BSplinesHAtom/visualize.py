import numpy as np
import matplotlib.pyplot as plt


data = np.genfromtxt("build/output/good_plot_points.npy").T

x, y = data[0], data[1]
mask = x <= 1

plt.plot(x[mask], y[mask])
plt.tight_layout()
plt.savefig("build/plots/good_plot.pdf")
