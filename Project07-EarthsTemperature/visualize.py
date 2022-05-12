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
plt.cla()

ns, iters, temps = np.genfromtxt("build/output/ns.npy", unpack=True)
plt.plot(ns, iters, "x")

plt.yscale("log")
plt.xscale("log")
plt.xlabel("Number of layers $N$")
plt.ylabel("Iterations till convergence")

plt.tight_layout()
plt.savefig("build/plots/ns_iters.pdf")
plt.cla()

plt.xlabel("Number of layers $N$")
plt.ylabel("Temperature $T/\mathrm{{}^oC}$")
plt.xscale("log")

plt.plot(ns, temps, "x")
plt.savefig("build/plots/ns_temps.pdf")
plt.cla()


sigma, iters, temps = np.genfromtxt("build/output/sigma_irs.npy", unpack=True)
plt.plot(sigma, iters, "x")

plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"$\alpha_{\mathrm{IR}}$")
plt.ylabel("Iterations till convergence")

plt.tight_layout()
plt.savefig("build/plots/sigma_ir_iters.pdf")
plt.cla()

plt.xlabel(r"$\alpha_{\mathrm{IR}}$")
plt.ylabel(r"Temperature $T/\mathrm{{}^oC}$")
plt.xscale("log")

plt.plot(sigma, temps, "x")
plt.savefig("build/plots/sigma_ir_temps.pdf")
plt.cla()


sigma, iters, temps = np.genfromtxt("build/output/sigma_viss.npy", unpack=True)
plt.plot(sigma, iters, "x")

plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"$\alpha_{\mathrm{vis}}$")
plt.ylabel("Iterations till convergence")

plt.tight_layout()
plt.savefig("build/plots/sigma_vis_iters.pdf")
plt.cla()

plt.xlabel(r"$\alpha_{\mathrm{vis}}$")
plt.ylabel(r"Temperature $T/\mathrm{{}^oC}$")
plt.xscale("log")

plt.plot(sigma, temps, "x")
plt.savefig("build/plots/sigma_vis_temps.pdf")
plt.cla()
