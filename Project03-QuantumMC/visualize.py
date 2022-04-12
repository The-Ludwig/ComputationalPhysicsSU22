import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from yaml import safe_load

sns.set_theme()

#############################
# 1D Distribution
#############################

test_distributions = np.genfromtxt("build/output/test_normal_distribution.npy").T

plt.hist(
    test_distributions[0], bins=30, density=True, histtype="step", label="Uniform Step"
)
plt.hist(
    test_distributions[1], bins=30, density=True, histtype="step", label="Gaussian Step"
)

x = np.linspace(np.min(test_distributions), np.max(test_distributions), 1000)
plt.plot(x, 1 / np.sqrt(2 * np.pi) * np.exp(-(x**2 / 2)), label="Expected")

plt.xlabel("x")
plt.ylabel("P(x)")

plt.legend()
plt.tight_layout()
plt.savefig("build/plots/test_1D_distributions.pdf")
plt.cla()


def plot_energy(file_input, file_output):

    alpha, energy, std = np.genfromtxt(file_input).T

    # sns.lineplot(
    #     x=alpha,
    #     y=energy,
    #     hue=std,
    #     # style="event",
    # )
    plt.plot(alpha, energy)
    plt.fill_between(
        alpha, (energy - 0.5 * std), (energy + 0.5 * std), color="b", alpha=0.1
    )

    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$\langle H \rangle / \hbar\omega$")

    plt.tight_layout()
    plt.savefig(file_output)
    plt.cla()


plot_energy("build/output/test_1d_energy.npy", "build/plots/1d_energy.pdf")

##############################
# 2D Distribution
##############################
x1, y1, x2, y2 = np.genfromtxt("build/output/test_2D_distribution.npy").T

N = 100
# plt.scatter(x1[:N], y1[:N])
# plt.scatter(x2[:N], y2[:N])

sns.kdeplot(
    x=x1[:N],
    y=y1[:N],
    fill=True,
    clip=(-3, 3),
    thresh=0,
    levels=50,
    cmap="rocket",
)

plt.tight_layout()
plt.savefig("build/plots/density_first_particle.pdf")
plt.cla()


sns.kdeplot(
    x=x2[:N],
    y=y2[:N],
    fill=True,
    clip=(-3, 3),
    thresh=0,
    levels=50,
    cmap="rocket",
)

plt.tight_layout()
plt.savefig("build/plots/density_second_particle.pdf")
plt.cla()


data_1 = pd.DataFrame(
    data={"x": x1[:N], "y": y1[:N], "particle": np.zeros(N, dtype=int)}
)
data_2 = pd.DataFrame(
    data={"x": x2[:N], "y": y2[:N], "particle": np.ones(N, dtype=int)}
)
data = pd.concat([data_1, data_2], ignore_index=True)
print(data)

# Show the joint distribution using kernel density estimation
g = sns.jointplot(
    data=data,
    x="x",
    y="y",
    hue="particle",
    kind="kde",
)

plt.tight_layout()
plt.savefig("build/plots/density_combined.pdf")
plt.clf()

plot_energy("build/output/test_2d_energy_0.npy", "build/plots/2d_energy_0.pdf")
plot_energy("build/output/test_2d_energy_1.npy", "build/plots/2d_energy_1.pdf")
plot_energy("build/output/test_2d_energy_1.npy", "build/plots/2d_energy_2.pdf")
plot_energy("build/output/test_2d_energy_8.npy", "build/plots/2d_energy_8.pdf")

r_particles = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
r1 = np.sqrt(x1**2 + y1**2)
r2 = np.sqrt(x2**2 + y2**2)

plt.hist(r_particles, bins=50, density=True)
plt.xlabel("$r_{ij}$")
plt.ylabel("$P(r_{ij})$")

plt.tight_layout()
plt.savefig("build/plots/r_distance.pdf")
plt.clf()


plt.hist(r1, bins=50, histtype="step", label="1", density=True)
plt.hist(r2, bins=50, histtype="step", label="2", density=True)

plt.xlabel("$r$")
plt.ylabel("$P(r)$")

plt.legend()
plt.tight_layout()
plt.savefig("build/plots/r_particles.pdf")
plt.clf()


sns.jointplot(x=r_particles, y=r1, kind="hex", color="#4CB391")
plt.xlabel("$r_{ij}$")
plt.ylabel("$r_1$")

plt.legend()
plt.tight_layout()
plt.savefig("build/plots/hex_r1.pdf")
plt.clf()


sns.jointplot(x=r_particles, y=r2, kind="hex", color="#4CB391")
plt.xlabel("$r_{ij}$")
plt.ylabel("$r_2$")

plt.legend()
plt.tight_layout()
plt.savefig("build/plots/hex_r2.pdf")
plt.clf()
