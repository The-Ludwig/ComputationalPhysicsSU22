import numpy as np
import matplotlib.pyplot as plt
from scipy.special import assoc_laguerre
import seaborn as sns
from mendeleev import element

sns.set_theme()


def integrate(x, y):
    return np.sum((y[1:] + y[:-1]) * (x[1:] - x[:-1])) / 2


def plot_rho(z, knots, ionized=False, rmax=5):
    symbol = element(z).symbol
    if ionized:
        symbol = symbol + "+"

    knots_mask = knots < rmax
    plt.plot(knots[knots_mask], [0] * knots[knots_mask].shape[0], "|", label="Knots")

    x, rho = np.genfromtxt(f"build/output/{symbol}_rho.npy", unpack=True)

    mask = x < rmax
    plt.plot(x[mask], rho[mask] * x[mask] ** 2 * np.pi * 4)

    plt.xlabel("$r/a_0$")
    plt.ylabel(r"$P(r) / a_0^{-1}$")

    plt.savefig(f"build/plots/{symbol}.pdf")

    plt.cla()


def plot_eigenfunction(
    basefilename, indices, savename, xmax=10, l=0, z=1, fakenorm=False, legend=True
):
    data = np.genfromtxt(f"build/output/{basefilename}_plot_points.npy").T
    eigvalues = np.genfromtxt(f"build/output/{basefilename}_eigenvalues.npy").T

    x, ys = data[0], data[1:]
    mask = x <= xmax

    iter = list(zip(ys, eigvalues))
    for idx in indices:
        y, eigvalue = iter[idx]
        y_ = normalize(x, y)

        ana = radial_wavefunction(x[mask], idx + 1 + l, l, z)

        if fakenorm:
            y_ = y_ / y_[mask][-10] * ana[-10]

        p = plt.plot(x[mask], y_[mask], label=f"$E_{{ {idx+1} }}={eigvalue}$")[0]

        if not legend:
            plt.text(x[mask][1] - 0.1, y_[mask][1] - 0.1, f"$E_{{ {idx+1} }}$")

        plt.plot(
            x[mask],
            ana,
            ".",
            alpha=0.5,
            color=p.get_color(),
        )

    plt.xlabel("$r/a_0$")
    plt.ylabel(r"$\psi_n(r) / a_0^{-3/2}$")

    if legend:
        plt.legend()

    plt.tight_layout()
    plt.savefig(f"build/plots/{savename}.pdf")
    plt.cla()


# plot_eigenfunction("good_wideplot", [0, 1, 2], "eigenfunctions", legend=False)
# plot_eigenfunction("good_wideplot", [10, 11, 12], "eigenfunctions_high", fakenorm=True)
# plot_eigenfunction("good_wideplot_l1", [1, 2, 3], "eigenfunctions_l1", l=1)
# plot_eigenfunction("good_wideplot_l2", [2, 3, 4], "eigenfunctions_l2", l=2)

# plot_eigenfunction(
#     "good_wideplot_l6", [6, 7, 8], "eigenfunctions_l6", fakenorm=True, l=6
# )

etas, iters = np.genfromtxt("build/output/etas.npy", unpack=True)
plt.plot(etas, iters)

plt.xlabel(r"$\eta$")
plt.ylabel(r"Iterations")

plt.tight_layout()
plt.savefig("build/plots/etas.pdf")
plt.cla()


zs, es, iterss, es_ion, iterss_ion = np.genfromtxt(
    "build/output/elements.npy", unpack=True
)

mask = np.abs(es - es_ion) < 1e6
mask = mask & (iterss < 100) & (iterss_ion < 100)

noble_gases = [2, 10, 18, 36, 54, 86]

syms = [element(int(z)).symbol for z in zs]
theo = [element(int(z)).ionenergies[1] / 27.211 for z in zs]

for z in noble_gases:
    idx = np.where(zs == z)[0][0]
    plt.annotate(syms[idx], (z, theo[idx] + 0.1))

plt.plot(zs, theo, label="Literature")

plt.plot(zs[mask], (es_ion - es)[mask], "x", label="Simulation")

plt.yscale("log")

plt.xlabel("Z")
plt.ylabel(r"$E_{\mathrm{ionization}} / \mathrm{Ht}$")

plt.legend()
plt.tight_layout()
plt.savefig("build/plots/ionization.pdf")
plt.cla()

knots = np.genfromtxt("build/output/knots_atom.npy", unpack=True)

plot_rho(10, knots)
plot_rho(10, knots, ionized=True)
plot_rho(2, knots)
plot_rho(18, knots)
plot_rho(8, knots)
plot_rho(36, knots)
