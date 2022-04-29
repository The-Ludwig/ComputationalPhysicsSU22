import numpy as np
import matplotlib.pyplot as plt
from yaml import load
from scipy.special import assoc_laguerre
import seaborn as sns

sns.set_theme()


def radial_wavefunction(x, n, l, z):
    prefac = (2.0 * z / n) ** 3 / 2.0 / n / np.prod(range(n - l, n + l + 1))
    return (
        -np.sqrt(prefac)
        * np.exp(-z * x / n)
        * (2.0 * z * x / n) ** l
        * assoc_laguerre(2.0 * z * x / n, n - l - 1, k=2.0 * l + 1.0)
    )


def integrate(x, y):
    return np.sum((y[1:] + y[:-1]) * (x[1:] - x[:-1])) / 2


def normalize(x, y):
    mask = x > 0
    prop_density = y[mask] ** 2
    N = integrate(x[mask], prop_density)
    if y[1] < 0:
        return y / x / np.sqrt(N)
    return -y / x / np.sqrt(N)


def plot_eigenfunction(
    basefilename,
    indices,
    savename,
    xmax=10,
    l=0,
    z=1,
    fakenorm=False,
    legend=True,
    analytical=True,
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

        if analytical:
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


plot_eigenfunction("good_wideplot", [0, 1, 2], "eigenfunctions", legend=False)

plot_eigenfunction("good_wideplot", [26], "free", analytical=False, xmax=100)


plot_eigenfunction("good_wideplot", [10, 11, 12], "eigenfunctions_high", fakenorm=True)
plot_eigenfunction("good_wideplot_l1", [1, 2, 3], "eigenfunctions_l1", l=1)
plot_eigenfunction("good_wideplot_l2", [2, 3, 4], "eigenfunctions_l2", l=2)
plot_eigenfunction(
    "good_wideplot_l6", [6, 7, 8], "eigenfunctions_l6", fakenorm=True, l=6
)


################
# Plot how many bound states we get in dependence of rmax
################

with open("config.yaml", "r") as f:
    configs = load(f)

rmax = []
nbound = []

for config in filter(lambda config: config["name"].startswith("rmax_"), configs):
    rmax += [float(config["knots"]["rmax"])]
    eigvalues = np.genfromtxt(f"build/output/{config['name']}_eigenvalues.npy").T

    nbound += [np.sum(eigvalues < 0)]

x = np.linspace(0, np.max(rmax), 1000)
# plt.plot(x, 3 * np.log(x) - 5)

plt.plot(rmax, nbound, "x")

# plt.yscale("log")
plt.xscale("log")

plt.xlabel(r"$r_{\mathrm{max}}/a_0$")
plt.ylabel(r"Number of found bound states")

plt.tight_layout()
plt.savefig("build/plots/rmax_v_nbound.pdf")
plt.cla()


#######################################
# plot eigenvalues
#######################################

eigvalues = np.genfromtxt(f"build/output/good_eigenvalues.npy").T
eigvalues_linear = np.genfromtxt(f"build/output/good_linear_eigenvalues.npy").T

# plt.ylim(-0.6, 0.1)
xmax = np.where(eigvalues > 0)[0][0] + 1

x = np.linspace(0.9, xmax + 1)
plt.plot(x, 0.5 / x**2, label=r"$-\frac{\mathrm{Ht}}{2N^2}$")

plt.plot(
    np.arange(eigvalues.shape[0])[eigvalues < 0] + 1,
    -eigvalues[eigvalues < 0],
    "x",
    label="Exponentially distributed knots",
)

plt.plot(
    np.arange(eigvalues_linear.shape[0])[eigvalues_linear < 0] + 1,
    -eigvalues_linear[eigvalues_linear < 0],
    "x",
    label="Linear distributed knots",
)

plt.xlim(0.9, xmax + 0.5)

plt.xlabel("N")
plt.ylabel("$-E_N/\mathrm{Hartree}$")

plt.legend()

plt.xscale("log")
plt.yscale("log")

plt.tight_layout()

plt.savefig("build/plots/eigvalues.pdf")
plt.cla()


########################
#########################
#########################
#########################
#########################
eigvalues_k7 = np.genfromtxt(f"build/output/k7_eigenvalues.npy").T

plt.plot(x, 0.5 / x**2, label=r"$-\frac{\mathrm{Ht}}{2N^2}$")

plt.plot(
    np.arange(eigvalues.shape[0])[eigvalues < 0] + 1,
    -eigvalues[eigvalues < 0],
    "x",
    label="k=4",
)


plt.plot(
    np.arange(eigvalues_k7.shape[0])[eigvalues_k7 < 0] + 1,
    -eigvalues_k7[eigvalues_k7 < 0],
    "x",
    label="k=7",
)

plt.xlim(0.9, xmax + 0.5)

plt.xlabel("N")
plt.ylabel("$-E_N/\mathrm{Hartree}$")

plt.legend()

plt.xscale("log")
plt.yscale("log")

plt.tight_layout()

plt.savefig("build/plots/ks.pdf")
plt.cla()


########################
#########################
#########################
#########################
#########################
eigvalues_few_knots_6 = np.genfromtxt(f"build/output/few_knots_6_eigenvalues.npy").T
eigvalues_few_knots_3 = np.genfromtxt(f"build/output/few_knots_3_eigenvalues.npy").T

plt.plot(x, 0.5 / x**2, label=r"$-\frac{\mathrm{Ht}}{2N^2}$")

plt.plot(
    np.arange(eigvalues_few_knots_3.shape[0])[eigvalues_few_knots_3 < 0] + 1,
    -eigvalues_few_knots_3[eigvalues_few_knots_3 < 0],
    "x",
    label="k=4",
)


plt.plot(
    np.arange(eigvalues_few_knots_6.shape[0])[eigvalues_few_knots_6 < 0] + 1,
    -eigvalues_few_knots_6[eigvalues_few_knots_6 < 0],
    "x",
    label="k=21",
)

plt.xlim(0.9, xmax + 0.5)

plt.xlabel("N")
plt.ylabel("$-E_N/\mathrm{Hartree}$")

plt.legend()

plt.xscale("log")
plt.yscale("log")

plt.tight_layout()

plt.savefig("build/plots/few_knots.pdf")
plt.cla()
