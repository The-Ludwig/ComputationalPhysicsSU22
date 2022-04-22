import numpy as np
import matplotlib.pyplot as plt
from yaml import safe_load

import seaborn as sns

sns.set_theme()

with open("config.yaml") as f:
    config = safe_load(f)


def plot_spline_test():
    data = np.genfromtxt("build/output/plot_splines.npy").T

    x, splines, dsplines, ddsplines = data[0], data[1:5], data[5:9], data[9:14]

    for spline in splines:
        plt.plot(x, spline, ".", color="blue")

    plt.tight_layout()
    plt.savefig("build/plots/spline.pdf")
    plt.cla()

    for spline in dsplines:
        plt.plot(x, spline, ".", color="blue")

    plt.tight_layout()
    plt.savefig("build/plots/dspline.pdf")
    plt.cla()

    for spline in ddsplines:
        plt.plot(x, spline, ".", color="blue")

    plt.tight_layout()
    plt.savefig("build/plots/ddspline.pdf")
    plt.cla()


plot_spline_test()


def get_knots(knot):
    if isinstance(knot, list):
        return knot
    return np.linspace(knot["start"], knot["end"], knot["N"])


for settings in config:
    name = settings["name"]
    knots = get_knots(settings["knots"])
    e = settings["e"]
    a_bohr = settings["r_bohr"]
    q = settings["q"]
    R = settings["R"]
    R1 = settings["R1"]
    R2 = settings["R2"]

    x, y = np.genfromtxt(f"build/output/solution_{name}_hydrogen.npy", unpack=True)
    plt.plot(x, y / x, label="Splines")
    plt.xlabel("$r/a_{\mathrm{Bohr}}$")
    plt.ylabel("$V(r)/\mathrm{a.u.}$")

    ana = e * (1 / (x / a_bohr) - np.exp(-2 * x / a_bohr) * (a_bohr / x + 1))
    ana = ana - ana[1] + y[1] / x[1]
    plt.plot(
        x,
        ana,
        alpha=0.5,
        label="Analytical Solution",
    )

    plt.plot(knots, [0] * len(knots), "|", label="Knots")
    plt.legend()

    plt.tight_layout()
    plt.savefig(f"build/plots/{name}_hydrogen.pdf")
    plt.cla()

    V = 4.0 * np.pi / 3.0 * (R2**3 - R1**3)

    def analytical_shell(x):
        if x < R1:
            return 4 * np.pi * q / V / 2 * (R2**2 - R1**2)
        elif x <= R2:
            return (
                4
                * np.pi
                * q
                / V
                * ((R2**2) / 2.0 - (R1**3 / x + x**2 / 2.0) / 3.0)
            )
        else:
            return q / x

    x, y = np.genfromtxt(f"build/output/solution_{name}_uniform_shell.npy", unpack=True)

    shell_np = np.vectorize(analytical_shell)
    ana = shell_np(x)
    plt.plot(x, y / x, label="Splines")

    plt.plot(
        x,
        ana - ana[1] + y[1] / x[1],
        alpha=0.5,
        label="Analytical Solution",
    )

    plt.xlabel("$r/R_2$")
    plt.ylabel("$V(r) / \mathrm{a.u.}$")

    plt.plot(knots, [0] * len(knots), "|", label="Knots")

    plt.legend()
    plt.tight_layout()
    plt.savefig(f"build/plots/{name}_shell.pdf")
    plt.cla()

    x, y = np.genfromtxt(
        f"build/output/solution_{name}_uniform_sphere.npy", unpack=True
    )
    plt.plot(x, y / x, label="Splines")

    def analytical_sphere(x):
        if x < R:
            return q / R / 2.0 * (3.0 - (x / R) ** 2)
        else:
            return q / x

    ana_np = np.vectorize(analytical_sphere)
    ana = ana_np(x)
    plt.plot(
        x,
        ana - ana[1] + y[1] / x[1],
        alpha=0.5,
        label="Analytical Solution",
    )

    plt.xlabel("$r/R$")
    plt.ylabel("$V(r) / \mathrm{a.u.}$")

    plt.plot(knots, [0] * len(knots), "|", label="Knots")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"build/plots/{name}_sphere.pdf")
    plt.cla()


# x, y = np.genfromtxt(f"build/output/{config['name']}.npy", unpack=True)

# plt.plot(x, y)
# plt.savefig(f"build/plots/{config['name']}.pdf")
