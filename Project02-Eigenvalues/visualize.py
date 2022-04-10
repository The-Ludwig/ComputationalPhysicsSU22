import numpy as np
import matplotlib.pyplot as plt
from yaml import safe_load
from itertools import combinations

with open("config.yaml") as f:
    config = safe_load(f)


def get_norm(x, wavefunction):
    pdf = wavefunction**2
    dx = x[1:] - x[:-1]
    N = pdf[:-1] @ dx
    return np.sqrt(N)


for setting in config:
    data = np.genfromtxt(f"build/output/{setting['name']}.npy").T
    x, V, functions = data[0], data[1], data[2:]

    for f1, f2 in combinations(functions, 2):
        print("________________________________________")
        product = np.abs(f1 @ f2)

        if product > 1e-5:
            print(f"Warning: two functions are not orthonormal! f1*f2={product}")

    # data_reference = np.genfromtxt(f"build/output/{setting['name']}_reference.npy").T

    fig, ax1 = plt.subplots()

    ax1.set_xlabel("$z/$unitless")
    ax1.set_ylabel(r"Potential / $\hbar\omega/2$")

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax2.set_ylabel("Wavefunctions")  # we already handled the x-label with ax1

    ax1.plot(x, V, label="Potential")

    v_max = np.max(V)

    for idx, function in enumerate(functions):
        norm = get_norm(x, function)
        ax2.plot(x, function / norm, label=f"Eigenfunction #{idx+1}")

    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"build/plots/{setting['name']}.pdf")
    plt.cla()
