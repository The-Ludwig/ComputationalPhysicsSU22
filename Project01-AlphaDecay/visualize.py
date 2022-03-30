import numpy as np
import matplotlib.pyplot as plt
from yaml import safe_load

with open("config.yaml") as f:
    config = safe_load(f)

for simulation in config["simulation_settings"]:
    for element in config["data"]:
        r, pot = np.genfromtxt(
            f"build/output/{simulation['Name']}-{element['Symbol']}-pot.npy",
            unpack=True,
        )

        r_pdf, pdf = np.genfromtxt(
            f"build/output/{simulation['Name']}-{element['Symbol']}-pdf.npy",
            unpack=True,
        )

        plt.step(r, pot, label="V", where="post")

        plt.xlabel("Radius r / fm")
        plt.ylabel("Potential V / MeV")
        plt.grid()

        plt.savefig(
            f"build/output/pdf/{simulation['Name']}-{element['Symbol']}-pot.pdf"
        )

        plt.plot(r_pdf, pdf * np.max(pot) / np.max(pdf), label="pdf")
        plt.ylabel("Potential V / MeV and Porbability Density / AU")
        plt.legend()

        plt.tight_layout()
        plt.savefig(
            f"build/output/pdf/{simulation['Name']}-{element['Symbol']}-pdf.pdf"
        )
        plt.cla()

    for element in config["data"]:
        r, pot = np.genfromtxt(
            f"build/output/{simulation['Name']}-{element['Symbol']}-pot.npy",
            unpack=True,
        )

        r_pdf, pdf = np.genfromtxt(
            f"build/output/{simulation['Name']}-{element['Symbol']}-pdf.npy",
            unpack=True,
        )

        plt.step(r, pot, label=element["Symbol"], where="post")

    plt.xlabel("Radius r / fm")
    plt.ylabel("Potential V / MeV")
    plt.grid()
    plt.legend()

    plt.tight_layout()
    plt.savefig(f"build/output/pdf/{simulation['Name']}-pot.pdf")
    plt.cla()
